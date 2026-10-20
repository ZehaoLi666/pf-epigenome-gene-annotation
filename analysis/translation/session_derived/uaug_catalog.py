"""Session-derived computation of uaug_catalog.tsv.gz.

Recovered from the record of the interactive analysis session in which uaug_catalog.tsv.gz was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: gene_boundaries_v3_all_genes.tsv, Lab/PlasmoDB/PlasmoDB-68/PlasmoDB-68_Pfalciparum3D7.gff, refined_utr_annotation.gff3, Lab/IGV_longread_transcript_models/PlasmoDB-68_Pfalciparum3D7_Genome.fasta
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import LAB_DIR, result  # paths: config/paths.yaml
import pandas as pd
import numpy as np
import re
import warnings
warnings.filterwarnings('ignore')

# Load gene boundaries
gb = pd.read_csv(result('gene_boundaries_v3_all_genes.tsv'), sep='\t')

# Load genome
def read_fa(p):
    seqs = {}; name = None; buf = []
    for l in open(p):
        if l.startswith('>'):
            if name: seqs[name] = ''.join(buf)
            name = l[1:].split()[0].split('|')[0].strip(); buf = []
        else: buf.append(l.strip().upper())
    seqs[name] = ''.join(buf); return seqs

genome = read_fa(f'{LAB_DIR}/IGV_longread_transcript_models/PlasmoDB-68_Pfalciparum3D7_Genome.fasta')
comp = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(comp)[::-1]

# Parse PlasmoDB GFF
G = f'{LAB_DIR}/PlasmoDB/PlasmoDB-68/PlasmoDB-68_Pfalciparum3D7.gff'
rows = []
for l in open(G):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if len(f) < 9: continue
    if f[2] in ('CDS', 'five_prime_UTR', 'three_prime_UTR', 'exon', 'mRNA'):
        a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
        rows.append((f[0], f[2], int(f[3]), int(f[4]), f[6], a.get('Parent', ''), a.get('ID', '')))
g68 = pd.DataFrame(rows, columns=['chrom', 'type', 'start', 'end', 'strand', 'parent', 'id'])

tx2gene = {}
for l in open(G):
    if '\tmRNA\t' in l:
        a = dict(x.split('=', 1) for x in l.rstrip().split('\t')[8].split(';') if '=' in x)
        tx2gene[a['ID']] = a['Parent']
g68['gene'] = g68.parent.map(tx2gene)
g68.loc[g68.type == 'mRNA', 'gene'] = g68.loc[g68.type == 'mRNA', 'parent']

genes = set(gb.gene_id)
gbi = gb.set_index('gene_id')

# Select transcript per gene
c1 = g68[(g68.type == 'CDS') & g68.gene.isin(genes) & g68.parent.str.endswith('.1')]
g68g = g68[g68.gene.isin(genes)]
cdsx = g68g[g68g.type == 'CDS'].groupby('parent').agg(gene=('gene', 'first'), s=('start', 'min'), e=('end', 'max'), strand=('strand', 'first'))
cdsx = cdsx.join(gbi[['atg', 'stop']], on='gene')
cdsx['ok'] = np.where(cdsx.strand == '+', (cdsx.s - 1 == cdsx.atg) & (cdsx.e - 1 == cdsx.stop), (cdsx.e - 1 == cdsx.atg) & (cdsx.s - 1 == cdsx.stop))
cdsx['is1'] = cdsx.index.str.endswith('.1')
sel = cdsx.reset_index().sort_values(['gene', 'ok', 'is1'], ascending=[True, False, False]).drop_duplicates('gene').set_index('gene')
tx_of = sel.parent.to_dict()
feat = g68g[g68g.parent.isin(set(tx_of.values()))]

# CDS intervals per gene
cds_sel = feat[feat.type == 'CDS']
cds_iv = {g: sorted(zip(d.start, d.end)) for g, d in cds_sel.groupby('gene')}

STOPS = {'TAA', 'TAG', 'TGA'}

def tx_seq(chrom, strand, ivs):
    ivs = sorted(ivs); s = ''; pos = []
    for a, b in ivs:
        s += genome[chrom][a-1:b]; pos.extend(range(a-1, b))
    if strand == '-': s = rc(s); pos = pos[::-1]
    return s, np.array(pos)

down_cache = {}
def downstream(g):
    if g in down_cache: return down_cache[g]
    r = gbi.loc[g]; ch, st = r.chrom, r.strand
    cs, _ = tx_seq(ch, st, cds_iv[g])
    if st == '+': d = genome[ch][int(r.stop)+1:int(r.stop)+2001]
    else: d = rc(genome[ch][max(0, int(r.stop)-2000):int(r.stop)])
    down_cache[g] = (cs, d); return cs, d

# Build ATG context PWM
ctx = []
for g, r in gbi.iterrows():
    ch = r.chrom; a = int(r.atg)
    up = genome[ch][a-6:a] if r.strand == '+' else rc(genome[ch][a+1:a+7])
    cs, _ = downstream(g); ctx.append(up + cs[:6])
ctx = np.array([list(x) for x in ctx])
pos_labels = [-6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6]
freq = pd.DataFrame({p: pd.Series(ctx[:, i]).value_counts(normalize=True) for i, p in enumerate(pos_labels)}).fillna(0).loc[list('ACGT')]

# Parse refined UTR GFF3
RG = result('refined_utr_annotation.gff3')
rr = []
for l in open(RG):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if f[2] in ('five_prime_UTR', 'CDS', 'three_prime_UTR'):
        a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
        rr.append((f[0], f[2], int(f[3]), int(f[4]), f[6], a['Parent'].rsplit('.', 1)[0]))
rg = pd.DataFrame(rr, columns=['chrom', 'type', 'start', 'end', 'strand', 'gene'])
u5r = rg[rg.type == 'five_prime_UTR']

# UTR background for PWM
refined_genes = set(gb.loc[gb.final_src.isin(['long-read', 'model']), 'gene_id'])
u5_ref_iv = {g: list(zip(d.start, d.end)) for g, d in u5r.groupby('gene')}
u5v = feat[feat.type == 'five_prime_UTR']
u5_v71_iv = {g: list(zip(d.start, d.end)) for g, d in u5v.groupby('gene')}

utrseqs = {g: tx_seq(gbi.loc[g, 'chrom'], gbi.loc[g, 'strand'], u5_ref_iv[g])[0] for g in refined_genes if g in u5_ref_iv}
allu = ''.join(utrseqs.values()); bg = {b: allu.count(b)/len(allu) for b in 'ACGT'}

PWM_POS = [-6, -5, -4, -3, -2, -1, 4, 5]
pwm = {p: {b: np.log2((freq.loc[b, p]*len(gbi)+0.5)/(len(gbi)+2)/bg[b]) for b in 'ACGT'} for p in PWM_POS}

# Kozak threshold
main_sc = []
for row in ctx:
    c = ''.join(row); main_sc.append(sum(pwm[p][c[p+6 if p < 0 else p+5]] for p in PWM_POS))
KZ_THR = float(np.percentile(main_sc, 25))

# Census function
def census(g, ivs):
    r = gbi.loc[g]; ch, st = r.chrom, r.strand
    u, pos = tx_seq(ch, st, ivs); L = len(u)
    cs, d = downstream(g); full = u + cs + d
    out = []
    for i in range(L - 2):
        if u[i:i+3] != 'ATG': continue
        j = i; stop_end = None
        while j + 3 <= len(full):
            if full[j:j+3] in STOPS: stop_end = j + 3; break
            j += 3
        inframe = ((L - i) % 3 == 0)
        if stop_end is not None and stop_end <= L: cls = 'contained'
        elif inframe and (stop_end is None or stop_end > L): cls = 'N_ext'
        elif stop_end is None: cls = 'overlap_nostop'
        else: cls = 'overlapping'
        ncod = ((stop_end - 3 - i) // 3) if stop_end is not None else np.nan
        c = ''.join(full[k] if 0 <= k < len(full) else 'N' for k in range(i-6, i+6)) if i >= 6 else 'N'*(6-i)+full[0:i+6]
        sc = sum(pwm[p].get(c[p+6 if p < 0 else p+5], 0) for p in PWM_POS)
        out.append(dict(gene_id=g, upos=i-L, gpos=int(pos[i]), cls=cls, ncodons=ncod,
                        m3=c[3], p4=c[9], kozak_pwm=sc, dist_to_cap=i))
    return L, u, out

recs = {}; uaug_rows = []
for ver, ivd, gs in [('refined', u5_ref_iv, refined_genes), ('v71', u5_v71_iv, set(u5_v71_iv))]:
    for g in gs:
        if g not in ivd: continue
        L, u, out = census(g, ivd[g])
        for o in out: o['version'] = ver
        uaug_rows += out
        recs[(ver, g)] = dict(L=L, AT=(u.count('A')+u.count('T'))/max(L, 1))
ua = pd.DataFrame(uaug_rows)

ua['is_uorf'] = ua.cls.isin(['contained', 'overlapping']) & (ua.ncodons >= 2)
ua['startstop'] = (ua.cls == 'contained') & (ua.ncodons == 1)

cat = ua[['version', 'gene_id', 'gpos', 'upos', 'dist_to_cap', 'cls', 'ncodons', 'is_uorf', 'm3', 'p4', 'kozak_pwm']].rename(columns={'gpos': 'uAUG_genomic_pos0', 'upos': 'pos_rel_mainATG', 'ncodons': 'uORF_codons_incl_AUG'})
cat.to_csv('uaug_catalog.tsv.gz', sep='\t', index=False, float_format='%.3f')
