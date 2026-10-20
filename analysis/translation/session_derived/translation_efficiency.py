"""Session-derived computation of translation_efficiency.tsv.

Recovered from the record of the interactive analysis session in which translation_efficiency.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: (none recorded)
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import LAB_DIR, result  # paths: config/paths.yaml
import pandas as pd
import numpy as np
import pyBigWig
from scipy import stats as ss
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

# Load gene boundaries
gb = pd.read_csv(result('gene_boundaries_v3_all_genes.tsv'), sep="\t")

# Load genome
def read_fa(p):
    seqs = {}; name = None; buf = []
    for l in open(p):
        if l.startswith('>'):
            if name: seqs[name] = ''.join(buf)
            name = l[1:].split()[0].split('|')[0].strip(); buf = []
        else: buf.append(l.strip().upper())
    seqs[name] = ''.join(buf); return seqs

genome = read_fa(f"{LAB_DIR}/IGV_longread_transcript_models/PlasmoDB-68_Pfalciparum3D7_Genome.fasta")
comp = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(comp)[::-1]

G = f"{LAB_DIR}/PlasmoDB/PlasmoDB-68/PlasmoDB-68_Pfalciparum3D7.gff"

# Parse GFF
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
cds = g68[(g68.type == 'CDS') & g68.gene.isin(genes)]
c1 = cds[cds.parent.str.endswith('.1')]

gbi = gb.set_index('gene_id')

# Build CDS intervals
cds_iv = {g: sorted(zip(d.start, d.end)) for g, d in cds.groupby('gene')}

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

# Load refined UTR annotation
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

# Load refined UTR annotation TSV
rt = pd.read_csv(result('refined_utr_annotation.tsv'), sep="\t").set_index('gene_id')

# Choose transcript per gene matching gb atg
g68g = g68[g68.gene.isin(genes)]
cdsx = g68g[g68g.type == 'CDS'].groupby('parent').agg(gene=('gene', 'first'), s=('start', 'min'), e=('end', 'max'), strand=('strand', 'first'))
cdsx = cdsx.join(gbi[['atg', 'stop']], on='gene')
cdsx['ok'] = np.where(cdsx.strand == '+', (cdsx.s-1 == cdsx.atg) & (cdsx.e-1 == cdsx.stop), (cdsx.e-1 == cdsx.atg) & (cdsx.s-1 == cdsx.stop))
cdsx['is1'] = cdsx.index.str.endswith('.1')
sel = cdsx.reset_index().sort_values(['gene', 'ok', 'is1'], ascending=[True, False, False]).drop_duplicates('gene').set_index('gene')
tx_of = sel.parent.to_dict()
feat = g68g[g68g.parent.isin(set(tx_of.values()))]
u5v = feat[feat.type == 'five_prime_UTR']

# Build UTR sequence dicts
u5_ref_iv = {g: list(zip(d.start, d.end)) for g, d in u5r.groupby('gene')}
u5_v71_iv = {g: list(zip(d.start, d.end)) for g, d in u5v.groupby('gene')}
refined_genes = set(gb.loc[gb.final_src.isin(['long-read', 'model']), 'gene_id'])

utrseqs = {g: tx_seq(gbi.loc[g, 'chrom'], gbi.loc[g, 'strand'], u5_ref_iv[g])[0] for g in refined_genes if g in u5_ref_iv}

# Kozak PWM
ctx = []
for g, r in gbi.iterrows():
    ch = r.chrom; a = int(r.atg)
    up = genome[ch][a-6:a] if r.strand == '+' else rc(genome[ch][a+1:a+7])
    cs, _ = downstream(g); ctx.append(up+cs[:6])
ctx = np.array([list(x) for x in ctx])
pos_labels = [-6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6]
freq = pd.DataFrame({p: pd.Series(ctx[:, i]).value_counts(normalize=True) for i, p in enumerate(pos_labels)}).fillna(0).loc[list('ACGT')]
allu = ''.join(utrseqs.values()); bg = {b: allu.count(b)/len(allu) for b in 'ACGT'}
PWM_POS = [-6, -5, -4, -3, -2, -1, 4, 5]
pwm = {p: {b: np.log2((freq.loc[b, p]*len(gbi)+0.5)/(len(gbi)+2)/bg[b]) for b in 'ACGT'} for p in PWM_POS}

main_sc = []
for row in ctx:
    c = ''.join(row); main_sc.append(sum(pwm[p][c[p+6 if p < 0 else p+5]] for p in PWM_POS))
KZ_THR = float(np.percentile(main_sc, 25))

# uAUG census
def census(g, ivs):
    r = gbi.loc[g]; ch, st = r.chrom, r.strand
    u, pos = tx_seq(ch, st, ivs); L = len(u)
    cs, d = downstream(g); full = u+cs+d
    out = []
    for i in range(L-2):
        if u[i:i+3] != 'ATG': continue
        j = i; stop_end = None
        while j+3 <= len(full):
            if full[j:j+3] in STOPS: stop_end = j+3; break
            j += 3
        inframe = ((L-i) % 3 == 0)
        if stop_end is not None and stop_end <= L: cls = 'contained'
        elif inframe and (stop_end is None or stop_end > L): cls = 'N_ext'
        elif stop_end is None: cls = 'overlap_nostop'
        else: cls = 'overlapping'
        ncod = ((stop_end-3-i)//3) if stop_end is not None else np.nan
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

def summarize(ver):
    d = ua[ua.version == ver]; gs = [g for (v, g) in recs if v == ver]
    base = pd.DataFrame({'gene_id': gs, 'utr5_len': [recs[(ver, g)]['L'] for g in gs], 'utr5_AT': [recs[(ver, g)]['AT'] for g in gs]}).set_index('gene_id')
    d = d.assign(c2=(d.cls == 'contained') & (d.ncodons >= 2), o2=(d.cls == 'overlapping') & (d.ncodons >= 2), nx=d.cls == 'N_ext',
                 m3A=d.m3 == 'A', rg=d.m3.isin(['A', 'G']) & (d.p4 == 'G'), strong=d.is_uorf & (d.kozak_pwm >= KZ_THR),
                 uolen=d.ncodons.where(d.is_uorf))
    grp = d.groupby('gene_id')
    s = pd.DataFrame({'n_uAUG': grp.size(), 'n_uORF': grp.is_uorf.sum(), 'n_uORF_contained': grp.c2.sum(), 'n_oORF_overlapping': grp.o2.sum(),
        'n_startstop': grp.startstop.sum(), 'n_Nterm_ext_inframe': grp.nx.sum(), 'max_uORF_codons': grp.uolen.max(), 'median_uORF_codons': grp.uolen.median(),
        'n_uAUG_m3A': grp.m3A.sum(), 'n_uAUG_m3R_p4G': grp.rg.sum(), 'max_uAUG_kozak_pwm': grp.kozak_pwm.max(), 'n_uORF_strongctx': grp.strong.sum(),
        'dist_uAUG_nearest_ATG': grp.upos.max().abs()})
    out = base.join(s)
    cnt = ['n_uAUG', 'n_uORF', 'n_uORF_contained', 'n_oORF_overlapping', 'n_startstop', 'n_Nterm_ext_inframe', 'n_uAUG_m3A', 'n_uAUG_m3R_p4G', 'n_uORF_strongctx']
    out[cnt] = out[cnt].fillna(0).astype(int); out['oORF_any'] = out.n_oORF_overlapping > 0
    return out

cen_ref = summarize('refined'); cen_v71 = summarize('v71')

# Build cov table
cov = gb.set_index('gene_id')[['chrom', 'strand', 'lr_group', 'final_src', 'utr5_final', 'utr3_final', 'utr5_annot', 'utr3_annot']].copy()
cov['cds_len'] = [len(downstream(g)[0]) for g in cov.index]
cov['cds_GC'] = [(lambda s: (s.count('G')+s.count('C'))/len(s))(downstream(g)[0]) for g in cov.index]
cov = cov.join(rt[['utr5_exonic_len', 'utr3_exonic_len', 'utr5_change_vs_v71', 'utr3_change_vs_v71', 'TSS_evidence', 'TTS_evidence']])

dt = pd.read_csv(result('gene_daft_expression_truth.tsv.gz'), sep="\t")
d2 = dt.set_index('gene_id'); cov = cov.join(d2[['DAFT_max', 'Qd', 'h0', 'h16', 'h32', 'h40']])
cov['daft_36'] = (cov.h32+cov.h40)/2

cen_cols = ['utr5_len', 'utr5_AT', 'n_uAUG', 'n_uORF', 'n_uORF_contained', 'n_oORF_overlapping', 'n_startstop', 'n_Nterm_ext_inframe', 'max_uORF_codons', 'median_uORF_codons', 'n_uAUG_m3A', 'n_uAUG_m3R_p4G', 'max_uAUG_kozak_pwm', 'n_uORF_strongctx', 'dist_uAUG_nearest_ATG', 'oORF_any']
for v, cen in [('ref', cen_ref), ('v71', cen_v71)]:
    c2 = cen.add_prefix(v+'_'); cov = cov.join(c2)

u3v = feat[feat.type == 'three_prime_UTR'].assign(L=lambda d: d.end-d.start+1).groupby('gene').L.sum()
cov['v71_utr3_len'] = u3v.reindex(cov.index)
cov['ref_utr3_len'] = cov.utr3_exonic_len.where(cov.TTS_evidence.isin(['long-read', 'model']))

# Build upstream intergenic and family
gl = []
for l in open(G):
    if l.startswith('#'): continue
    f = l.split('\t')
    if len(f) > 8 and f[2] in ('protein_coding_gene', 'ncRNA_gene', 'pseudogene'):
        gl.append((f[0], int(f[3]), int(f[4]), f[8].split('ID=')[1].split(';')[0].strip()))
gl = pd.DataFrame(gl, columns=['chrom', 's', 'e', 'id'])
upd = {}
for ch, dg in gl.groupby('chrom'):
    dg = dg.sort_values('s'); S_ = dg.s.values; E_ = dg.e.values; ids = dg.id.values
    for i, g in enumerate(ids):
        if g not in gbi.index: continue
        r_ = gbi.loc[g]; a = int(r_.atg)
        if r_.strand == '+':
            prev = E_[:i][E_[:i] < a+1]; upd[g] = a-(prev.max() if len(prev) else 0)
        else:
            nxt = S_[i+1:][S_[i+1:] > a+1]; upd[g] = (nxt.min()-1-a) if len(nxt) else np.nan
cov['upstream_intergenic'] = pd.Series(upd)
cov['family'] = np.select([cov.index.isin(gb.gene_id[gb.desc.str.contains('erythrocyte membrane protein 1|rifin|stevor|Pfmc-2TM|surfin', case=False, na=False)])], ['variant_family'], 'other')

# Polysome bigwig signal
TR = f"{LAB_DIR}/lncRNA project tracks/translatome/"
bwf = {'0h': 'Polysome_Ring_0h_CPM_final.bw', '18h': 'Polysome_Trophozoite_18h_CPM_final.bw', '36h': 'Polysome_Schizont_36h_CPM_final.bw'}
psig = {}
for st, f in bwf.items():
    b = pyBigWig.open(TR+f); vals = {}
    for g, ivs in cds_iv.items():
        ch = gbi.loc[g, 'chrom']; tot = 0; L = 0
        for a, e in ivs:
            v = b.stats(ch, a-1, e, type='mean', exact=True)[0]; tot += (v or 0)*(e-a+1); L += e-a+1
        vals[g] = tot/L
    psig[st] = pd.Series(vals); b.close()
psig = pd.DataFrame(psig)

# Proxy TE
proxy = pd.DataFrame(index=cov.index)
STAGE_DAFT = {'0h': 'h0', '18h': 'h16', '36h': 'daft_36'}
for st in ['0h', '18h', '36h']:
    dd = cov[STAGE_DAFT[st]]
    ok = dd >= dd.quantile(0.2)
    lt = np.log2(psig[st]+0.05)-np.log2(dd+0.01)
    proxy['log2TE_'+st] = (lt-lt[ok].median()).where(ok); proxy['ss_CPM_'+st] = dd

# Load count data
LIB = {'SRR836068': 'poly_0h', 'SRR836069': 'poly_18h', 'SRR836070': 'poly_36h_r1', 'SRR836071': 'poly_36h_r2',
       'SRR835998': 'ss_0h', 'SRR836066': 'ss_18h', 'SRR836067': 'ss_36h'}

def load_counts(path):
    fc = pd.read_csv(path, sep='\t', comment='#').set_index('Geneid')
    cols = [c for c in fc.columns if c.endswith('.bam')]
    fc = fc[cols]
    fc.columns = [LIB[c.split('/')[-1].split('.')[0]] for c in cols]
    return fc

def te_table(cnt, min_ss=20):
    cpm = cnt / cnt.sum() * 1e6
    te = pd.DataFrame(index=cnt.index)
    for st, p, s in [('0h', ['poly_0h'], 'ss_0h'), ('18h', ['poly_18h'], 'ss_18h'), ('36h', ['poly_36h_r1', 'poly_36h_r2'], 'ss_36h')]:
        pc = cpm[p].mean(axis=1)
        pr = cnt[p].sum(axis=1)
        ps_p = np.mean([0.5 / cnt[x].sum() * 1e6 for x in p]); ps_s = 0.5 / cnt[s].sum() * 1e6
        lt = np.log2((pc + ps_p) / (cpm[s] + ps_s))
        ok = cnt[s] >= min_ss
        lt = lt.where(ok)
        te['poly_reads_' + st] = pr; te['ss_reads_' + st] = cnt[s]
        te['poly_CPM_' + st] = pc; te['ss_CPM_' + st] = cpm[s]
        te['log2TE_raw_' + st] = lt
        te['log2TE_' + st] = lt - lt.median()
        te['TE_pass_' + st] = ok
    return te

cnt0 = load_counts('counts/cds_counts_s0.txt')
cnt0 = cnt0.loc[[g for g in cnt0.index if g in gbi.index]]
te = te_table(cnt0, min_ss=20)

# Build final translation_efficiency.tsv
tef = te.copy()
tef = tef.join(cnt0.add_prefix('cds_pairs_'))
for st in ['0h', '18h', '36h']:
    tef['proxy_log2TE_polyBW_over_DAFT_'+st] = proxy['log2TE_'+st]
    tef['DAFT_nearest_'+st] = cov[STAGE_DAFT[st]]
tef = tef.join(cov[['lr_group', 'final_src', 'cds_len', 'cds_GC', 'upstream_intergenic', 'DAFT_max', 'Qd']])
tef['stage_quality_note'] = '0h(ring) polysome library low quality: 22.4% alignment, 1.05M CDS pairs, r(log CPM) with matched ss_0h 0.336; interpret 0h TE as unreliable'
tef.index.name = 'gene_id'
tef.to_csv('translation_efficiency.tsv', sep='\t', float_format='%.5g')
