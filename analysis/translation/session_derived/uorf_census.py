"""Session-derived computation of uorf_census.tsv.

Recovered from the record of the interactive analysis session in which uorf_census.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: gene_boundaries_v3_all_genes.tsv, Lab/PlasmoDB/PlasmoDB-68/PlasmoDB-68_Pfalciparum3D7.gff, refined_utr_annotation.gff3, refined_utr_annotation.tsv, Lab/IGV_longread_transcript_models/PlasmoDB-68_Pfalciparum3D7_Genome.fasta
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import LAB_DIR, result  # paths: config/paths.yaml
import pandas as pd, numpy as np, re, warnings
warnings.filterwarnings('ignore')

gb = pd.read_csv(result('gene_boundaries_v3_all_genes.tsv'), sep="\t")

def read_fa(p):
    seqs={}; name=None; buf=[]
    for l in open(p):
        if l.startswith('>'):
            if name: seqs[name]=''.join(buf)
            name=l[1:].split()[0].split('|')[0].strip(); buf=[]
        else: buf.append(l.strip().upper())
    seqs[name]=''.join(buf); return seqs

genome = read_fa(f'{LAB_DIR}/IGV_longread_transcript_models/PlasmoDB-68_Pfalciparum3D7_Genome.fasta')
comp = str.maketrans('ACGTN','TGCAN')
def rc(s): return s.translate(comp)[::-1]

G = f'{LAB_DIR}/PlasmoDB/PlasmoDB-68/PlasmoDB-68_Pfalciparum3D7.gff'
rows = []
for l in open(G):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if len(f) < 9: continue
    if f[2] in ('CDS','five_prime_UTR','three_prime_UTR','exon','mRNA'):
        a = dict(x.split('=',1) for x in f[8].split(';') if '=' in x)
        rows.append((f[0],f[2],int(f[3]),int(f[4]),f[6],a.get('Parent',''),a.get('ID','')))
g68 = pd.DataFrame(rows, columns=['chrom','type','start','end','strand','parent','id'])

tx2gene = {}
for l in open(G):
    if '\tmRNA\t' in l:
        a = dict(x.split('=',1) for x in l.rstrip().split('\t')[8].split(';') if '=' in x)
        tx2gene[a['ID']] = a['Parent']
g68['gene'] = g68.parent.map(tx2gene)
g68.loc[g68.type=='mRNA','gene'] = g68.loc[g68.type=='mRNA','parent']

genes = set(gb.gene_id)
cds = g68[(g68.type=='CDS') & g68.gene.isin(genes)]
c1 = cds[cds.parent.str.endswith('.1')]

gbi = gb.set_index('gene_id')

g68g = g68[g68.gene.isin(genes)]
cdsx = g68g[g68g.type=='CDS'].groupby('parent').agg(gene=('gene','first'),s=('start','min'),e=('end','max'),strand=('strand','first'))
cdsx = cdsx.join(gbi[['atg','stop']], on='gene')
cdsx['ok'] = np.where(cdsx.strand=='+', (cdsx.s-1==cdsx.atg)&(cdsx.e-1==cdsx.stop), (cdsx.e-1==cdsx.atg)&(cdsx.s-1==cdsx.stop))
cdsx['is1'] = cdsx.index.str.endswith('.1')
sel = cdsx.reset_index().sort_values(['gene','ok','is1'],ascending=[True,False,False]).drop_duplicates('gene').set_index('gene')
tx_of = sel.parent.to_dict()
feat = g68g[g68g.parent.isin(set(tx_of.values()))]
u5v = feat[feat.type=='five_prime_UTR']

cds_iv = {g:sorted(zip(d.start,d.end)) for g,d in g68g[g68g.type=='CDS'].groupby('gene')}

def tx_seq(chrom, strand, ivs):
    ivs = sorted(ivs); s = ''; pos = []
    for a,b in ivs:
        s += genome[chrom][a-1:b]; pos.extend(range(a-1,b))
    if strand == '-': s = rc(s); pos = pos[::-1]
    return s, np.array(pos)

down_cache = {}
STOPS = {'TAA','TAG','TGA'}

def downstream(g):
    if g in down_cache: return down_cache[g]
    r = gbi.loc[g]; ch,st = r.chrom, r.strand
    cs,_ = tx_seq(ch, st, cds_iv[g])
    if st == '+': d = genome[ch][int(r.stop)+1:int(r.stop)+2001]
    else: d = rc(genome[ch][max(0,int(r.stop)-2000):int(r.stop)])
    down_cache[g] = (cs,d); return cs,d

ctx = []
for g,r in gbi.iterrows():
    ch = r.chrom; a = int(r.atg)
    up = genome[ch][a-6:a] if r.strand=='+' else rc(genome[ch][a+1:a+7])
    cs,_ = downstream(g); ctx.append(up+cs[:6])
ctx = np.array([list(x) for x in ctx])
pos_labels = [-6,-5,-4,-3,-2,-1,1,2,3,4,5,6]
freq = pd.DataFrame({p:pd.Series(ctx[:,i]).value_counts(normalize=True) for i,p in enumerate(pos_labels)}).fillna(0).loc[list('ACGT')]

RG = result('refined_utr_annotation.gff3')
rr = []
for l in open(RG):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if f[2] in ('five_prime_UTR','CDS','three_prime_UTR'):
        a = dict(x.split('=',1) for x in f[8].split(';') if '=' in x)
        rr.append((f[0],f[2],int(f[3]),int(f[4]),f[6],a['Parent'].rsplit('.',1)[0]))
rg = pd.DataFrame(rr, columns=['chrom','type','start','end','strand','gene'])
u5r = rg[rg.type=='five_prime_UTR']

refined_genes = set(gb.loc[gb.final_src.isin(['long-read','model']),'gene_id'])

u5_ref_iv = {g:list(zip(d.start,d.end)) for g,d in u5r.groupby('gene')}
u5_v71_iv = {g:list(zip(d.start,d.end)) for g,d in u5v.groupby('gene')}

utrseqs = {g:tx_seq(gbi.loc[g,'chrom'],gbi.loc[g,'strand'],u5_ref_iv[g])[0] for g in refined_genes if g in u5_ref_iv}
allu = ''.join(utrseqs.values()); bg = {b:allu.count(b)/len(allu) for b in 'ACGT'}

PWM_POS = [-6,-5,-4,-3,-2,-1,4,5]
pwm = {p:{b:np.log2((freq.loc[b,p]*len(gbi)+0.5)/(len(gbi)+2)/bg[b]) for b in 'ACGT'} for p in PWM_POS}

main_sc = []
for row in ctx:
    c = ''.join(row); main_sc.append(sum(pwm[p][c[p+6 if p<0 else p+5]] for p in PWM_POS))
KZ_THR = float(np.percentile(main_sc, 25))

def census(g, ivs):
    r = gbi.loc[g]; ch,st = r.chrom, r.strand
    u,pos = tx_seq(ch, st, ivs); L = len(u)
    cs,d = downstream(g); full = u+cs+d
    out = []
    for i in range(L-2):
        if u[i:i+3] != 'ATG': continue
        j = i; stop_end = None
        while j+3 <= len(full):
            if full[j:j+3] in STOPS: stop_end = j+3; break
            j += 3
        inframe = ((L-i)%3==0)
        if stop_end is not None and stop_end <= L: cls = 'contained'
        elif inframe and (stop_end is None or stop_end > L): cls = 'N_ext'
        elif stop_end is None: cls = 'overlap_nostop'
        else: cls = 'overlapping'
        ncod = ((stop_end-3-i)//3) if stop_end is not None else np.nan
        c = ''.join(full[k] if 0<=k<len(full) else 'N' for k in range(i-6,i+6)) if i>=6 else 'N'*(6-i)+full[0:i+6]
        sc = sum(pwm[p].get(c[p+6 if p<0 else p+5],0) for p in PWM_POS)
        out.append(dict(gene_id=g, upos=i-L, gpos=int(pos[i]), cls=cls, ncodons=ncod,
                        m3=c[3], p4=c[9], kozak_pwm=sc, dist_to_cap=i))
    return L, u, out

recs = {}; uaug_rows = []
for ver,ivd,gs in [('refined',u5_ref_iv,refined_genes),('v71',u5_v71_iv,set(u5_v71_iv))]:
    for g in gs:
        if g not in ivd: continue
        L,u,out = census(g, ivd[g])
        for o in out: o['version'] = ver
        uaug_rows += out
        recs[(ver,g)] = dict(L=L, AT=(u.count('A')+u.count('T'))/max(L,1))
ua = pd.DataFrame(uaug_rows)

ua['is_uorf'] = ua.cls.isin(['contained','overlapping']) & (ua.ncodons>=2)
ua['startstop'] = (ua.cls=='contained') & (ua.ncodons==1)

def summarize(ver):
    d = ua[ua.version==ver]; gs = [g for (v,g) in recs if v==ver]
    base = pd.DataFrame({'gene_id':gs,'utr5_len':[recs[(ver,g)]['L'] for g in gs],'utr5_AT':[recs[(ver,g)]['AT'] for g in gs]}).set_index('gene_id')
    d = d.assign(c2=(d.cls=='contained')&(d.ncodons>=2), o2=(d.cls=='overlapping')&(d.ncodons>=2), nx=d.cls=='N_ext',
                 m3A=d.m3=='A', rg=d.m3.isin(['A','G'])&(d.p4=='G'), strong=d.is_uorf&(d.kozak_pwm>=KZ_THR),
                 uolen=d.ncodons.where(d.is_uorf))
    grp = d.groupby('gene_id')
    s = pd.DataFrame({'n_uAUG':grp.size(),'n_uORF':grp.is_uorf.sum(),'n_uORF_contained':grp.c2.sum(),'n_oORF_overlapping':grp.o2.sum(),
      'n_startstop':grp.startstop.sum(),'n_Nterm_ext_inframe':grp.nx.sum(),'max_uORF_codons':grp.uolen.max(),'median_uORF_codons':grp.uolen.median(),
      'n_uAUG_m3A':grp.m3A.sum(),'n_uAUG_m3R_p4G':grp.rg.sum(),'max_uAUG_kozak_pwm':grp.kozak_pwm.max(),'n_uORF_strongctx':grp.strong.sum(),
      'dist_uAUG_nearest_ATG':grp.upos.max().abs()})
    out = base.join(s)
    cnt = ['n_uAUG','n_uORF','n_uORF_contained','n_oORF_overlapping','n_startstop','n_Nterm_ext_inframe','n_uAUG_m3A','n_uAUG_m3R_p4G','n_uORF_strongctx']
    out[cnt] = out[cnt].fillna(0).astype(int); out['oORF_any'] = out.n_oORF_overlapping > 0
    return out

cen_ref = summarize('refined'); cen_v71 = summarize('v71')

rt = pd.read_csv(result('refined_utr_annotation.tsv'), sep="\t").set_index('gene_id')

R = ua[ua.version=='refined']; V = ua[ua.version=='v71']
both_genes = sorted(refined_genes)
Rd = {g:d for g,d in R.groupby('gene_id')}
Vd = {g:d for g,d in V.groupby('gene_id')}
comp_rows = []
for g in both_genes:
    r = gbi.loc[g]; st = r.strand
    rset = Rd.get(g, R.iloc[:0]); vset = Vd.get(g, V.iloc[:0])
    rp = set(rset.gpos); vp = set(vset.gpos)
    ru = set(rset.loc[rset.is_uorf,'gpos']); vu = set(vset.loc[vset.is_uorf,'gpos'])
    tss = int(r.final_TSS)
    upstream = lambda p: (p<tss) if st=='+' else (p>tss)
    lost = vp-rp; lostu = vu-ru
    comp_rows.append(dict(gene_id=g, n_uAUG_shared=len(rp&vp), n_uAUG_gained=len(rp-vp), n_uAUG_lost=len(lost),
        n_uAUG_lost_trimmed=sum(upstream(p) for p in lost),
        n_uORF_shared=len(ru&vu), n_uORF_gained=len(ru-vu), n_uORF_lost=len(lostu), n_uORF_lost_trimmed=sum(upstream(p) for p in lostu),
        n_oORF_gained=len(set(rset.loc[rset.is_uorf&(rset.cls=='overlapping'),'gpos'])-set(vset.loc[vset.is_uorf&(vset.cls=='overlapping'),'gpos'])),
        n_oORF_lost=len(set(vset.loc[vset.is_uorf&(vset.cls=='overlapping'),'gpos'])-set(rset.loc[rset.is_uorf&(rset.cls=='overlapping'),'gpos']))))
cmpdf = pd.DataFrame(comp_rows).set_index('gene_id')
cmpdf = cmpdf.join(rt[['utr5_change_vs_v71']])
cmpdf['utr5_len_ref'] = cen_ref.utr5_len; cmpdf['utr5_len_v71'] = cen_v71.utr5_len.reindex(cmpdf.index).fillna(0)

cen_cols = ['utr5_len','utr5_AT','n_uAUG','n_uORF','n_uORF_contained','n_oORF_overlapping','n_startstop','n_Nterm_ext_inframe','max_uORF_codons','median_uORF_codons','n_uAUG_m3A','n_uAUG_m3R_p4G','max_uAUG_kozak_pwm','n_uORF_strongctx','dist_uAUG_nearest_ATG','oORF_any']
out = gb.set_index('gene_id')[['chrom','strand','atg','lr_group','final_src','final_TSS','annot_tss']].copy()
out['utr5_change_vs_v71'] = rt.utr5_change_vs_v71
out['refined_utr5_source'] = rt.TSS_evidence
for v,cen in [('ref',cen_ref),('v71',cen_v71)]:
    out = out.join(cen[cen_cols].add_prefix(v+'_'))
out = out.join(cmpdf.drop(columns=['utr5_change_vs_v71','utr5_len_ref','utr5_len_v71']))
out['ref_uORF_status_vs_v71'] = np.select(
    [out.ref_n_uORF.isna(),
     (out.ref_n_uORF>0)&(out.v71_n_uORF.fillna(0)>0),
     (out.ref_n_uORF>0),
     (out.v71_n_uORF.fillna(0)>0)],
    ['no_refined_UTR','uORF_both','uORF_gained','uORF_lost'],
    'uORF_neither')
out.index.name = 'gene_id'
out.to_csv('uorf_census.tsv', sep='\t', float_format='%.4g')
