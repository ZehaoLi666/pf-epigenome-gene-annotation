from pfconfig import LAB_DIR  # paths: config/paths.yaml
import pandas as pd, numpy as np, pyBigWig, re
LAB=LAB_DIR
IGV=f"{LAB}/IGV_longread_transcript_models"
GFF=f"{LAB}/PlasmoDB/PlasmoDB-68/PlasmoDB-68_Pfalciparum3D7.gff"
CAGE=f"{LAB}/lncRNA project tracks/cage_tss"
def load_gff():
    g=pd.read_csv(GFF,sep="\t",comment="#",header=None,names="chrom src type start end score strand phase attr".split(),dtype={'start':int,'end':int})
    g=g[g.type.isin(['protein_coding_gene','ncRNA_gene','mRNA','CDS','five_prime_UTR','three_prime_UTR','exon','gene'])].copy()
    g['ID']=g.attr.str.extract(r'ID=([^;]+)')
    g['Parent']=g.attr.str.extract(r'Parent=([^;]+)')
    g['gene']=np.where(g.type.str.endswith('gene'),g.ID,None)
    g.loc[g.type=='mRNA','gene']=g.loc[g.type=='mRNA','Parent']
    tx2g=dict(zip(g[g.type=='mRNA'].ID,g[g.type=='mRNA'].Parent))
    m=~g.type.str.endswith('gene')&(g.type!='mRNA')
    g.loc[m,'gene']=g.loc[m,'Parent'].map(lambda p: tx2g.get(str(p).split(',')[0],p))
    return g
def read_bed12(path):
    rows=[]
    for ln in open(path):
        if ln.startswith(('track','#')): continue
        f=ln.rstrip('\n').split('\t')
        r=dict(chrom=f[0],start=int(f[1]),end=int(f[2]),name=f[3],strand=f[5] if len(f)>5 else '.')
        if len(f)>=12:
            r['rgb']=f[8]; sz=[int(x) for x in f[10].strip(',').split(',')]; st=[int(x) for x in f[11].strip(',').split(',')]
            r['blocks']=[(r['start']+s,r['start']+s+z) for s,z in zip(st,sz)]
        else: r['blocks']=[(r['start'],r['end'])]
        rows.append(r)
    return pd.DataFrame(rows)
def bw_vals(path,chrom,a,b,nb):
    bw=pyBigWig.open(path); v=np.array(bw.stats(chrom,a,b,type='mean',nBins=nb),dtype=float); bw.close()
    return np.nan_to_num(v)
def daft(chrom,a,b,strand,nb=400):
    return bw_vals(f"{IGV}/DAFTseq_pool_{'plus' if strand=='+' else 'minus'}_CPM_10bp.bw",chrom,a,b,nb)
def cage(chrom,a,b,strand,nb=400):
    s='pos' if strand=='+' else 'neg'
    return sum(bw_vals(f"{CAGE}/CTSS_{st}_{s}_CPM.bw",chrom,a,b,nb) for st in ['ring','trophozoite','schizont'])
