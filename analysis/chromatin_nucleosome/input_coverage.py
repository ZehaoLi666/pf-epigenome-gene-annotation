from pfconfig import LAB_DIR  # paths: config/paths.yaml
import numpy as np, pandas as pd, pyBigWig, glob
TR=f"{LAB_DIR}/lncRNA project tracks"; IGV=f"{LAB_DIR}/IGV_longread_transcript_models"
ST=['Ring','Trophozoite','Schizont']
TRACKS={'Stranded RNA (DAFT-seq), gene strand':('daft_same',None),'Stranded RNA (DAFT-seq), opposite strand':('daft_opp',None),
 'Polysome profiling':('st',f"{TR}/translatome/Polysome_{{s}}_*_CPM_final.bw"),'Nascent RNA (GRO-seq)':('st',f"{TR}/nascent_transcription/GRO_{{s}}_mean_normalized.bw"),
 'H2A.Z':('st',f"{TR}/histone_marks/H2AZ_{{s}}_minusInput.bw"),'H2B.Z':('st',f"{TR}/histone_marks/H2BZ_{{s}}_minusInput.bw"),
 'H3K9ac':('st',f"{TR}/histone_marks/H3K9ac_{{s}}_minusInput.bw"),'H3K27ac':('st',f"{TR}/histone_marks/H3K27ac_{{s}}_minusInput.bw"),
 'H3K4me3':('st',f"{TR}/histone_marks/H3K4me3_{{s}}_minusInput.bw"),'H3K4me1':('st',f"{TR}/histone_marks/H3K4me1_{{s}}_minusInput.bw"),
 'H3':('st',f"{TR}/histone_marks/H3_{{s}}_minusInput.bw"),'H3K9me3':('st',f"{TR}/histone_marks/H3K9me3_{{s}}_minusInput.bw"),
 'HP1':('st',f"{TR}/machinery/HP1_{{s}}_minusInput.bw"),'Pol II':('pol',None),
 'ATAC-seq':('st',f"{TR}/accessibility_nucleosome/ATAC_{{s}}_minusInput.bw"),'MNase-seq':('st',f"{TR}/accessibility_nucleosome/MNase_{{s}}_CPM.bw"),
 'ApiAP2 (mean of 13)':('ap2',None),'CAGE (validation only)':('cage',None)}
def files(kind,pat):
    if kind=='st': return [glob.glob(pat.format(s=s))[0] for s in ST]
    if kind=='pol': return sorted(glob.glob(f"{TR}/machinery/PolII_*_minusInput.bw"))
    if kind=='ap2': return sorted(glob.glob(f"{TR}/apiap2_tf/AP2-*_minusInput.bw"))
def mean_over(paths,iv):
    out=np.zeros(len(iv))
    for p in paths:
        bw=pyBigWig.open(p); v=np.array([bw.stats(c,int(a),int(b),type='mean')[0] or 0.0 for c,a,b in iv]); bw.close(); out+=v
    return out/len(paths)
def stranded(prefix_map,iv,strands,same=True):
    out=np.zeros(len(iv)); bws={k:pyBigWig.open(v) for k,v in prefix_map.items()}
    for i,((c,a,b),s) in enumerate(zip(iv,strands)):
        k=s if same else ('-' if s=='+' else '+'); out[i]=bws[k].stats(c,int(a),int(b),type='mean')[0] or 0.0
    return out
def run(W):
    iv=list(zip(W.chrom,W.a,W.b)); res={}
    daft={'+':f"{IGV}/DAFTseq_pool_plus_CPM_10bp.bw",'-':f"{IGV}/DAFTseq_pool_minus_CPM_10bp.bw"}
    for name,(kind,pat) in TRACKS.items():
        if kind=='daft_same': res[name]=stranded(daft,iv,W.strand,True)
        elif kind=='daft_opp': res[name]=stranded(daft,iv,W.strand,False)
        elif kind=='cage':
            v=np.zeros(len(iv))
            for st in ['ring','trophozoite','schizont']:
                v+=stranded({'+':f"{TR}/cage_tss/CTSS_{st}_pos_CPM.bw",'-':f"{TR}/cage_tss/CTSS_{st}_neg_CPM.bw"},iv,W.strand,True)
            res[name]=v/3
        else: res[name]=mean_over(files(kind,pat),iv)
    return pd.DataFrame(res,index=W.index)
