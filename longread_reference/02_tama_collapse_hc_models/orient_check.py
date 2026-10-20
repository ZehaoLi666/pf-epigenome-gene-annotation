from pfconfig import HPCC_BIGDATA, HPCC_LAB  # paths: config/paths.yaml

import pysam, bisect, collections, sys, glob, os
GTF=f"{HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7.gtf"
genes=collections.defaultdict(dict)
for ln in open(GTF):
    if ln.startswith("#"): continue
    f=ln.rstrip("\n").split("\t")
    if len(f)<9 or f[2] not in ("exon","CDS"): continue
    gid=f[8].split('gene_id "')[1].split('"')[0]
    g=genes[f[0]].get(gid)
    s,e=int(f[3])-1,int(f[4])
    if g is None: genes[f[0]][gid]=[s,e,f[6]]
    else: g[0]=min(g[0],s); g[1]=max(g[1],e)
idx={}
for ch,d in genes.items():
    L=sorted(d.values()); idx[ch]=(L,[x[0] for x in L])
def strand_at(ch,a,b):
    if ch not in idx: return None
    L,st=idx[ch]; i=bisect.bisect_right(st,b)
    hits=set()
    for j in range(max(0,i-60),i):
        g=L[j]
        if g[1]>a and g[0]<b and min(g[1],b)-max(g[0],a)>=0.5*(b-a): hits.add(g[2])
    return hits.pop() if len(hits)==1 else None
P=f"{HPCC_BIGDATA}/GSE186109_longread"
bams=sorted(glob.glob(P+"/bam/*/*.sorted.bam"))+[
 f"{HPCC_LAB}/SRA_toolkit/Long_read_RNA_datasets/LeeVV_2021_nanopore/SRR11094274.no_long_introns.sorted.bam",
 f"{HPCC_LAB}/SRA_toolkit/Long_read_RNA_datasets/Marcus_Nanopore/ERR9746080.no_long_introns.sorted.bam",
 f"{HPCC_LAB}/SRA_toolkit/Long_read_RNA_datasets/Mengquan_long_read/merged.no_long_introns.sorted.bam"]
extra=sorted(glob.glob(P+"/bam_raw/*/*.sorted.bam"))
bams+=extra
out=open(P+"/diag_orientation.tsv","w")
out.write("bam\tn_assigned\tfrac_sense\tfrac_antisense\tn_spliced\tfrac_spliced_sense\tts_tag_frac_plus\n")
N=200000
for b in bams:
    if not os.path.exists(b): continue
    h=pysam.AlignmentFile(b,"rb")
    n=s=sp=sps=tsn=tsp=0; seen=0
    for r in h.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        seen+=1
        if seen%5: continue   # 1-in-5 sample across the file
        g=strand_at(r.reference_name,r.reference_start,r.reference_end)
        if g is None: continue
        rs="-" if r.is_reverse else "+"
        n+=1; s+= (rs==g)
        if "N" in r.cigarstring:
            sp+=1; sps+=(rs==g)
        if r.has_tag("ts"):
            tsn+=1; tsp+= (r.get_tag("ts")=="+")
        if n>=N: break
    name=b.replace(P+"/","").replace(f"{HPCC_LAB}/SRA_toolkit/Long_read_RNA_datasets/","")
    out.write(f"{name}\t{n}\t{s/max(n,1):.4f}\t{1-s/max(n,1):.4f}\t{sp}\t{sps/max(sp,1):.4f}\t{tsp/max(tsn,1):.4f}\n"); out.flush()
out.close(); open(P+"/logs/ORIENT_DONE","w").write("ok\n")
