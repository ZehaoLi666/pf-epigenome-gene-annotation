from pfconfig import HPCC_BIGDATA, HPCC_LAB  # paths: config/paths.yaml

import pysam, re, collections, bisect, sys
P=f"{HPCC_BIGDATA}/GSE186109_longread"
GTF=f"{HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7.gtf"
genes=collections.defaultdict(dict)
for ln in open(GTF):
    if ln.startswith("#"): continue
    f=ln.rstrip("\n").split("\t")
    if len(f)<9 or f[2]!="exon": continue
    gid=f[8].split('gene_id "')[1].split('"')[0]
    s,e=int(f[3])-1,int(f[4]); g=genes[f[0]].get(gid)
    if g is None: genes[f[0]][gid]=[s,e,f[6]]
    else: g[0]=min(g[0],s); g[1]=max(g[1],e)
idx={ch:(sorted(d.values()),[x[0] for x in sorted(d.values())]) for ch,d in genes.items()}
def gstrand(ch,a,b):
    if ch not in idx: return None
    L,st=idx[ch]; i=bisect.bisect_right(st,b); hits=set()
    for j in range(max(0,i-60),i):
        g=L[j]
        if g[1]>a and g[0]<b and min(g[1],b)-max(g[0],a)>=0.5*(b-a): hits.add(g[2])
    return hits.pop() if len(hits)==1 else None
comp=str.maketrans("ACGTNacgtn","TGCANtgcan")
polyA=re.compile("A{10,}"); polyT=re.compile("T{10,}")
out=open(P+"/by_sample/ont_mixed/oriented.fastq","w")
stats=collections.Counter()
for run in ("run1","run2"):
    h=pysam.AlignmentFile(f"{P}/bam/ont_mixed_ub/{run}.sorted.bam","rb")
    for r in h.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        seq=r.query_sequence; q=r.qual
        if seq is None: continue
        aln="-" if r.is_reverse else "+"
        spliced="N" in r.cigarstring
        ts=r.get_tag("ts") if r.has_tag("ts") else None
        t_ts=None
        if spliced and ts in ("+","-"):
            t_ts = aln if ts=="+" else ("-" if aln=="+" else "+")
        a3=bool(polyA.search(seq[-60:])); t5=bool(polyT.search(seq[:60]))
        t_pa = "+" if (a3 and not t5) else ("-" if (t5 and not a3) else None)
        g=gstrand(r.reference_name,r.reference_start,r.reference_end)
        if t_ts and t_pa: stats["ts_vs_polyA_agree" if t_ts==t_pa else "ts_vs_polyA_disagree"]+=1
        if t_ts and g: stats["ts_vs_gene_agree" if t_ts==g else "ts_vs_gene_disagree"]+=1
        if t_pa and g: stats["polyA_vs_gene_agree" if t_pa==g else "polyA_vs_gene_disagree"]+=1
        t = t_ts if t_ts else t_pa
        if t is None: stats["undetermined"]+=1; continue
        stats["oriented_by_"+("ts" if t_ts else "polyA")]+=1
        if t=="+": s2,q2=seq,q
        else: s2,q2=seq.translate(comp)[::-1], (q[::-1] if q else None)
        out.write(f"@{r.query_name}\n{s2}\n+\n{q2 if q2 else 'I'*len(s2)}\n")
out.close()
with open(P+"/diag_ont_orientation.tsv","w") as f:
    for k,v in sorted(stats.items()): f.write(f"{k}\t{v}\n")
open(P+"/logs/ONT_ORIENT_DONE","w").write("ok\n")
