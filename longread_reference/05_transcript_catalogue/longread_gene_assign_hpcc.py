
import sys,gzip,csv,collections,bisect
genes=collections.defaultdict(list)
for r in csv.DictReader(open(sys.argv[1]),delimiter='\t'):
    genes[(r['chrom'],r['strand'])].append((int(r['lo']),int(r['hi']),r['gene_id']))
idx={}
for k,v in genes.items():
    v.sort(); idx[k]=( [x[0] for x in v], v )
maxlen=max(h-l for v in genes.values() for l,h,_ in v)
fl={}
if sys.argv[3]!='none':
    for r in csv.DictReader(open(sys.argv[3]),delimiter='\t'):
        if r['sample']==sys.argv[4] and r['tier']==sys.argv[5]:
            fl[r['read']]=(r['n_fl'],r['n_partial'])
out=gzip.open(sys.argv[6],'wt')
nin=nout=0
for line in open(sys.argv[2]):
    f=line.rstrip('\n').split('\t'); nin+=1
    ch,s,e,name,strand=f[0],int(f[1]),int(f[2]),f[3],f[5]
    k=(ch,strand)
    if k not in idx: continue
    sizes=[int(x) for x in f[10].strip(',').split(',')]; starts=[int(x) for x in f[11].strip(',').split(',')]
    blocks=[(s+a,s+a+b) for a,b in zip(starts,sizes)]
    los,v=idx[k]
    i0=bisect.bisect_left(los,s-maxlen); i1=bisect.bisect_right(los,e)
    best=None;bo=0
    for lo,hi,gid in v[i0:i1]:
        if hi<=s or lo>=e: continue
        o=sum(max(0,min(b1,hi)-max(b0,lo)) for b0,b1 in blocks)
        if o>bo: bo=o;best=gid
    if best is None or bo<50: continue
    tss,tes=(s,e-1) if strand=='+' else (e-1,s)
    nf,npart=fl.get(name,('1','0'))
    out.write('\t'.join(map(str,[best,name,tss,tes,len(blocks),bo,nf,npart]))+'\n'); nout+=1
out.close()
print(sys.argv[4],sys.argv[5],nin,nout)
