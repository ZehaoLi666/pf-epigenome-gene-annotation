exec(open("v4/common.py").read())
import sys; sys.path.insert(0,"fig1"); import tracks
S=pd.read_csv("fig1/v71_vs_reference_ends.tsv",sep="\t").set_index('gene_id'); NREF=len(S); assert NREF==2389
c5,c3=S.c5_100,S.c3_100; FIX=['over-long','too short','no annotated UTR','between ends']
nboth=int(((c5=='agree')&(c3=='agree')).sum()); ndiff=int((c5.isin(FIX)|c3.isin(FIX)).sum())
FUN=pd.read_csv("fig2/funnel.tsv",sep="\t")
fig=plt.figure(figsize=(7.2,10.9))
gs=GridSpec(5,6,figure=fig,height_ratios=[2.0,1.0,0.36,1.45,2.55],hspace=0.48,wspace=1.2,left=0.105,right=0.985,top=0.97,bottom=0.025)
# ---------- A data + reference funnel ----------
sA=gs[0,0:4].subgridspec(2,2,height_ratios=[0.42,1.6],width_ratios=[0.34,1],hspace=0.25,wspace=0.0)
axA=fig.add_subplot(sA[0,:]); axA.set_xlim(0,100); axA.set_ylim(0,100); axA.axis('off'); plabel(axA,'A',dx=-0.03,dy=1.05)
axA.text(0,118,"Long-read data and reference construction",fontsize=7,va='bottom')
ins=[("Shaw 2022","PacBio Iso-Seq + ONT\n5′-cap enriched · 3 stages"),("Lee 2021","nanopore direct RNA\nmixed stages"),
     ("Hoshizaki 2022","nanopore\nmixed stages"),("Yang 2021","PacBio Iso-Seq\nmixed stages")]
for i,(t1,t2) in enumerate(ins):
    x=1+i*25; box(axA,x,22,23,76,ec='#888888')
    axA.text(x+11.5,78,t1,fontsize=6.0,fontweight='bold',ha='center',va='center'); axA.text(x+11.5,43,t2,fontsize=4.9,ha='center',va='center',color='#555555',linespacing=1.05)
    arrow(axA,(x+11.5,20),(x+11.5,-6))
FUN=pd.read_csv("fig2/funnel.tsv",sep="\t")
steps=[("All transcript models",10835,316897,None,0),
       ("Supported (Gold + Silver)",int(FUN.loci[0]),int(FUN.models[0]),None,0)]
steps+=[(n,int(l),int(m),None,g) for n,l,m,g in zip(["On protein-coding genes","Detected in a stage","+ structure & stage junctions","+ 3′ end ≥2 reads, not A-primed","Covers the complete CDS"],FUN.loci[1:],FUN.models[1:],[0,1,1,1,2])]
GRP=["Step 1  Long-read evidence (datasets pooled)","Step 2  Stage confirmation","Step 3  Complete gene model"]
axF0=fig.add_subplot(sA[1,1]); yy=[]; yc=0.0; prev=None
for st in steps:
    if st[4]!=prev: yc-=0.55 if prev is not None else 0; yy.append(None); prev=st[4]
    yc-=1; yy.append(yc)
# yy has None markers for group headers; header y = first bar y + 0.78
hdr=[]; ys=[]
for k,v in enumerate(yy):
    if v is None: hdr.append(yy[k+1]+0.72)
    else: ys.append(v)
y=np.array(ys)
cols=['#cfcfcf','#9e9e9e','#9e9e9e']+['#8fa3c7']*3+[COL['ref']]
axF0.barh(y,[s_[1] for s_ in steps],color=cols,height=0.62)
for yi,(nm,l,m,note,g) in zip(y,steps):
    axF0.text(l+150,yi,f"{l:,} {'loci' if nm in ('All transcript models','Supported (Gold + Silver)') else 'genes'} · {m:,} models",va='center',fontsize=5.1)
    if note: axF0.text(l+150,yi-0.42,note,va='center',fontsize=4.6,color='#666666')
for hy,t in zip(hdr,GRP): axF0.text(-0.55,hy,t,transform=axF0.get_yaxis_transform(),fontsize=5.3,fontweight='bold',va='center',ha='left',color='#333333')
axF0.set_yticks(y); axF0.set_yticklabels([s_[0] for s_ in steps],fontsize=5.3); axF0.set_xlim(0,17000); axF0.set_xticks([0,5000,10000]); axF0.spines['bottom'].set_bounds(0,10000)
axF0.set_ylim(y.min()-0.6,hdr[0]+0.5)
axF0.set_xlabel("Loci (first two rows) or genes",fontsize=6); axF0.tick_params(axis='y',length=0)
# ---------- B genome-wide long-read coverage (read level) ----------
GC=pd.read_csv("cov/genome_coverage.tsv",sep="\t"); GN=pd.read_csv("cov/gene_coverage_summary.tsv",sep="\t")
sB=gs[0,4:6].subgridspec(2,1,height_ratios=[1,1],hspace=0.95)
axB=fig.add_subplot(sB[0]); plabel(axB,'B',dx=-0.42,dy=1.22)
lab={'ALL':'All pooled','Shaw2022_PacBio':'Shaw 2022 PacBio','Shaw2022_ONT':'Shaw 2022 ONT','Lee2021':'Lee 2021','Hoshizaki2022':'Hoshizaki 2022','Yang2021':'Yang 2021'}
order=['ALL','Shaw2022_PacBio','Shaw2022_ONT','Lee2021','Hoshizaki2022','Yang2021']; yb=np.arange(len(order))[::-1]
g1=GC[GC.min_depth==1].set_index('source').genome; g5=GC[GC.min_depth==5].set_index('source').genome
axB.barh(yb,[g1[o]*100 for o in order],color='#c6d3e8',height=0.66)
axB.barh(yb,[g5[o]*100 for o in order],color=COL['ref'],height=0.66)
for yi,o in zip(yb,order): axB.text(g1[o]*100+1.5,yi,f"{g1[o]:.0%}",va='center',fontsize=4.8)
axB.set_yticks(yb); axB.set_yticklabels([lab[o] for o in order],fontsize=5.0); axB.tick_params(axis='y',length=0)
axB.get_yticklabels()[0].set_fontweight('bold')
axB.set_xlim(0,112); axB.set_xticks([0,50,100]); axB.spines['bottom'].set_bounds(0,100); axB.set_xlabel("Nuclear genome covered (%)",fontsize=5.6,labelpad=1)
axB.text(0.99,1.02,"≥1 read",transform=axB.transAxes,color='#7f96bd',fontsize=4.8,ha='right',va='bottom'); axB.text(0.70,1.02,"≥5 reads",transform=axB.transAxes,color=COL['ref'],fontsize=4.8,ha='right',va='bottom')
title(axB,"Long reads cover most of the genome",y=1.13)
axB2=fig.add_subplot(sB[1])
PG=pd.read_csv("cov/gene_longread_coverage.tsv",sep="\t",index_col=0)
grp=[('confirmed','Reference',COL['ref']),('truncated-only','Only truncated models','#8c8c8c'),('HC model, unconfirmed','Not stage-confirmed','#a9a9a9'),('no HC model','No supported model','#cfcfcf')]
yg=np.arange(len(grp))[::-1]
for yi,(k,l,cc) in zip(yg,grp):
    v=np.log10(PG.loc[PG.lr_group==k,'ALL_reads']+1)
    bp=axB2.boxplot([v],positions=[yi],orientation='horizontal',widths=0.58,showfliers=False,patch_artist=True,medianprops=dict(color='white' if cc==COL['ref'] else INK,lw=0.9),
                    boxprops=dict(fc=cc,ec=INK,lw=0.5),whiskerprops=dict(lw=0.5),capprops=dict(lw=0.5))
    axB2.text(4.62,yi,f"{int(10**np.median(v)-0.5):,}",va='center',fontsize=4.9,ha='right')
axB2.set_yticks(yg); axB2.set_yticklabels([f"{l} ({(PG.lr_group==k).sum():,})" for k,l,_ in grp],fontsize=5.0); axB2.tick_params(axis='y',length=0)
axB2.set_xlim(0,4.7); axB2.set_xticks([0,1,2,3,4]); axB2.set_xticklabels(['1','10','100','1k','10k']); axB2.spines['bottom'].set_bounds(0,4)
axB2.set_xlabel("Long reads on the CDS per gene (pooled)",fontsize=5.6,labelpad=1)
title(axB2,"Non-reference genes also have reads",y=1.02)
# ---------- C-E long-read limitations ----------
CM=pd.read_csv("fig2/model_completeness.tsv",sep="\t").set_index('cls'); tot=int(CM.n.sum())
axC=fig.add_subplot(gs[1,0:2]); plabel(axC,'C',dx=-0.42,dy=1.0)
rowsC=[('full CDS','Complete CDS',COL['ref']),('5′ truncated','5′ end inside CDS',COL['trunc']),('3′ truncated','3′ end inside CDS','#c8c8c8')]
yc_=np.arange(3)[::-1]
for yi,(k,l,cc) in zip(yc_,rowsC):
    n=int(CM.loc[k,'n']); axC.barh(yi,n/tot*100,color=cc,height=0.6)
    axC.text(n/tot*100+1,yi,f"{n/tot:.0%}",va='center',fontsize=5.3)
axC.set_yticks(yc_); axC.set_yticklabels([r[1] for r in rowsC],fontsize=5.1); axC.tick_params(axis='y',length=0)
axC.set_xlim(0,55); axC.set_xticks([0,25,50]); axC.set_xlabel(f"Stage-confirmed mRNA models (%, n = {tot:,})",fontsize=5.6,labelpad=1)
title(axC,"Most models are truncated",y=1.02)
TG=pd.read_csv("fig2/tss_capG_groups.tsv",sep="\t").set_index('grp').loc[['Upstream, CAGE ≥ q25','Upstream, CAGE < q25','Internal, other','Internal-weak']]
axD=fig.add_subplot(gs[1,2:4]); plabel(axD,'D',dx=-0.30,dy=1.0)
xx=np.arange(4); axD.bar(xx,TG.ONT_extraG_frac,color=[COL['tss'],COL['tss'],COL['art'],COL['art']],width=0.62)
for xi,v in enumerate(TG.ONT_extraG_frac): axD.text(xi,v+0.015,f"{v:.0%}",ha='center',fontsize=5.3)
axD.set_xticks(xx); axD.set_xticklabels(['Upstream,\nCAGE high','Upstream,\nCAGE low','Inside\nCDS','Inside CDS,\nweak'],fontsize=4.8)
axD.set_ylim(0,0.68); axD.set_yticks([0,0.2,0.4,0.6]); axD.set_yticklabels(['0','20%','40%','60%']); axD.set_ylabel("ONT reads with a\nnon-templated 5′ G",fontsize=5.6)
title(axD,"Degraded 5′ ends lack the cap",y=1.02)
IP=pd.read_csv("cov/internal_priming_v71.tsv",sep="\t")
axE=fig.add_subplot(gs[1,4:6]); plabel(axE,'E',dx=-0.42,dy=1.0)
ye=np.arange(len(IP))[::-1]
axE.barh(ye,IP.frac*100,color=[COL['art'],COL['art'],'#e8b27a',COL['utr'],'#cfcfcf'],height=0.62)
for yi,v in zip(ye,IP.frac): axE.text(v*100+0.6,yi,f"{v:.0%}",va='center',fontsize=5.3)
axE.set_yticks(ye); axE.set_yticklabels(["Long reads, inside CDS","Long reads, beyond CDS","All long-read ends","PlasmoDB v71","Random sites"],fontsize=5.0); axE.tick_params(axis='y',length=0)
axE.set_xlim(0,34); axE.set_xticks([0,10,20,30]); axE.set_xlabel("3′ ends with ≥80% A in the next 20 nt (%)",fontsize=5.6,labelpad=1)
title(axE,"False 3′ ends from internal priming",y=1.02)
# ---------- C comparison rules (horizontal) ----------
axR=fig.add_subplot(gs[2,:]); axR.set_xlim(0,100); axR.set_ylim(0,23); axR.axis('off'); plabel(axR,'F',dx=-0.06,dy=0.95)
axR.text(0,24.5,"Each v71 UTR end of the 2,389 reference genes was compared with every reference end cluster (main and alternative, ±100 bp); 5′ end shown",fontsize=6.3,va='bottom')
cases=[("Agrees",'agree',3,[(3.2,14),(5.0,14),(3.1,13.8)],None),("Over-long",'ol_none',1,[(4.0,14),(5.2,14),(4.1,13.8)],(1,4.0)),
       ("Too short",'short',4.5,[(2.5,14),(2.8,14),(2.6,13.8)],(2.5,4.5)),("No annotated UTR",'noutr',6,[(3.5,14),(3.8,14),(3.6,13.8)],(3.5,6)),
       ("Between ends",'between',3.6,[(2.2,14),(5.0,14),(2.4,13.8)],None)]
for i,(nm,ck,u5,reads,shade) in enumerate(cases):
    x0=9+i*18.4
    axR.text(x0+8.5,19.5,nm,fontsize=6.3,fontweight='bold',ha='center',color=COL[ck] if ck!='ol_none' else '#c0640f')
    yy=13.0; axR.plot([x0,x0+17],[yy,yy],color='#bbbbbb',lw=0.5)
    if u5<6: axR.add_patch(Rectangle((x0+u5,yy-1.1),6-u5,2.2,fc=COL['utr'],ec='none'))
    axR.add_patch(Rectangle((x0+6,yy-2.0),6,4.0,fc=COL['cds'],ec='none')); axR.add_patch(Rectangle((x0+12,yy-1.1),2,2.2,fc=COL['utr'],ec='none'))
    for kk,(p,q) in enumerate(reads): axR.plot([x0+p,x0+q],[7.5-kk*1.9]*2,color='#8c8c8c',lw=1.0,solid_capstyle='butt')
    if shade: axR.add_patch(Rectangle((x0+shade[0],2.0),shade[1]-shade[0],13.5,fc=COL[ck],alpha=0.4,ec='none',zorder=0))
axR.text(0.5,13.0,"PlasmoDB v71",fontsize=5.5,va='center'); axR.text(0.5,5.6,"reference\nend clusters",fontsize=5.5,va='center',linespacing=0.9)
# ---------- D heatmap ----------
axG=fig.add_subplot(gs[3,0:2]); plabel(axG,'G',dx=-0.30)
R5=['agree','over-long','too short','no annotated UTR','between ends']; C3=R5+['not assessable']
v=pd.crosstab(pd.Categorical(c5,R5),pd.Categorical(c3,C3)).values
axG.imshow(np.log10(v+1),cmap='Greys',vmin=0,vmax=np.log10(v.max()+1)*1.25,aspect='auto')
for i in range(v.shape[0]):
    for j in range(v.shape[1]): axG.text(j,i,f"{v[i,j]:,}",ha='center',va='center',fontsize=5.2,color='white' if v[i,j]>250 else INK)
lab=['Agrees','Over-\nlong','Too\nshort','No\nUTR','Between\nends','Not\nassessed']
axG.set_xticks(range(6)); axG.set_xticklabels(lab,fontsize=4.8); axG.set_yticks(range(5)); axG.set_yticklabels(lab[:5],fontsize=5)
axG.set_xlabel("3′ end vs reference"); axG.set_ylabel("5′ end vs reference")
for k in ['top','right']: axG.spines[k].set_visible(True)
axG.add_patch(Rectangle((-0.5,-0.5),1,1,fill=False,ec=COL['agree'],lw=1.2))
title(axG,f"{ndiff:,} of {NREF:,} reference genes ({ndiff/NREF:.0%})\ndiffer from v71 at ≥1 end")
# ---------- H, I UTR length scatter ----------
maxd=lambda q: q.apply(lambda x: max(map(float,str(x).split(';'))) if isinstance(x,str) and x else np.nan)
gcol={'agree':COL['agree'],'over-long':'#d9822b','too short':COL['short'],'no annotated UTR':COL['noutr'],'between ends':COL['between']}
tr=lambda x: np.log10(1+np.clip(np.asarray(x,float),0,None)); ticks=[0,100,300,1000,3000]
for k,(Lv,Dm,cc,lbl) in enumerate([(S.L5_v71,maxd(S.ref_TSS_dists),c5,"5′ UTR"),(S.L3_v71,maxd(S.ref_TTS_dists),c3,"3′ UTR")]):
    ax=fig.add_subplot(gs[3,2+2*k:4+2*k]); plabel(ax,'HI'[k],dx=-0.24 if k==0 else -0.10)
    for gname in ['agree','between ends','no annotated UTR','too short','over-long']:
        m=(cc==gname)&Dm.notna()
        jit=np.where(Lv[m]==0,np.random.default_rng(1).uniform(-0.06,0.06,m.sum()),0)
        ax.scatter(tr(Lv[m])+jit,tr(Dm[m]),s=1.3,lw=0,color=gcol[gname],alpha=0.7,rasterized=True)
    lim=(-0.15,np.log10(6001)); ax.plot(lim,lim,color='#777777',lw=0.5,ls=(0,(2,2)),zorder=0)
    ax.set_xlim(lim); ax.set_ylim(lim); tl=['0','100','300','1k','3k']
    ax.set_xticks(tr(ticks)); ax.set_xticklabels(tl); ax.set_yticks(tr(ticks)); ax.set_yticklabels(tl if k==0 else [])
    ax.set_xlabel(f"PlasmoDB v71 {lbl} length (bp)",fontsize=5.8)
    if k==0: ax.set_ylabel("Farthest reference end (bp)",fontsize=5.8,labelpad=1)
    title(ax,lbl,y=1.02)
    if k==1:
        for i,g_ in enumerate(['agree','over-long','too short','no annotated UTR','between ends']):
            ax.text(0.98,0.40-i*0.08,g_,transform=ax.transAxes,ha='right',fontsize=5.0,color=gcol[g_])
# ---------- F, G tracks ----------
gff=tracks.load_gff()
LR=tracks.read_bed12(f"{tracks.IGV}/confirmed_transcripts_by_category.bed"); LR['locus']=LR.name.str.split('|').str[2]; LR['cat']=LR.name.str.split('|').str[1]
LNC=tracks.read_bed12(f"{tracks.IGV}/lncRNA_catalogue_3datasets.bed"); RB=pd.read_csv(BOUND,sep='\t').set_index('gene_id')
EX=[dict(gene="PF3D7_0424600",win=(1112950,1117900),title="5′ UTR absorbs an upstream lncRNA",mode='ol5'),
    dict(gene="PF3D7_0627500",win=(1100350,1103250),title="3′ UTR extends past the long-read 3′ ends",mode='ol3'),
    dict(gene="PF3D7_1334800",win=(1412200,1415300),title="Missing 5′ UTR defined by long reads",mode='new5')]
sub=gs[4,:].subgridspec(5,4,height_ratios=[0.8,1.5,0.75,1.2,0.85],width_ratios=[0.13,1,1,1],hspace=0.14,wspace=0.1)
rowlab=["PlasmoDB v71","Long-read\nmodels","lncRNA","Stranded RNA\n(DAFT-seq)","CAGE"]
for j,ex in enumerate(EX):
    gname=ex['gene']; a,b=ex['win']; r=gff[gff.gene==gname]; chrom=r.chrom.iloc[0]; s=r.strand.iloc[0]; os_='+' if s=='-' else '-'
    axs=[fig.add_subplot(sub[i,j+1]) for i in range(5)]; rb=RB.loc[gname]; sg_=1 if s=='+' else -1
    if ex['mode']=='ol5': m_=max(map(float,S.loc[gname,'ref_TSS_dists'].split(';'))); shade=tuple(sorted((rb.annot_tss, rb.atg-sg_*m_))); shc=COL['ol_unit']; al=0.22
    elif ex['mode']=='ol3': m_=max(map(float,S.loc[gname,'ref_TTS_dists'].split(';'))); shade=tuple(sorted((rb.annot_tts, rb.stop+sg_*m_))); shc=COL['ol_anti']; al=0.45
    else: m_=max(map(float,S.loc[gname,'ref_TSS_dists'].split(';'))); shade=tuple(sorted((rb.atg, rb.atg-sg_*m_))); shc=COL['noutr']; al=0.22
    for ax in axs:
        ax.axvspan(*shade,color=shc,alpha=al,lw=0,zorder=0); ax.set_xlim((b,a) if s=='-' else (a,b)); ax.set_xticks([]); ax.set_yticks([])
        for sp in ['left','bottom','top','right']: ax.spines[sp].set_visible(False)
    ax=axs[0]; ax.set_ylim(-1.0,1.5)
    gg=gff[(gff.chrom==chrom)&(gff.end>a)&(gff.start<b)&gff.type.isin(['CDS','five_prime_UTR','three_prime_UTR','exon'])]
    for gn_,sb in gg.groupby('gene'):
        foc=gn_==gname; alp=1 if foc else 0.4; ex_=sb[sb.type=='exon']
        if len(ex_)>1: ax.plot([ex_.start.min(),ex_.end.max()],[0,0],color=COL['cds'],lw=0.5,alpha=alp)
        for _,f in sb[sb.type!='exon'].iterrows():
            h=0.9 if f.type=='CDS' else 0.5; ax.add_patch(Rectangle((f.start,-h/2),f.end-f.start,h,fc=COL['cds'] if f.type=='CDS' else COL['utr'],ec='none',alpha=alp))
        if foc: ax.text((max(sb.start.min(),a)+min(sb.end.max(),b))/2,0.85,gn_,ha='center',fontsize=5.4,fontweight='bold')
    ax=axs[1]; lr=LR[(LR.chrom==chrom)&(LR.end>a)&(LR.start<b)]; foc=lr[(lr.locus==gname)&~lr.cat.str.contains('trunc')].head(5)
    for k,(_,t_) in enumerate(foc.iterrows()):
        ax.plot([t_.start,t_.end],[-k,-k],color='#555555',lw=0.4)
        for bb in t_.blocks: ax.plot(bb,[-k,-k],color='#555555',lw=2.0,solid_capstyle='butt')
    ax.set_ylim(-max(len(foc),5)+0.3,0.8)
    ax=axs[2]; units=[(t_.start,t_.end,t_.strand) for _,t_ in lr[(lr.locus!=gname)&lr.cat.str.contains('lncRNA|ncRNA|antisense|no_v68')].iterrows()]
    units+=[(t_.start,t_.end,t_.strand) for _,t_ in LNC[(LNC.chrom==chrom)&(LNC.end>a)&(LNC.start<b)].iterrows()]
    for k,(u0,u1,us) in enumerate(units): ax.add_patch(Rectangle((u0,-k*0.9-0.3),u1-u0,0.6,fc=COL['lnc'] if us==s else 'none',ec=COL['lnc'],lw=0.5))
    ax.set_ylim(-max(len(units),2)*0.9-0.1,0.6)
    ax=axs[3]; nb=int((b-a)/10); x=np.linspace(a,b,nb)
    ds=np.log10(1+tracks.daft(chrom,a,b,s,nb)); do=np.log10(1+tracks.daft(chrom,a,b,os_,nb))
    ax.fill_between(x,0,ds,color=COL['daft'],lw=0); ax.fill_between(x,0,-do,color='#b8c7e0',lw=0); ax.axhline(0,color='#999999',lw=0.3)
    ym=max(ds.max(),do.max())*1.05; ax.set_ylim(-ym,ym)
    ax.text(0.01,0.93,"gene strand",transform=ax.transAxes,fontsize=4.8,va='top',color=COL['daft']); ax.text(0.01,0.07,"opposite",transform=ax.transAxes,fontsize=4.8,va='bottom',color='#8aa0c8')
    ax=axs[4]; cg=tracks.cage(chrom,a,b,s,nb); ax.fill_between(x,0,cg,color=COL['cage'],lw=0.6); ax.plot(x,cg,color=COL['cage'],lw=0.5); ax.set_ylim(0,cg.max()*1.1); ax.axhline(0,color='#999999',lw=0.3)
    ax.text(0.99,0.9,f"max {cg.max():.0f} CPM",transform=ax.transAxes,fontsize=4.8,ha='right',va='top',color=COL['cage'])
    sb0=b-150 if s=='-' else a+150; sb1=sb0-500 if s=='-' else sb0+500
    ax.plot([sb0,sb1],[-0.18*cg.max()]*2,color=INK,lw=0.8,clip_on=False); ax.text((sb0+sb1)/2,-0.3*cg.max(),"500 bp",fontsize=5,ha='center',va='top',clip_on=False)
    ax.text(0.99,-0.2,f"{chrom.replace('Pf3D7_','chr').replace('_v3','')}:{a:,}–{b:,} ({s})",transform=ax.transAxes,fontsize=4.8,ha='right',va='top',color='#666666')
    axs[0].text(0,1.45,ex['title'],transform=axs[0].transAxes,fontsize=6.5,va='bottom'); axs[0].text(-0.02,1.47,"JKL"[j],transform=axs[0].transAxes,fontsize=9,fontweight='bold',ha='right',va='bottom')
    if j==0:
        for i,l in enumerate(rowlab): axs[i].text(-0.04,0.5,l,transform=axs[i].transAxes,fontsize=5.4,ha='right',va='center',linespacing=0.95)
fig.savefig("v4/Figure1.png",dpi=300); fig.savefig("v4/Figure1.pdf",dpi=300); OV=overlaps(fig)
