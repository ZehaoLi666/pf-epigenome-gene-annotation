"""
Figure 1 - data preparation.

Gathers every input of Figure 1 and writes one set of source tables per panel into ../data/.
The panel scripts (panel_*.py) read ONLY those tables.

Inputs
  * intermediate result tables (RESULTS_DIR in config/paths.yaml; sha256-checked against the versions used for the paper)
  * local track files under LAB_DIR (PlasmoDB GFF, long-read BED, DAFT-seq / CAGE bigWigs)
This is the only script that needs non-standard packages (pyBigWig) or external paths.
The derivations are copied from the generating scripts kept in code/original/ (make_fig1.py, tracks.py, common.py).
Run:  python prep_data.py      (env with pandas, numpy, pyBigWig)
"""
from pfconfig import LAB_DIR, result  # paths: config/paths.yaml
import os, hashlib, numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data'); os.makedirs(DATA, exist_ok=True)

# ---------------------------------------------------------------- sources
SOURCES = {
    'v71_vs_reference_ends.tsv': dict(sha256='e7fc980da1114f2edf6d4527663f663319f524f193a8702c65eb644e9530940a',
        path=result('v71_vs_reference_ends.tsv')),
    'funnel.tsv': dict(sha256='fca5030fa5c82b120d16e2b7b2c7517f6738d5cb44888dc5dfa741dcc018f132',
        path=result('funnel.tsv')),
    'genome_coverage.tsv': dict(sha256='825ddcc1c4a921d38d9e13d9045f62ce7c083aba6c5d6c48d35c1314d648b870',
        path=result('genome_coverage.tsv')),
    'gene_longread_coverage.tsv.gz': dict(sha256='23ae9a964952b74dd19cf4f7658cf9383ed3d68d19c8c97ce9254f948d16a7e6',
        path=result('gene_longread_coverage.tsv.gz')),
    'model_completeness.tsv': dict(sha256='5f0cb9d3bf9ddc9f6c7e73839db5217424e8d1bc9f7afe893bfdfa5a924fc7e9',
        path=result('model_completeness.tsv')),
    'tss_capG_groups.tsv': dict(sha256='1b41d2a28f23e76ab5db655d582bbff4afe28a03077448759c572bff055b4d92',
        path=result('tss_read_evidence_by_group.tsv')),
    'internal_priming_v71.tsv': dict(sha256='96b04fa54dbdb7c1716fd52040b3cc528701509ae35de9a8a20fc4324f22c75a',
        path=result('internal_priming_v71.tsv')),
    'gene_boundaries_v3_all_genes.tsv': dict(sha256='730d59b8a2c05228533b3ca2dda38b8fdd8468ceb922f1782102ec14905274ef',
        path=result('gene_boundaries_v3_all_genes.tsv')),
}

LAB = LAB_DIR
IGV = f"{LAB}/IGV_longread_transcript_models"
GFF = f"{LAB}/PlasmoDB/PlasmoDB-68/PlasmoDB-68_Pfalciparum3D7.gff"   # file used by the approved script for the gene-model track
CAGE = f"{LAB}/lncRNA project tracks/cage_tss"

def src(name):
    """Path of an input table (RESULTS_DIR/<name>); checks sha256."""
    s = SOURCES[name]
    p = s['path']
    h = hashlib.sha256(open(p, 'rb').read()).hexdigest()
    assert h == s['sha256'], f"checksum mismatch for {name}"
    return p

def out(df, name):
    df.to_csv(os.path.join(DATA, name), sep='\t', index=False)
    print(f"wrote data/{name}  ({len(df)} rows)")

# ---------------------------------------------------------------- track helpers (fig1/tracks.py)
def load_gff():
    g = pd.read_csv(GFF, sep="\t", comment="#", header=None, names="chrom src type start end score strand phase attr".split(), dtype={'start': int, 'end': int})
    g = g[g.type.isin(['protein_coding_gene', 'ncRNA_gene', 'mRNA', 'CDS', 'five_prime_UTR', 'three_prime_UTR', 'exon', 'gene'])].copy()
    g['ID'] = g.attr.str.extract(r'ID=([^;]+)')
    g['Parent'] = g.attr.str.extract(r'Parent=([^;]+)')
    g['gene'] = np.where(g.type.str.endswith('gene'), g.ID, None)
    g.loc[g.type == 'mRNA', 'gene'] = g.loc[g.type == 'mRNA', 'Parent']
    tx2g = dict(zip(g[g.type == 'mRNA'].ID, g[g.type == 'mRNA'].Parent))
    m = ~g.type.str.endswith('gene') & (g.type != 'mRNA')
    g.loc[m, 'gene'] = g.loc[m, 'Parent'].map(lambda p: tx2g.get(str(p).split(',')[0], p))
    return g

def read_bed12(path):
    rows = []
    for ln in open(path):
        if ln.startswith(('track', '#')): continue
        f = ln.rstrip('\n').split('\t')
        r = dict(chrom=f[0], start=int(f[1]), end=int(f[2]), name=f[3], strand=f[5] if len(f) > 5 else '.')
        if len(f) >= 12:
            r['rgb'] = f[8]; sz = [int(x) for x in f[10].strip(',').split(',')]; st = [int(x) for x in f[11].strip(',').split(',')]
            r['blocks'] = [(r['start'] + s, r['start'] + s + z) for s, z in zip(st, sz)]
        else: r['blocks'] = [(r['start'], r['end'])]
        rows.append(r)
    return pd.DataFrame(rows)

def bw_vals(path, chrom, a, b, nb):
    import pyBigWig
    bw = pyBigWig.open(path); v = np.array(bw.stats(chrom, a, b, type='mean', nBins=nb), dtype=float); bw.close()
    return np.nan_to_num(v)

def daft(chrom, a, b, strand, nb):
    return bw_vals(f"{IGV}/DAFTseq_pool_{'plus' if strand == '+' else 'minus'}_CPM_10bp.bw", chrom, a, b, nb)

def cage(chrom, a, b, strand, nb):
    # NB: as in the approved script, the three stage tracks are SUMMED (ring + trophozoite + schizont CPM)
    s = 'pos' if strand == '+' else 'neg'
    return sum(bw_vals(f"{CAGE}/CTSS_{st}_{s}_CPM.bw", chrom, a, b, nb) for st in ['ring', 'trophozoite', 'schizont'])

# ================================================================ panel A
FUN = pd.read_csv(src('funnel.tsv'), sep="\t")
GRP = ["Step 1  Long-read evidence (datasets pooled)", "Step 2  Stage confirmation", "Step 3  Complete gene model"]
rows = [dict(step=1, step_label=GRP[0], row_label="All transcript models", n_loci_or_genes=10835, unit='loci', n_models=316897,
             source="end-aware merge of TAMA models (constant in approved script; see legend A)"),
        dict(step=1, step_label=GRP[0], row_label="Supported (Gold + Silver)", n_loci_or_genes=int(FUN.loci[0]), unit='loci', n_models=int(FUN.models[0]), source=f"funnel.tsv: {FUN.step[0]}")]
for nm, l, m, g, st in zip(["On protein-coding genes", "Detected in a stage", "+ structure & stage junctions", "+ 3′ end ≥2 reads, not A-primed", "Covers the complete CDS"],
                           FUN.loci[1:], FUN.models[1:], [1, 2, 2, 2, 3], FUN.step[1:]):
    rows.append(dict(step=g, step_label=GRP[g - 1], row_label=nm, n_loci_or_genes=int(l), unit='genes', n_models=int(m), source=f"funnel.tsv: {st}"))
out(pd.DataFrame(rows), "Figure1_A_reference_funnel.tsv")
out(pd.DataFrame(dict(dataset=["Shaw 2022", "Lee 2021", "Hoshizaki 2022", "Yang 2021"],
                      line1=["PacBio Iso-Seq + ONT", "nanopore direct RNA", "nanopore", "PacBio Iso-Seq"],
                      line2=["5′-cap enriched · 3 stages", "mixed stages", "mixed stages", "mixed stages"])), "Figure1_A_datasets.tsv")

# ================================================================ panel B
GC = pd.read_csv(src('genome_coverage.tsv'), sep="\t")
lab = {'ALL': 'All pooled', 'Shaw2022_PacBio': 'Shaw 2022 PacBio', 'Shaw2022_ONT': 'Shaw 2022 ONT', 'Lee2021': 'Lee 2021', 'Hoshizaki2022': 'Hoshizaki 2022', 'Yang2021': 'Yang 2021'}
GC.insert(1, 'label', GC.source.map(lab))
out(GC, "Figure1_B_genome_coverage.tsv")
PG = pd.read_csv(src('gene_longread_coverage.tsv.gz'), sep="\t", index_col=0)
out(PG[['lr_group', 'ALL_reads']].reset_index().rename(columns={'gene': 'gene_id', 'ALL_reads': 'pooled_reads_on_CDS'}), "Figure1_B_gene_longread_counts.tsv")

# ================================================================ panel C
CM = pd.read_csv(src('model_completeness.tsv'), sep="\t")
CM['label'] = CM.cls.map({'full CDS': 'Complete CDS', '5′ truncated': '5′ end inside CDS', '3′ truncated': '3′ end inside CDS'})
CM['percent'] = CM.n / CM.n.sum() * 100
out(CM.rename(columns={'n': 'n_models', 'loci': 'n_genes_total'}), "Figure1_C_model_completeness.tsv")

# ================================================================ panel D
TG = pd.read_csv(src('tss_capG_groups.tsv'), sep="\t").set_index('grp').loc[['Upstream, CAGE ≥ q25', 'Upstream, CAGE < q25', 'Internal, other', 'Internal-weak']].reset_index()
TG.insert(1, 'label', ['Upstream, CAGE high', 'Upstream, CAGE low', 'Inside CDS', 'Inside CDS, weak'])
out(TG, "Figure1_D_tss_capG.tsv")

# ================================================================ panel E
IP = pd.read_csv(src('internal_priming_v71.tsv'), sep="\t")
IP['label_source'] = IP.label.str.replace('\n', ' ')
IP['label'] = ["Long reads, inside CDS", "Long reads, beyond CDS", "All long-read ends", "PlasmoDB v71", "Random sites"]
out(IP[['label', 'label_source', 'frac', 'n']], "Figure1_E_internal_priming.tsv")

# ================================================================ panels G, H, I
S = pd.read_csv(src('v71_vs_reference_ends.tsv'), sep="\t").set_index('gene_id'); assert len(S) == 2389
out(S[['c5_100', 'c3_100', 'c5_200', 'c3_200']].reset_index(), "Figure1_G_gene_classes.tsv")
R5 = ['agree', 'over-long', 'too short', 'no annotated UTR', 'between ends']; C3 = R5 + ['not assessable']
M = pd.crosstab(pd.Categorical(S.c5_100, R5), pd.Categorical(S.c3_100, C3))
M.index.name = 'class_5prime (rows) / class_3prime (columns)'; out(M.reset_index(), "Figure1_G_class_matrix.tsv")
maxd = lambda q: q.apply(lambda x: max(map(float, str(x).split(';'))) if isinstance(x, str) and x else np.nan)
for k, (L, D, c, nm) in enumerate([(S.L5_v71, S.ref_TSS_dists, S.c5_100, 'H_5UTR'), (S.L3_v71, S.ref_TTS_dists, S.c3_100, 'I_3UTR')]):
    out(pd.DataFrame({'gene_id': S.index, 'class': c.values, 'v71_UTR_length_bp': L.values,
                      'farthest_reference_end_bp': maxd(D).values, 'all_reference_end_dists_bp': D.values}), f"Figure1_{nm}_lengths.tsv")

# ================================================================ panels J, K, L (genome tracks)
COLSH = dict(ol_unit='#a6400a', ol_anti='#f5b041', noutr='#2a9d8f')
gff = load_gff()
LR = read_bed12(f"{IGV}/confirmed_transcripts_by_category.bed"); LR['locus'] = LR.name.str.split('|').str[2]; LR['cat'] = LR.name.str.split('|').str[1]
LNC = read_bed12(f"{IGV}/lncRNA_catalogue_3datasets.bed"); RB = pd.read_csv(src('gene_boundaries_v3_all_genes.tsv'), sep='\t').set_index('gene_id')
EX = [dict(panel='J', gene="PF3D7_0424600", win=(1112950, 1117900), title="5′ UTR absorbs an upstream lncRNA", mode='ol5'),
      dict(panel='K', gene="PF3D7_0627500", win=(1100350, 1103250), title="3′ UTR extends past the long-read 3′ ends", mode='ol3'),
      dict(panel='L', gene="PF3D7_1334800", win=(1412200, 1415300), title="Missing 5′ UTR defined by long reads", mode='new5')]
for ex in EX:
    P_ = ex['panel']; gname = ex['gene']; a, b = ex['win']; r = gff[gff.gene == gname]; chrom = r.chrom.iloc[0]; s = r.strand.iloc[0]; os_ = '+' if s == '-' else '-'
    rb = RB.loc[gname]; sg_ = 1 if s == '+' else -1
    if ex['mode'] == 'ol5': m_ = max(map(float, S.loc[gname, 'ref_TSS_dists'].split(';'))); shade = tuple(sorted((rb.annot_tss, rb.atg - sg_ * m_))); shc = 'ol_unit'; al = 0.22
    elif ex['mode'] == 'ol3': m_ = max(map(float, S.loc[gname, 'ref_TTS_dists'].split(';'))); shade = tuple(sorted((rb.annot_tts, rb.stop + sg_ * m_))); shc = 'ol_anti'; al = 0.45
    else: m_ = max(map(float, S.loc[gname, 'ref_TSS_dists'].split(';'))); shade = tuple(sorted((rb.atg, rb.atg - sg_ * m_))); shc = 'noutr'; al = 0.22
    out(pd.DataFrame([dict(gene_id=gname, chrom=chrom, chrom_label=chrom.replace('Pf3D7_', 'chr').replace('_v3', ''), window_start=a, window_end=b, strand=s,
                           mode=ex['mode'], title=ex['title'], shade_start=float(shade[0]), shade_end=float(shade[1]), shade_color=COLSH[shc], shade_alpha=al,
                           farthest_reference_end_bp=m_)]), f"Figure1_{P_}_meta.tsv")
    gg = gff[(gff.chrom == chrom) & (gff.end > a) & (gff.start < b) & gff.type.isin(['CDS', 'five_prime_UTR', 'three_prime_UTR', 'exon'])]
    out(gg[['gene', 'type', 'start', 'end', 'strand']].rename(columns={'gene': 'gene_id'}), f"Figure1_{P_}_v71_gene_models.tsv")
    lr = LR[(LR.chrom == chrom) & (LR.end > a) & (LR.start < b)]
    foc = lr[(lr.locus == gname) & ~lr.cat.str.contains('trunc')].head(5)
    out(pd.DataFrame(dict(row=range(len(foc)), name=foc.name.values, start=foc.start.values, end=foc.end.values,
                          blocks=[','.join(f"{x}-{y}" for x, y in t) for t in foc.blocks])), f"Figure1_{P_}_longread_models.tsv")
    units = [(t_.start, t_.end, t_.strand, 'long-read lncRNA unit: ' + t_['name']) for _, t_ in lr[(lr.locus != gname) & lr.cat.str.contains('lncRNA|ncRNA|antisense|no_v68')].iterrows()]
    units += [(t_.start, t_.end, t_.strand, 'lncRNA catalogue: ' + t_['name']) for _, t_ in LNC[(LNC.chrom == chrom) & (LNC.end > a) & (LNC.start < b)].iterrows()]
    out(pd.DataFrame(units, columns=['start', 'end', 'strand', 'source']), f"Figure1_{P_}_lncRNA.tsv")
    nb = int((b - a) / 10); x = np.linspace(a, b, nb)
    out(pd.DataFrame(dict(position=x, daft_gene_strand_cpm=daft(chrom, a, b, s, nb), daft_opposite_strand_cpm=daft(chrom, a, b, os_, nb),
                          cage_ctss_cpm_sum_3_stages=cage(chrom, a, b, s, nb))), f"Figure1_{P_}_coverage_10bp_bins.tsv")
print("done")
