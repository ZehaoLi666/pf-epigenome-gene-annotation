import pandas as pd
import os
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'common.py')).read())
TV = 'to verify'
COLS = ['dataset_id', 'category', 'study', 'year', 'journal', 'PMID', 'DOI', 'accession', 'runs_or_samples_used', 'assay', 'library_notes',
        'stages_timepoints_used', 'strain', 'replicates_used', 'input_control', 'processing_summary', 'role', 'used_in', 'verification_notes']
R = []
def add(**k): R.append({c: k.get(c, '') for c in COLS})

# ---------------- long-read transcriptomes (training labels / reference) ----------------
add(dataset_id='LR1', category='Long-read RNA', study='Shaw et al.', year=2022, journal='PLoS One 17:e0276956', PMID=TV, DOI='10.1371/journal.pone.0276956',
    accession='GEO GSE186109', runs_or_samples_used='30 SRA runs (SRR164801xx); 7 PacBio Sequel samples + ONT SRR18188665, SRR18188666',
    assay='PacBio Iso-Seq (cDNA) and ONT cDNA; 5′-cap-enriched (HseIF4E–eIF4G)', library_notes='PacBio isoform tiers HQ/LQ per sample; ONT reads re-oriented',
    stages_timepoints_used='PacBio: ring, trophozoite, schizont (total and nuclear RNA each) + schizont/trophozoite cytoplasmic; ONT: mixed asexual stages',
    strain='3D7 (to verify)', replicates_used='1 library per stage × fraction (no biological replicates)', input_control='none',
    processing_summary='minimap2 2.30 splice mode vs PlasmoDB v48 genome (PacBio -ax splice -uf -C5; ONT -ax splice -uf -k14); primary alignments to BED12; TAMA end-aware merge with the other three long-read sets; read-end QC (non-templated 5′ G, downstream A content, poly(A) tails); stage confirmation per model.',
    role='training label (long-read reference); stage confirmation', used_in='Fig 1A–L; Fig 2 labels; Table S2, S3',
    verification_notes='Only staged long-read set; primary reads per stage in TAMA collapse: ring 174,388 / troph 221,325 / schizont 270,156 / cytoplasmic 1,000; ONT 802,520 oriented primary reads (reference_v1_report.md). PMID and strain not checked.')
add(dataset_id='LR2', category='Long-read RNA', study='Lee VV et al.', year=2021, journal='mSystems 6:e01081-20', PMID='33688018', DOI='10.1128/mSystems.01081-20',
    accession='SRA SRR11094274 (P. falciparum direct RNA run cited by Shaw 2022; to verify that this is the only run used)', runs_or_samples_used=TV,
    assay='ONT direct RNA sequencing (poly(A)+ mRNA)', library_notes='3′-anchored (poly(A)); 5′ ends frequently truncated',
    stages_timepoints_used='mixed asexual stages', strain='3D7 (to verify)', replicates_used='1 (to verify)', input_control='none',
    processing_summary='Re-processed on HPCC ($HPCC_BIGDATA/longread_reprocess_2026/); collapsed and merged with TAMA; used for TES and junction support (130,267 primary reads).',
    role='training label (long-read reference; 3′ end and junction support)', used_in='Fig 1A–C; Table S2, S3', verification_notes='Lab directory LeeVV_2021_nanopore.')
add(dataset_id='LR3', category='Long-read RNA', study='Hoshizaki J et al.', year=2022, journal='BMC Genomics 23:780', PMID=TV, DOI='10.1186/s12864-022-09017-2',
    accession='ArrayExpress E-MTAB-11766; ENA ERR9746080', runs_or_samples_used='ERR9746080 (TEX+ and TEX− libraries are deposited merged)',
    assay='ONT direct RNA sequencing', library_notes='TEX+/TEX− cannot be separated in the deposited file, so TEX was not used as 5′-end evidence',
    stages_timepoints_used='mixed asexual stages', strain='3D7 (to verify)', replicates_used='1 merged run', input_control='none',
    processing_summary='As LR2 (87,998 primary reads).', role='training label (long-read reference; 3′ end and junction support)', used_in='Fig 1A–C; Table S2, S3',
    verification_notes='Lab directory Marcus_Nanopore (M.C.S. Lee group).')
add(dataset_id='LR4', category='Long-read RNA', study='Yang M et al.', year=2021, journal='Front Cell Infect Microbiol 11:631545', PMID='33708645', DOI='10.3389/fcimb.2021.631545',
    accession=TV, runs_or_samples_used=TV, assay='PacBio Iso-Seq (SMARTer cDNA)', library_notes='',
    stages_timepoints_used='mixed asexual stages (pool of 6 time points per the publication)', strain='3D7 (to verify)', replicates_used=TV, input_control='none',
    processing_summary='As LR2 (33,623 primary reads).', role='training label (long-read reference; 3′ end and junction support)', used_in='Fig 1A–C; Table S2, S3',
    verification_notes='Lab directory Mengquan_long_read. Accession not recorded in project files.')

# ---------------- short-read RNA inputs ----------------
add(dataset_id='RNA1', category='Short-read RNA', study='López-Barragán MJ et al.', year=2011, journal='BMC Genomics 12:587 (to verify)', PMID='22129310 (to verify)', DOI=TV,
    accession='SRA SRR364834–SRR364849 (lab directory RNAseq_7stages)', runs_or_samples_used='ring, trophozoite and schizont runs (exact run-to-stage map to verify)',
    assay='poly(A)+ RNA-seq', library_notes='Treated as unstranded for coverage; strand inferred from GT–AG splice motif; poly(A) sites from non-templated A soft-clips',
    stages_timepoints_used='ring, trophozoite, schizont (of 7 sequenced stages)', strain=TV, replicates_used=TV, input_control='none',
    processing_summary='Existing lab BAMs; spliced-read exon coverage, intron coverage and poly(A)-site channels per stage and inferred strand (18 channels: 3 types × 3 stages × 2 strands); raw coverage itself was dropped from the final (v3) inputs.',
    role='model input (18 inferred stranded RNA channels); stage confirmation (splice-junction support, to verify)', used_in='Fig 2C; Fig 3E (L_nostrand, N_noshortread); Table S6',
    verification_notes='Identity of this series as López-Barragán 2011 comes from project transcripts (SRR3648xx). Confirm that the stage-confirmation junction counts (junction_uniq_counts_by_stage.tsv) come from the same series.')
add(dataset_id='RNA2', category='Short-read RNA', study='Chappell L et al.', year=2020, journal='BMC Genomics 21:395', PMID='32513207', DOI='10.1186/s12864-020-06787-5',
    accession='ENA ERR2234550–ERR2234556 and ERR2234564–ERR2234570 (study ERP001570 per project notes; to verify)', runs_or_samples_used='14 runs: 0, 8, 16, 24, 32, 40, 48 h × 2 lanes',
    assay='DAFT-seq (amplification-free, directional RNA-seq)', library_notes='Reverse-stranded (read 1 antisense to RNA; 3–5% same-strand)',
    stages_timepoints_used='0–48 h every 8 h; model channels: ring = 8/16 h, trophozoite = 24/32 h, schizont = 40 h (per project notes)', strain='3D7',
    replicates_used='2 sequencing lanes per time point (technical)', input_control='none',
    processing_summary='HISAT2 alignment (HPCC flank_coverage_2026/bam); stranded coverage per stage; featureCounts (-p --countReadPairs -B -s 2) for gene expression; 0 h and 48 h align ~49%, others 75–86%.',
    role='model input (6 stranded channels); biology readout (expression, co-expression r, TSS-switch decomposition, expression strata)', used_in='Fig 2C; Fig 3E–F; Fig 4F; Fig 5E–H; Tables S6–S8',
    verification_notes='Runs from chappell_3D7_daftseq_runs.tsv.')
add(dataset_id='RNA3', category='Short-read RNA', study='Bunnik EM et al.', year=2013, journal='Genome Biol 14:R128', PMID='24267660', DOI=TV,
    accession='SRA SRP021890', runs_or_samples_used='Polysome-associated mRNA SRR836068 (0 h), SRR836069 (18 h), SRR836070 + SRR836071 (36 h)',
    assay='Polysome profiling (polysome-associated mRNA-seq)', library_notes='Unstranded; aligner not splice-aware; rRNA/structural-RNA reads removed before CPM',
    stages_timepoints_used='0 h (ring), 18 h (trophozoite), 36 h (schizont)', strain='3D7 (to verify)', replicates_used='ring 1, trophozoite 1, schizont 2', input_control='none',
    processing_summary='Re-downloaded and re-mapped paired-end (bowtie2), MAPQ ≥10, proper pairs, rRNA removed; 10-bp CPM bigWig. 0 h library aligns 22.4% (low quality).',
    role='model input (3 channels); biology readout (translation efficiency numerator)', used_in='Fig 2C; Fig 3E (L_nopoly); Fig S10; Table S7 (TE)', verification_notes='')
add(dataset_id='RNA4', category='Short-read RNA', study='Bunnik EM et al.', year=2013, journal='Genome Biol 14:R128', PMID='24267660', DOI=TV,
    accession='SRA SRP021890', runs_or_samples_used='Steady-state mRNA SRR835998 (0 h), SRR836066 (18 h), SRR836067 (36 h)', assay='mRNA-seq (steady-state, matched to polysome libraries)',
    library_notes='Different kit from polysome libraries (NuGEN Encore vs NEBNext)', stages_timepoints_used='0, 18, 36 h', strain='3D7 (to verify)', replicates_used='1 per time point',
    input_control='none', processing_summary='Downloaded from ENA (md5 checked), processed identically to RNA3; featureCounts on v68 CDS.',
    role='biology readout only (TE denominator)', used_in='Fig S10; Table S7 (TE)', verification_notes='')

# ---------------- nascent transcription / Pol II ----------------
add(dataset_id='NT1', category='Nascent transcription', study='Lu XM et al.', year=2017, journal='Nucleic Acids Res 45:7825–7840', PMID='28531310', DOI='10.1093/nar/gkx464',
    accession='GEO GSE85478', runs_or_samples_used='Author-normalised coverage files for ER, LR, ET, LT, ES, LS', assay='GRO-seq', library_notes='Deposited files carry no strand; treated as unstranded',
    stages_timepoints_used='ring = mean(early, late ring); trophozoite = mean(ET, LT); schizont = mean(ES, LS)', strain='3D7 (to verify)', replicates_used='2 time points per stage (no replicate within time point)',
    input_control='none', processing_summary='Author-supplied per-base coverage (scale factor already applied) averaged per stage; chromosome names harmonised. Early ring contributes only 9.5% of the ring total.',
    role='model input (3 channels)', used_in='Fig 2C; Fig 3E (L_nogro)', verification_notes='')
add(dataset_id='NT2', category='Nascent transcription', study='Lu XM et al.', year=2017, journal='Nucleic Acids Res 45:7825–7840', PMID='28531310', DOI='10.1093/nar/gkx464',
    accession='GEO GSE85478', runs_or_samples_used='PolII_ER, PolII_ET, PolII_LS; input_PolII_LS (GSM IDs to verify)', assay='ChIP-seq, RNA Pol II CTD phospho-Ser2 (Abcam ab5095)',
    library_notes='Elongating Pol II; gene-body enriched', stages_timepoints_used='early ring, early trophozoite, late schizont', strain='3D7 (to verify)', replicates_used='1 IP per stage',
    input_control='single late-schizont input used for all three stages', processing_summary='deeptools bamCoverage, 10 bp, CPM; input subtracted; negatives clipped to 0.',
    role='model input (3 channels; one of the 8 chromatin tracks)', used_in='Fig 2C; Fig 3E; Table S6 (TL_PolIIS2P)', verification_notes='')

# ---------------- chromatin ----------------
add(dataset_id='CH1', category='Chromatin', study='Tang J et al. (Duffy lab)', year=2020, journal=TV, PMID='33225957', DOI=TV,
    accession='BioProject PRJNA612099; SRA SRP252482', runs_or_samples_used='H2A.Z: R SRR11292863, SRR11292895; T SRR11292878, SRR11292885; S SRR11292862, SRR11292870. H3K27ac: R SRR11292854, SRR11292874; T SRR11292879, SRR11292886; S SRR11292864, SRR11292871. Inputs: R SRR11292857, SRR11292897; T SRR11292891, SRR11292892; S SRR11292867, SRR11292875',
    assay='ChIP-seq H2A.Z, H3K27ac (H3, H3K4me1 used only in ablations)', library_notes='', stages_timepoints_used='ring, trophozoite, schizont', strain=TV,
    replicates_used='2 IP / 2 input per stage', input_control='matched input per stage', processing_summary='Lab BAMs; deeptools bamCoverage 10 bp CPM; replicates averaged, then input subtracted, negatives clipped to 0.',
    role='model input (H2A.Z, H3K27ac); ablation only (H3, H3K4me1)', used_in='Fig 2C; Fig 3E; Table S6', verification_notes='Lab directory Jingyi_Chip.')
add(dataset_id='CH2', category='Chromatin', study='Hoeijmakers WAM et al.', year=2013, journal=TV, PMID='23320541', DOI=TV, accession='GEO GSE39702; SRA SRP014595',
    runs_or_samples_used='R SRR527826 / input SRR527829; T SRR527827 / SRR527830; S SRR527828 / SRR527831', assay='ChIP-seq H2B.Z', library_notes='', stages_timepoints_used='ring, trophozoite, schizont',
    strain=TV, replicates_used='1 IP / 1 input per stage', input_control='matched input', processing_summary='As CH1.', role='model input', used_in='Fig 2C; Fig 3E; Table S6 (TL_H2BZ)', verification_notes='')
add(dataset_id='CH3', category='Chromatin', study='Bártfai R et al.', year=2010, journal='PLoS Pathog (to verify)', PMID='21187892', DOI=TV, accession='GEO GSE23867; SRA SRP003508',
    runs_or_samples_used='H3K4me3: R SRR065656; T SRR065657, SRR065658; S SRR065659. H3K9ac: R SRR065652; T SRR065653, SRR065654; S SRR065655. Inputs: R SRR065660; T SRR065661, SRR065662; S SRR065664',
    assay='ChIP-seq H3K4me3, H3K9ac', library_notes='', stages_timepoints_used='ring, trophozoite, schizont', strain=TV, replicates_used='ring 1, trophozoite 2, schizont 1 (IP and input)',
    input_control='matched input', processing_summary='As CH1.', role='model input', used_in='Fig 2C; Fig 3E; Table S6 (TL_H3K4me3, TL_H3K9ac)', verification_notes='Lab directory Richard_2010_Chip.')
add(dataset_id='CH4', category='Chromatin', study='Toenhake CG et al.', year=2018, journal='Cell Host Microbe 23:557 (to verify)', PMID='29649445', DOI=TV, accession='GEO GSE104075; SRA SRP118379',
    runs_or_samples_used='R SRR6055327–29, SRR6689444–46; T SRR6055330–31, SRR6689447–48; S SRR6055332–34, SRR6689449–51; gDNA controls SRR6055335, SRR6055336',
    assay='ATAC-seq', library_notes='', stages_timepoints_used='ring, trophozoite, schizont (time points pooled per stage; time-point map to verify)', strain=TV,
    replicates_used='ring 6, trophozoite 4, schizont 6 libraries', input_control='2 genomic-DNA (naked DNA) libraries subtracted', processing_summary='As CH1 (10 bp CPM, control subtracted, clipped at 0).',
    role='model input; biology readout (Fig 5A, divergent-pair NDR calls)', used_in='Fig 2C; Fig 3E; Fig 5A; Fig S7–S8; Table S6, S8', verification_notes='')
add(dataset_id='CH5', category='Chromatin', study='Bunnik EM et al.', year=2014, journal='BMC Genomics (to verify)', PMID='24885191', DOI=TV, accession='SRA SRP026365',
    runs_or_samples_used='MNase: SRR923440 (0 h), SRR923449 (18 h), SRR923450 (36 h); sonication reference SRR923437–SRR923439', assay='MNase-seq (paired-end 2×51 bp)',
    library_notes='Read coverage without fragment extension (dyad resolution is low)', stages_timepoints_used='0 h ring, 18 h trophozoite, 36 h schizont', strain='3D7 (to verify)',
    replicates_used='1 per stage', input_control='none subtracted (sonication libraries kept as reference only)', processing_summary='bamCoverage 10 bp CPM, no extension, no control subtraction, no clipping (canonical IP-only tracks; manifest_mnase_canonical.csv).',
    role='model input', used_in='Fig 2C; Fig 3E; Table S6 (TL_MNase)', verification_notes='Not the same data as Kensche 2016 (BR2).')
add(dataset_id='CH6', category='Chromatin (ablation only)', study='H3K9me3: three studies (lab directories Myriam_2022, Riward_2022 = Morillo RC et al., Subudhi_2023)', year='2022–2023', journal=TV,
    PMID='37013531 (ring); 36379668 (schizont); bioRxiv preprint (trophozoite)', DOI=TV, accession='GEO GSE202214 (ring); GSE210062 (schizont); GSE230206 (trophozoite, 20 hpi DMSO arm)',
    runs_or_samples_used='R SRR19082629, -36, -37, SRR19082703; S SRR20710201, -04, -07; T SRR24234015, -16 (inputs in FILE_CATALOG.csv)', assay='ChIP-seq / CUT&RUN H3K9me3',
    library_notes='Cross-study; amplitudes not comparable across stages', stages_timepoints_used='ring, trophozoite, schizont', strain='differs by study (to verify)', replicates_used='R 4/6, S 3/3, T 2/2 (IP/control)',
    input_control='matched per study', processing_summary='As CH1.', role='ablation only (TA_H3K9me3; not in the final model)', used_in='Table S6', verification_notes='')
add(dataset_id='CH7', category='Chromatin / TF', study='Shang X et al.', year=2022, journal='Nucleic Acids Res (to verify)', PMID='35288749', DOI=TV, accession='GEO GSE184658, GSE184659; SRA SRP339069',
    runs_or_samples_used='13 ApiAP2 factors × ring/trophozoite/schizont = 39 tracks (AP2-06B, -08, -11A, -11B, -12, -G2, -G3, -G5, -HC, -HS, -L, -O5, -exp); HP1; GFP mock IP (schizont only)',
    assay='ChIP-seq of GFP-tagged ApiAP2 factors and HP1', library_notes='', stages_timepoints_used='ring, trophozoite, schizont', strain='3D7 transgenic lines (to verify)',
    replicates_used='2 IP / 1 input per track', input_control='matched input per IP pair (meta/ip_control_pairs.tsv)', processing_summary='Lab BAMs (MAPQ ≥30, deduplicated); bamCoverage 10 bp CPM, IP replicates averaged, input subtracted, clipped at 0. Track classes (29 gene-body, 10 promoter-type) and GC-bias QC in Table S9.',
    role='biology readout only (10 promoter-type tracks, Fig 5B); ablation only (TA_AP2*, TA_HP1)', used_in='Fig 5B; Fig S6–S7; Table S6, S9', verification_notes='AP2-11B gene ID is inferred (see Table S9).')

# ---------------- independent 5' / 3' end validation ----------------
add(dataset_id='V1', category='5′-end validation', study='Adjalley SH et al.', year=2016, journal='Cell Rep 14:2463 (to verify)', PMID='26947071', DOI='10.1016/j.celrep.2016.02.025',
    accession='GEO GSE68982 (GSM1689659–GSM1689670)', runs_or_samples_used='repA/repB × 6 time points', assay='CAGE-like 5′-cap sequencing (cap-seq)', library_notes='Stranded, single-nucleotide 5′ ends',
    stages_timepoints_used='2, 10, 18, 26, 34, 42 hpi (2/10 → ring, 18/26 → trophozoite, 34/42 → schizont)', strain='3D7', replicates_used='2 biological replicates per time point', input_control='none',
    processing_summary='CAGEr; replicates/time points pooled per stage; CTSS bigWigs CPM-normalised per stage (one denominator for both strands); distclu clusters at 20 bp. Per-gene dominant TSS for benchmarking (4,096 genes).',
    role='validation only (never used for labels or inputs in the final model); biology readout (stage TSS usage)', used_in='Fig 1D; Fig 3C; Fig 4A, C, F; Fig 5H–I; Tables S5, S7',
    verification_notes='Earlier (pre-v3) reference tiers used CAGE; the final v3 reference and labels are CAGE-free (see manuscript).')
add(dataset_id='V2', category='5′-end validation', study='Chappell L et al.', year=2020, journal='BMC Genomics 21:395', PMID='32513207', DOI='10.1186/s12864-020-06787-5',
    accession='Supplementary tag-cluster table (file chappell2020_5UTRseq_tagclusters.tsv); raw runs "Pf3D7_FATEseq_v1" (14 runs; accessions to verify)', runs_or_samples_used='7 time points × 2 lanes',
    assay='5′UTR-seq (amplification-free, Smart-seq TSO-anchored)', library_notes='TSS tags = reads containing the TSO (GTATCAACGCAGAGTAC)', stages_timepoints_used='0–48 h every 8 h (TP1–TP7 assumed = 0…48 h)', strain='3D7',
    replicates_used='2 lanes per time point (technical)', input_control='none', processing_summary='Benchmark: published tag clusters (Quantile_TPM), per-gene dominant cluster (3,029 genes). Decomposition: raw FASTQ re-processed, TSO-containing reads only (37k–95k TSS tags per lane).',
    role='validation only; biology readout (stage TSS usage)', used_in='Fig 4A, C; Fig 5H–I; Tables S5, S7', verification_notes='PlasmoDB v71 UTRs derive from Chappell 2020 coverage, so this set is not fully independent of v71.')
add(dataset_id='V3', category='5′-end validation', study='Shaw PJ et al.', year=2021, journal='PeerJ 9:e11983', PMID=TV, DOI='10.7717/peerj.11983', accession='GEO GSE103036',
    runs_or_samples_used='5′-capped nucleotide (5CN) calls (file shaw2021_5CAPture_5CN_clusters.tsv)', assay='5CAPture-seq (eIF4E-enriched full-length cDNA vs unenriched control)', library_notes='',
    stages_timepoints_used='ring, trophozoite, schizont', strain='3D7 (to verify)', replicates_used='replicate experiments combined by the authors', input_control='unenriched cDNA control',
    processing_summary='5CN sites with p < 0.01, merged within 20 bp into clusters; per-gene dominant cluster (1,806 genes).', role='validation only', used_in='Fig 4A, C; Table S5', verification_notes='')
add(dataset_id='V4', category='5′-end validation', study='Consensus TSS (this study)', year=2026, journal='—', PMID='—', DOI='—', accession='derived from V1–V3',
    runs_or_samples_used='—', assay='derived', stages_timepoints_used='—', strain='—', replicates_used='—', input_control='—',
    processing_summary='Gene-level TSS supported by ≥2 of V1–V3 within 100 bp (1,649 genes).', role='validation only (primary TSS truth set)', used_in='Fig 4A, F; Table S5', verification_notes='')
add(dataset_id='V5', category='3′-end validation', study='Siegel TN et al.', year=2014, journal='BMC Genomics 15:150', PMID='24559473', DOI='10.1186/1471-2164-15-150', accession=TV,
    runs_or_samples_used='Published poly(A)-site table (supplementary; project file 12864_2013_7013_MOESM10_ESM.xlsx, to verify)', assay='Strand-specific RNA-seq (RNA-ligation based); poly(A)-tailed reads define PAS',
    library_notes='', stages_timepoints_used='blood-stage time points pooled (not stage-resolved)', strain='3D7 (to verify)', replicates_used=TV, input_control='none',
    processing_summary='Per-gene PAS clusters; major PAS used as TTS truth (3,439 genes with ≥1 cluster in our gene set).', role='validation only (TTS truth set)', used_in='Fig 4B, C; Table S5', verification_notes='')
add(dataset_id='V6', category='3′-end validation (supplementary)', study=TV, year=TV, journal=TV, PMID=TV, DOI=TV, accession='BioProject PRJNA1176112', runs_or_samples_used=TV, assay='PAT-seq (poly(A)-tail anchored 3′-end sequencing)',
    library_notes='93% of anchored sites have A-rich downstream sequence (internal priming)', stages_timepoints_used=TV, strain=TV, replicates_used=TV, input_control='none',
    processing_summary='Anchored cleavage sites; signal at predicted TTS (patseq_signal_at_tts_summary.tsv).', role='validation only (supplementary; not used for main claims)', used_in='Table S5 (PATseq sheet)',
    verification_notes='Study identity not recorded in project files.')

# ---------------- biology readouts ----------------
add(dataset_id='BR1', category='Translation', study='Caro F et al.', year=2014, journal='eLife 3:e04106', PMID='25493618', DOI='10.7554/eLife.04106', accession='GEO GSE58402 (GSM1410296–GSM1410299)',
    runs_or_samples_used='SRR1378560–63 and SRR2014726–29 (two sequencing runs of the same 4 samples, merged)', assay='Ribosome profiling (Ribo-seq)', library_notes='No 3-nt periodicity in these data',
    stages_timepoints_used='ring 11 h; early/late trophozoite 21/31 h (merged); schizont 45 h', strain='W2', replicates_used='1 per time point', input_control='—',
    processing_summary='P-site counts in distal-unique leaders, proximal UTRs, uORFs and CDS.', role='biology readout only', used_in='Fig S10; Table S7 (decomposition ribo_* columns)', verification_notes='')
add(dataset_id='BR2', category='Nucleosomes', study='Kensche PR et al.; re-processed by Holzinger et al. (nucDetective)', year='2016; 2026', journal='Nucleic Acids Res 44:2110 (to verify); PLoS Comput Biol (2026)',
    PMID='26578577 (Kensche)', DOI='10.1093/nar/gkv1214; 10.1371/journal.pcbi.1014557; Zenodo 10.5281/zenodo.16779899 (v1.1)', accession='GEO GSE66185',
    runs_or_samples_used='nucDetective output_nucDetective.zip and Pf_nucleosome_annotation (md5 checked); PLOS S1 Table (+1 nucleosome annotation)', assay='Paired-end MNase-seq',
    library_notes='', stages_timepoints_used='T5, T10, T15, T20, T25, T30, T35, T40 (8 IDC time points; averaged)', strain='3D7', replicates_used=TV, input_control='—',
    processing_summary='qNorm coverage, DANPOS dyads and reference nucleosome map from nucDetective; PlasmoDB-57 coordinates (sequence identical to v68). NDR = smoothed occupancy < 0.5 × genome median, ≥50 bp.',
    role='biology readout only (independent of model inputs)', used_in='Fig 5C, F, G; Fig S8; Table S8', verification_notes='')
add(dataset_id='BR3', category='Chromatin states', study='Brown AS Jr, Llinás M, Mahony S', year=2026, journal='BMC Genomics 27:142', PMID='41501628', DOI='10.1186/s12864-025-12455-3',
    accession='GitHub seqcode/pfal-chrom-states-trackhub (chrom-states-IGV/{10,30,40}hpi_chrom-states.bed)', runs_or_samples_used='11-state segmentations at 10, 30, 40 hpi', assay='HMM chromatin-state segmentation (200-bp bins)',
    library_notes='', stages_timepoints_used='10, 30, 40 hpi', strain=TV, replicates_used='—', input_control='—', processing_summary='Used as published (no liftover needed); TSS-associated states 2/5/6/8 scored as a TSS predictor.',
    role='benchmark comparison / biology readout only', used_in='Table S5 (HMM_states sheet); supplementary chromatin-state analysis', verification_notes='bioRxiv preprint PMID 40909583.')

# ---------------- reference resources ----------------
add(dataset_id='REF1', category='Reference', study='PlasmoDB / VEuPathDB P. falciparum 3D7', year='release 48 / 68 / 71', journal='—', PMID='—', DOI='—', accession='PlasmoDB-48 genome FASTA; PlasmoDB-68 GFF (gene structures identical to v71)',
    runs_or_samples_used='5,285 nuclear protein-coding genes', assay='genome and annotation', stages_timepoints_used='—', strain='3D7', replicates_used='—', input_control='—',
    processing_summary='v48 and v68 nuclear sequences are identical; v71 annotation used as the comparator (v71 = v68 gene coordinates).', role='reference', used_in='all', verification_notes='')
add(dataset_id='REF2', category='Model resource', study='Schiff Y et al. (Caduceus)', year=2024, journal='ICML 2024 (to verify)', PMID='—', DOI=TV, accession=TV,
    runs_or_samples_used='pretrained checkpoint (name to verify)', assay='DNA language model embeddings', stages_timepoints_used='—', strain='—', replicates_used='—', input_control='—',
    processing_summary='Embeddings reduced to 16 PCs, forward + reverse complement (32 channels).', role='model input (gLM channels)', used_in='Fig 2C; Fig 3E (L_noglm)', verification_notes='')

# ---------------- considered but not used ----------------
add(dataset_id='NU1', category='Not used', study='Reers AB et al.', year=2023, journal=TV, PMID='37322481', DOI=TV, accession='GEO GSE215426 / GSE215429', assay='ChIP-seq (merozoite/schizont/ring) and chromatin states',
    processing_summary='Segmentation not public; not re-run.', role='not used', used_in='—', verification_notes='Listed for completeness.')
add(dataset_id='NU2', category='Not used', study='Hi-C (lab directory LLM-study/HiC_data_preprocess2)', year=TV, journal=TV, PMID='24671853', DOI=TV, accession='GEO GSE50199; SRA SRP029967', assay='Hi-C',
    processing_summary='Processed tracks exist in the lab track folder but are not model inputs and are not analysed in this manuscript.', role='not used', used_in='—', verification_notes='')

S1 = pd.DataFrame(R, columns=COLS)
desc = {'dataset_id': 'Identifier used in this table (LR = long-read, RNA = short-read RNA, NT = nascent transcription, CH = chromatin, V = validation, BR = biology readout, REF = reference, NU = not used).',
        'category': 'Data type.', 'study': 'First author / source.', 'year': 'Publication year.', 'journal': 'Journal, volume and article number where known.', 'PMID': 'PubMed ID ("to verify" = not confirmed from project records).',
        'DOI': 'DOI of the publication or data deposit.', 'accession': 'Repository accession(s) (GEO/SRA/ENA/ArrayExpress/Zenodo/GitHub).', 'runs_or_samples_used': 'Runs or samples actually used in this study.',
        'assay': 'Assay.', 'library_notes': 'Library properties relevant to the analysis (strandedness, biases).', 'stages_timepoints_used': 'Stages / time points used and how they were mapped to ring, trophozoite, schizont.',
        'strain': 'Parasite strain.', 'replicates_used': 'Replicates used (IP/control counts for ChIP).', 'input_control': 'Input or control used for normalisation.',
        'processing_summary': 'Processing applied in this study.', 'role': 'Role in this study: model input / training label / validation only / biology readout only / reference / not used.',
        'used_in': 'Figures and supplementary tables that use the dataset.', 'verification_notes': 'Caveats and fields still to verify.'}
intro = ['Every external dataset used in the study, its role, and how it was processed.',
         'Roles: "model input" = channel of the final boundary model (V3_top8: 90 channels; see Model_inputs sheet); "training label" = used to build the long-read reference and labels; "validation only" = never used for labels or inputs; "biology readout only" = used only in the TSS-anchored biology analyses.',
         'Fields marked "to verify" could not be confirmed from project records (DATA_DOCUMENTATION.md, FILE_CATALOG.csv, manifest_mnase_canonical.csv, chappell_3D7_daftseq_runs.tsv, analysis reports) and should be checked against the publications before submission.',
         'Sources: $LAB_DIR/lncRNA project tracks/DATA_DOCUMENTATION.md, FILE_CATALOG.csv, assay_summary.csv, accessibility_nucleosome/manifest_mnase_canonical.csv, cage_tss/manifest_cage_cpm.csv; project reports final_model_labels_inputs.md, UTR_BOUNDARY_MODEL_DESIGN.md, reference_v1_report.md, README_data_layout.md, report_tss_expression_decomposition.md, report_uorf_translation.md, report_divergent_mnase_ndr.md, report_chromatin_states_nucleosomes.md, report_apiap2_genebody.md; manuscript draft v1 (Manuscript_draft_CN_v1.md). Accessions/DOIs for Lee 2021, Yang 2021, Shaw 2021, Caro 2014 and Siegel 2014 were checked against publisher/PubMed pages (Oct 2026).']
MI = pd.DataFrame([
    ('Chromatin (unstranded)', 'H2A.Z, H2B.Z, H3K9ac, H3K27ac, H3K4me3, ATAC-seq, MNase-seq, Pol II Ser2P × ring/trophozoite/schizont', 24, 'CH1–CH5, NT2'),
    ('Nascent / translated RNA (unstranded)', 'GRO-seq × 3 stages; polysome profiling × 3 stages', 6, 'NT1, RNA3'),
    ('Stranded RNA, measured', 'DAFT-seq × 3 stages × 2 strands', 6, 'RNA2'),
    ('Stranded RNA, inferred from short-read mRNA-seq', 'spliced-read exon coverage, intron coverage (strand from GT–AG motif), poly(A) sites (non-templated A soft-clips) × 3 stages × 2 strands', 18, 'RNA1'),
    ('DNA', 'one-hot base fraction per 20-bp bin', 4, 'REF1'),
    ('DNA language model', 'Caduceus embeddings, 16 PCs, forward + reverse complement', 32, 'REF2'),
    ('Not inputs', 'CAGE (validation only); raw mRNA-seq coverage (dropped in v3); H3K4me1, H3, H3K9me3, HP1 and 13 ApiAP2 factors (tested in ablations only)', 0, 'V1, CH1, CH6, CH7')],
    columns=['group', 'tracks', 'n_channels', 'dataset_ids'])
mid = {'group': 'Input group of the final model (V3_top8).', 'tracks': 'Tracks in the group.', 'n_channels': 'Number of input channels (total 90).', 'dataset_ids': 'dataset_id(s) in the Data_sources sheet.'}
path, miss = write_table('TableS1_data_sources.xlsx', 'Table S1. Data sources', intro,
    [dict(name='Data_sources', df=S1, desc=desc, source='compiled from the files listed above', note=f'{len(S1)} datasets/resources'),
     dict(name='Model_inputs', df=MI, desc=mid, source='final_model_labels_inputs.md', note='90 channels in total')], coord=False)
print(path, miss, (S1 == TV).sum().sum() + S1.apply(lambda c: c.astype(str).str.contains('to verify')).sum().sum())
