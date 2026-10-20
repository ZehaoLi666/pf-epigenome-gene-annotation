"""Session-derived computation of hc_reference.bed.gz.

Recovered from the record of the interactive analysis session in which hc_reference.bed.gz was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: (none recorded)
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
# hc_reference.bed.gz is produced on a remote HPC cluster (ssh:hpcc) and cannot be
# reproduced in a standalone Python session. The file is created by:
#
# 1. Running TAMA collapse + merge jobs on HPCC (SLURM array jobs tama_collapse_up.sbatch,
#    tama_collapse_relaxed.sbatch, ont_orient_align.sbatch, merge_evidence.sbatch)
# 2. Running evidence_tiers.py on the cluster which writes hc_reference.bed, then
#    gzip -kf hc_reference.bed
#
# The resulting file was transferred to the local workspace via base64 chunked download
# from $HPCC_BIGDATA/GSE186109_longread/reference_v1/pack.tgz
# and extracted with:

import tarfile, os

os.makedirs("ref_v1", exist_ok=True)
tarfile.open("ref_v1/pack.tgz").extractall("ref_v1")

import shutil
shutil.copy("ref_v1/hc_reference.bed.gz", "hc_reference.bed.gz")
