source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
# common environment for ablation jobs
cd ${PF_HPCC_BIGDATA}/ablation_2026/work
export PY=${PF_HPCC_BIGDATA}/.conda/envs/enformer_pure/bin/python
export PYTHONPATH=.:boundary_model POOL=1 ML_LABELS=data/Lml3soft40_pool.npz PE_LABELS=data/Lml3_pool.npz PE_POS=data/ml3_end_positives.tsv
export TMPDIR=${PF_HPCC_BIGDATA}/ablation_2026/tmp OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
mkdir -p $TMPDIR logs
$PY -c "import torch, pandas, sklearn; print('env ok', torch.__version__)" || exit 1
