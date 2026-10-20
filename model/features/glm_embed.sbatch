#!/bin/bash
#SBATCH --account=lerochlab
#SBATCH --partition=short_gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --job-name=pf_glm_embed
#SBATCH --output=logs/glm_embed_%j.log
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
export HF_HOME=${PF_HPCC_HOME}/.cache/huggingface HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1
export TMPDIR=${PF_HPCC_BIGDATA}/boundary_seg/tmp; mkdir -p $TMPDIR
${PF_HPCC_BIGDATA}/.conda/envs/caduceus_min_env/bin/python ${PF_HPCC_BIGDATA}/boundary_seg/scripts/glm_embed.py
echo EXIT $?
