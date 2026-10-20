#!/bin/bash
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
source /opt/linux/rhel/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/etc/profile.d/conda.sh
E=${PF_HPCC_BIGDATA}/.conda/envs/tama_py27
conda create -y -p $E -c conda-forge -c bioconda --override-channels python=2.7 biopython pysam > ${PF_HPCC_BIGDATA}/GSE186109_longread/logs/tama_env.log 2>&1
echo "RC=$?" >> ${PF_HPCC_BIGDATA}/GSE186109_longread/logs/tama_env.log
$E/bin/python -c "import Bio, pysam; print('OK', Bio.__version__, pysam.__version__)" >> ${PF_HPCC_BIGDATA}/GSE186109_longread/logs/tama_env.log 2>&1
touch ${PF_HPCC_BIGDATA}/GSE186109_longread/logs/TAMA_ENV_DONE
