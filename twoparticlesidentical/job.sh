#!/bin/bash
#SBATCH --job-name=Hgraph_hardcore_all_alphas_varying
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:0:0
#SBATCH --mem=20G
#SBATCH --account=MATH022462
#SBATCH --array=0-14

. ~/initConda.sh                    

source activate /user/work/vq22287/anaconda3/envs/py312

python savefuncs.py ${SLURM_ARRAY_TASK_ID}

conda deactivate