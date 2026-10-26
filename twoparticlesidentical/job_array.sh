#!/bin/bash
#SBATCH --job-name=energylevels
#SBATCH --partition=compute
#SBATCH --array=1-11
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --account=MATH022462
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

export PYTHONUNBUFFERED=1

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

export I_MPI_HYDRA_BOOTSTRAP=fork
export FI_PROVIDER=tcp

. ~/initConda.sh
source activate /user/work/vq22287/anaconda3/envs/py312

# Pass array ID to energylevelsalpha.py
mpirun -np 14 python -u -m mpi4py.futures "$SLURM_SUBMIT_DIR/energylevelsalpha.py" $SLURM_ARRAY_TASK_ID

conda deactivate