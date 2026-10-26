#!/bin/bash
#SBATCH --job-name=energylevelsalpha
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --account=MATH022462
#SBATCH --output=logs/job_%j.out
#SBATCH --error=logs/job_%j.err

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

# Run with 14 parallel workers (4.5GB RAM per process)
mpirun -np 14 python -u -m mpi4py.futures energylevelsalpha.py

conda deactivate