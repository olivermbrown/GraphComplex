#!/bin/bash
#SBATCH --job-name=test_run
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --mem=16G
#SBATCH --account=MATH022462
#SBATCH --output=logs/test_%j.out
#SBATCH --error=logs/test_%j.err

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Use TCP fabric provider verified on BluePebble
export FI_PROVIDER=tcp

. ~/initConda.sh
source activate /user/work/vq22287/anaconda3/envs/py312

echo "=== SLURM BATCH DEBUG HEADER ==="
echo "Running on host: $(hostname)"
echo "Current directory: $(pwd)"
echo "Python binary path: $(which python)"
echo "================================="

# Launch 28 MPI workers on single node
mpirun -np 28 python -u -m mpi4py.futures test_energylevels.py

conda deactivate