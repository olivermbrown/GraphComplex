#!/bin/bash
#SBATCH --job-name=ygraph_test
#SBATCH --partition=veryshort
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:0:0
#SBATCH --mem=100M

echo python save_funcs.py