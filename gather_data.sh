#!/bin/bash
#SBATCH -N 1                    # number of nodes
#SBATCH -n 18
#SBATCH -t 6:00:00
#SBATCH -p gpu
#SBATCH -A hpc_bb_karki2
#SBATCH -o  output.out
#SBATCH -e  err.out

PYTHONPATH=/project/ashaky3/dpmd3/bin/python


source /scratch/ashaky3/.env/bin/activate

python standard_femo_5.py -n 33280 -s 68 -b 512 -l out.dump -o o -k parse_dpmd_data -p 1 -e 0
