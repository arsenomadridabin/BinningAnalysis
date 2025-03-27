#!/bin/bash
#SBATCH -N 1                    # number of nodes
#SBATCH -n 20
##SBATCH -c 12                  # 6 threads per MPI process
#SBATCH -t 10:00:00
#SBATCH -p workq
#SBATCH -A hpc_bb_karki2
#SBATCH -o  gb.out
#SBATCH -e  error.out


PYTHONPATH=/scratch/ashaky3/dpMD/cpu-2.0.3/bin/python
source /scratch/ashaky3/.env/bin/activate

python standard_femo_5.py -n 33280 -s 68 -b 512 -i fe.json -o fe -k merge_analysis -p 1 -e 0

