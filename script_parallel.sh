#!/bin/bash
#SBATCH -N 1                    # Number of nodes
#SBATCH -n 20                   # Number of tasks (MPI processes)
#SBATCH -t 10:00:00             # Max run time
#SBATCH -p workq
#SBATCH -A hpc_bb_karki2
#SBATCH -o gb1.out
#SBATCH -e error1.out

source /scratch/ashaky3/.env/bin/activate

# Do not run in this in parallel with beow jobs ...dependency issue
#srun python process_atom.py $SLURM_ARRAY_TASK_ID

python to_json_parallel.py

# Run below jobs in parallel

python standard_femo_5.py -n 32768 -s 68 -b 512 -i fe.json -o fe -k average_distribution -p 1 -e 0 &
python standard_femo_5.py -n 32768 -s 68 -b 512 -i mg.json -o mg -k average_distribution -p 1 -e 0 &
python standard_femo_5.py -n 32768 -s 68 -b 512 -i si.json -o si -k average_distribution -p 1 -e 0 &
python standard_femo_5.py -n 32768 -s 68 -b 512 -i o.json -o o -k average_distribution -p 1 -e 0 &
#python standard_femo_5.py -n 32768 -s 68 -b 512 -i h.json -o h -k average_distribution -p 1 -e 0 &

python standard_femo_5.py -n 32768 -s 68 -b 512 -i fe.json -j mg.json -o mg_in_fe -p 1 -k get_sub_atom_count_in_fe -e 0 -x 5 &
python standard_femo_5.py -n 32768 -s 68 -b 512 -i fe.json -j si.json -o si_in_fe -p 1 -k get_sub_atom_count_in_fe -e 0 -x 5 &
python standard_femo_5.py -n 32768 -s 68 -b 512 -i fe.json -j o.json -o o_in_fe -p 1 -k get_sub_atom_count_in_fe -e 0 -x 5 &
python standard_femo_5.py -n 32768 -s 68 -b 512 -i fe.json -j fe.json -o fe_in_fe -p 1 -k get_sub_atom_count_in_fe -e 0 -x 5 &
#srun -n 1 python standard_femo_5.py -n 32768 -s 68 -b 512 -i fe.json -j h.json -o h_in_fe -p 1 -k get_sub_atom_count_in_fe -e 0 -x 5 &

wait  # Ensure all parallel jobs finish before running this final script

# Run the final script
python script_generate_composition.py --config config.json
