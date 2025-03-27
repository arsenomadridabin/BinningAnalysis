#!/bin/bash
#SBATCH -N 1                    # number of nodes
#SBATCH -t 1:00:00
#SBATCH -n 20
#SBATCH -p single
#SBATCH -A hpc_bb_karki2
#SBATCH -o  debug.out
#SBATCH -e  err.out


source /scratch/ashaky3/.env/bin/activate

python to_json.py

python standard_femo_5.py -n 32768 -s 68 -b 512 -i fe.json -o fe -k average_distribution -p 1 -e 0 
python standard_femo_5.py -n 32768 -s 68 -b 512 -i mg.json -o mg -k average_distribution -p 1 -e 0 
python standard_femo_5.py -n 32768 -s 68 -b 512 -i si.json -o si -k average_distribution -p 1 -e 0 
python standard_femo_5.py -n 32768 -s 68 -b 512 -i o.json -o o -k average_distribution -p 1 -e 0 
python standard_femo_5.py -n 32768 -s 68 -b 512 -i h.json -o h -k average_distribution -p 1 -e 0

python standard_femo_5.py -n 32768 -s 68 -b 512 -i fe.json -j mg.json -o mg_in_fe -p 1 -k get_sub_atom_count_in_fe -e 0 -x 5
python standard_femo_5.py -n 32768 -s 68 -b 512 -i fe.json -j si.json -o si_in_fe -p 1 -k get_sub_atom_count_in_fe -e 0 -x 5
python standard_femo_5.py -n 32768 -s 68 -b 512 -i fe.json -j o.json -o o_in_fe -p 1 -k get_sub_atom_count_in_fe -e 0 -x 5
python standard_femo_5.py -n 32768 -s 68 -b 512 -i fe.json -j fe.json -o fe_in_fe -p 1 -k get_sub_atom_count_in_fe -e 0 -x 5
python standard_femo_5.py -n 32768 -s 68 -b 512 -i fe.json -j h.json -o h_in_fe -p 1 -k get_sub_atom_count_in_fe -e 0 -x 5

python script_generate_composition.py --low_cut_off_metal 38 --high_cut_off_silicate 8

wait
