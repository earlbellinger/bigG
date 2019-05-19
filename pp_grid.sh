#!/bin/bash
#SBATCH --job-name=pp_grid
#SBATCH --partition=q8
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --array=1-11%11
#SBATCH --output=./slurm_logs/%A_%a.out

echo "========= Job started  at `date` =========="

echo "My jobid: $SLURM_JOB_ID"
echo "My array id: $SLURM_ARRAY_TASK_ID"

module load anaconda3

counter=0
for ii in `seq 0.95 0.01 1.05`; do
    counter=$(echo "$counter + 1" | bc -l)
    if [[ ! $SLURM_ARRAY_TASK_ID -eq $counter ]]; then continue; fi 
    #if [[ "$ii" = "-0.00" ]]; then ii=0.00; fi
    ./astec_track.sh -d pp_grid -n $ii -M 0.75 -Y 0.25 -Z 0.001 -a 1.8 -t 11.e9 -r1 $ii
done 

echo "========= Job finished at `date` =========="
#
