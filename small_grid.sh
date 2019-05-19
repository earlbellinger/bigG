#!/bin/bash
#SBATCH --job-name=small_grid
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
for ii in `seq -0.05 0.01 0.05`; do
    counter=$(echo "$counter + 1" | bc -l)
    if [[ ! $SLURM_ARRAY_TASK_ID -eq $counter ]]; then continue; fi 
    if [[ "$ii" = "-0.00" ]]; then ii=0.00; fi
    ./astec_track.sh -d small_grid -n $ii -M 0.8 -Y 0.25 -Z 0.005 -a 1.8 -t 11.e9 -b $ii
done 

echo "========= Job finished at `date` =========="
#
