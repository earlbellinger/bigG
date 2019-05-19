#!/bin/bash
#SBATCH --job-name=sobol_dispatcher
#SBATCH --partition=q8
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --array=1-4096%10
#SBATCH --output=./slurm_logs/%A_%a.out

echo "========= Job started  at `date` =========="

echo "My jobid: $SLURM_JOB_ID"
echo "My array id: $SLURM_ARRAY_TASK_ID"

models_per_job=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
s=$(echo "20000 + $models_per_job * ($SLURM_ARRAY_TASK_ID - 1)" | bc -l)

#python3 sobol_dispatcher.py -s $s -N $models_per_job -d AGSS09 \
#    -M 0.6 1 -Y 0.2 0.34 -Z 0.001 0.02 -a 0.5 3 -t 6 13.799 -b -0.1 0.1 -r -A 6

module load anaconda3

python3 sobol_dispatcher.py -s $s -N $models_per_job -d GS98f -r

echo "========= Job finished at `date` =========="
#
