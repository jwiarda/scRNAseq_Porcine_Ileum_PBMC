#!/bin/bash
#SBATCH -N 1
#SBATCH -p short
#SBATCH --ntasks-per-node=40
#SBATCH -t 24:00:00
#SBATCH -J PBMC197
#SBATCH -o OUT/PBMC197.out
#SBATCH -e ERR/PBMC197.err
#SBATCH --mail-user=sathesh@iastate.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
scontrol show job $SLURM_JOB_ID
module purge
cellranger count --id=PBMC1 \
                   --transcriptome=ssc97 \
                   --fastqs=Pig-1-PBMC \
                   --sample=Pig-1-PBMC-1,Pig-1-PBMC-2,Pig-1-PBMC-3,Pig-1-PBMC-4 \
		   --chemistry=SC3Pv2 \
                   --localcores=38
