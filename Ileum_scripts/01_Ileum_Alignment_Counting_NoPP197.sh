#!/bin/bash
#SBATCH -N 1
#SBATCH -p short
#SBATCH --ntasks-per-node=40
#SBATCH -t 24:00:00
#SBATCH -J NoPP197
#SBATCH -o OUT/NoPP197.out
#SBATCH -e ERR/NoPP197.err
#SBATCH --mail-user=sathesh@iastate.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
scontrol show job $SLURM_JOB_ID
module purge
cellranger count --id=NoPP1 \
                   --transcriptome=ssc97 \
                   --fastqs=Pig1-No-PP \
                   --sample=Pig1-No-PP-1,Pig1-No-PP-2,Pig1-No-PP-3,Pig1-No-PP-4 \
		   --chemistry=SC3Pv2 \
                   --localcores=38
