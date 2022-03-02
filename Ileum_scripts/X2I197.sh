#!/bin/bash
#SBATCH -N 1
#SBATCH -p short
#SBATCH --ntasks-per-node=40
#SBATCH -t 24:00:00
#SBATCH -J X2I197
#SBATCH -o OUT/X2I197.out
#SBATCH -e ERR/X2I197.err
#SBATCH --mail-user=sathesh@iastate.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
scontrol show job $SLURM_JOB_ID
module purge
cellranger count --id=X2I1 \
                   --transcriptome=ssc97 \
                   --fastqs=Pig-1-2l \
                   --sample=Pig-1-2I-1,Pig-1-2I-2,Pig-1-2I-3,Pig-1-2I-4 \
		   --chemistry=SC3Pv2 \
                   --localcores=38
