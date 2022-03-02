#!/bin/bash
#SBATCH -N 1
#SBATCH -p short
#SBATCH --ntasks-per-node=40
#SBATCH -t 24:00:00
#SBATCH -J IPP197
#SBATCH -o OUT/IPP197.out
#SBATCH -e ERR/IPP197.err
#SBATCH --mail-user=sathesh@iastate.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
scontrol show job $SLURM_JOB_ID
module purge
cellranger count --id=IPP1 \
                   --transcriptome=ssc97 \
                   --fastqs=Pig-1-IPP \
                   --sample=Pig-1-IPP-1,Pig-1-IPP-2,Pig-1-IPP-3,Pig-1-IPP-4 \
		   --chemistry=SC3Pv2 \
                   --localcores=38
