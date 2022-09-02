#!/bin/bash
## Initial SBATCH commands
#SBATCH --job-name=bcl2fastq
#SBATCH --mail-type=END
#SBATCH --mail-user=itamayou@unav.es
#SBATCH --time=23:59:00
#SBATCH --output=Logs/bcl2fastq.log
#SBATCH --mem=10G
#SBATCH --cpus-per-task=10
#SBATCH -p short


## Run script

date "+Hasi: %Y-%m-%d %H:%M:%S"

module load bcl2fastq2

bcl2fastq -r 10 \
	-p 20 \
	--no-lane-splitting \
	--minimum-trimmed-read-length 60 \
	--mask-short-adapter-read 10 \
	--create-fastq-for-index-reads \
	-w 10

