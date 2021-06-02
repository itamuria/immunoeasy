#!/bin/bash
## Initial SBATCH commands
#SBATCH --job-name=multiqc
#SBATCH --mail-type=END
#SBATCH --mail-user=itamayou@unav.es
#SBATCH --time=23:59:00
#SBATCH --output=Multiqc.log
#SBATCH --mem=2G
#SBATCH --cpus-per-task=2
#SBATCH -p short

date "+Hasi: %Y-%m-%d %H:%M:%S"

module load FastQC/0.11.8-Java-1.8
module load  MultiQC/1.7-foss-2018b-Python-2.7.15

cd /home/itamayou/Hepamut/Hepamut_5xPac

multiqc .

cd /home/itamayou/Hepamut/Hepamut_7xPac

multiqc .
