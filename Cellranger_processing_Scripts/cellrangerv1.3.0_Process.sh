#!/bin/bash

#################### Install Cellranger [v1.3.0]
#################### Install bcl2fastq  [v2.18.0]
#################### Install samtools   [v1.8.0]
#################### cellranger.csv is required to run cellranger
#################### PATH to transcriptome files is required to run cellranger
#################### PATH-FASTQ path to fastq files in required to run cellranger

cellranger mkfastq --run=$DIR/bcl-files --samplesheet=cellranger.csv --delete-undetermined

cellranger count --id=$SAMPLENAME --transcriptome=PATH --fastqs=PATH-FASTQ --localmem 132