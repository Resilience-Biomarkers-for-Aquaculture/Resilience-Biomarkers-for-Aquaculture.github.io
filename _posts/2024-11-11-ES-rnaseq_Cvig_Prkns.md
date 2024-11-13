---
layout: post
title: Run RNAseq on perkinsus datasets
tags: nextflow Klone RNAseq
---

## RNASeq workflow for meta analysis with Cvig immunity 

All datasets work with the Eastern Oyster, *C. virginica*, so I need to download the [genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002022765.2/), [link2](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/022/765/GCF_002022765.2_C_virginica-3.0/). 

```
cd /mmfs1/gscratch/scrubbed/elstrand/genomes
mkdir C_gigas

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/022/765/GCF_002022765.2_C_virginica-3.0/GCF_002022765.2_C_virginica-3.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/022/765/GCF_002022765.2_C_virginica-3.0/GCF_002022765.2_C_virginica-3.0_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/022/765/GCF_002022765.2_C_virginica-3.0/GCF_002022765.2_C_virginica-3.0_genomic.gtf.gz
```

### Datasets

See [ABOUT](https://resilience-biomarkers-for-aquaculture.github.io/about/) page for table of datasets  

| C. virginica         | T   | disease | Perkinsus | infection tolerance | resilience  | Proestou et al. 2023       | https://doi.org/10.3389/fgene.2023.1054558        | https://www.ncbi.nlm.nih.gov/sra?LinkName=bioproject_sra_all&from_uid=894694 | SraRunTable.csv    |   |   |
|----------------------|-----|---------|-----------|---------------------|-------------|----------------------------|---------------------------------------------------|------------------------------------------------------------------------------|--------------------|---|---|
| C. virginica         | E,T | disease | Perkinsus | infected            | sensitivity | Johnson et al. 2020        | https://doi.org/10.3389/fmars.2020.00598          | https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP246310&o=acc_s%3Aa         | SraRunTable(1).csv |   |   |
| C. gigas & virginica | T   | disease | Perkinsus | infection tolerance | resilience  | Chan et al. 2021           | 10.3389/fgene.2021.795706                         | https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN11031730&o=acc_s%3Aa      | SraRunTable(2).csv |   |   |
| C. virginica         | T   | disease | Perkinsus | infection           | sensitivity | Sullivan and Proestou 2021 | https://doi.org/10.1016/j.aquaculture.2021.736831 | https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP301630&o=acc_s%3Aa         | SraRunTable(3).csv |   |   |


### Workflow 

This workflow is fully documented in a separate post [here](https://resilience-biomarkers-for-aquaculture.github.io/ES-RNAseq-pipeline-construction-post/).

Using the mamba environment I created for [CGigas thermal tolerance dataset #1](https://resilience-biomarkers-for-aquaculture.github.io/ES-RNAseq_with_reference_dataset1/). 

Path to fastq files: `/gscratch/scrubbed/strigg/analyses/20241108_RNAseq`   
Path to samplesheet created by fetchNGS: `/gscratch/scrubbed/strigg/analyses/20241108_RNAseq/samplesheet/samplesheet.csv` 

Path to Emma's analysis: `/mmfs1/gscratch/scrubbed/elstrand/Cvir_disease_2024Nov11`
Path to C.gigas genome: `/mmfs1/gscratch/scrubbed/elstrand/genomes/C_virginica/GCF_002022765.2_C_virginica-3.0_genomic.fna.gz`

Create new samplesheet to fit rnaseq workflow: 

```
scp xxx@klone.hyak.uw.edu:/gscratch/scrubbed/strigg/analyses/20241108_RNAseq/samplesheet/samplesheet.csv C:\Users\EmmaStrand\Downloads

scp C:\Users\EmmaStrand\MyProjects\Resilience_Biomarkers_Aquaculture\samplesheet_rnaseq_Cvir_disease_set1.csv xxx@klone.hyak.uw.edu:/mmfs1/gscratch/scrubbed/elstrand/Cvir_disease_2024Nov11/
```

R script on my computer to fix this sample sheet. I realized after that I over-wrote the .csv file but just moving on to uploading.

```
library(tidyverse)

fetchNGS_samplehseet <- 
  read.csv("C:/Users/EmmaStrand/MyProjects/Resilience_Biomarkers_Aquaculture/samplesheet_rnaseq_Cvir_disease_set1.csv") 

rnaseq_samplesheet <- fetchNGS_samplehseet %>%
  dplyr::select(sample, fastq_1, fastq_2) %>%
  mutate(strandedness = "auto")

rnaseq_samplesheet %>% 
  write.csv("C:/Users/EmmaStrand/MyProjects/Resilience_Biomarkers_Aquaculture/samplesheet_rnaseq_Cvir_disease_set1.csv",
            row.names = FALSE)
```


Nextflow RNAsew workflow (`01-rnaseq.sh`):

Start conda environment before running sbatch. 

```
#!/bin/bash
#SBATCH --job-name=rnaseq
#SBATCH --account=srlab
#SBATCH --error=output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=cpu-g2-mem2x
#SBATCH --nodes=1
#SBATCH --time=1-20:00:00
#SBATCH --mem=100G
#SBATCH --ntasks=7
#SBATCH --cpus-per-task=2
#SBATCH --chdir=/mmfs1/gscratch/scrubbed/elstrand/Cvir_disease_2024Nov11/scripts

#conda activate nextflow

samplesheet="/mmfs1/gscratch/scrubbed/elstrand/Cvir_disease_2024Nov11/samplesheet_rnaseq_Cvir_disease_set1.csv"
output="/mmfs1/gscratch/scrubbed/elstrand/Cvir_disease_2024Nov11/"
genome="/mmfs1/gscratch/scrubbed/elstrand/genomes/C_virginica/GCF_002022765.2_C_virginica-3.0_genomic.fna.gz"
gff="/mmfs1/gscratch/scrubbed/elstrand/genomes/C_virginica/GCF_002022765.2_C_virginica-3.0_genomic.gff.gz"
gtf="/mmfs1/gscratch/scrubbed/elstrand/genomes/C_virginica/GCF_002022765.2_C_virginica-3.0_genomic.gtf.gz"

nextflow run nf-core/rnaseq \
    -resume \
    -profile singularity \
    --input ${samplesheet} \
    --outdir ${output} \
    --gtf ${gtf} \
    --gff ${gff} \
    --fasta ${genome} \
    --trimmer fastp \
    --extra_fastp_args '--cut_mean_quality 30 --trim_front1 10 --trim_front2 10' \
    --aligner star_salmon \
    --skip_pseudo_alignment \
    --multiqc_title Cvir_disease_2024Nov11 \
    --deseq2_vst
```

### 2024-11-11 / 2024-11-13

First attempt to sbatch rnaseq.sh above with all four datasets at once. Should we also think about running these separately? UW was doing maintenance so checked output on 11-13. I hadn't activated the nextflow conda environment so I activated that (`conda activate nextflow`) and used `sbatch 01-rnaseq.sh`. This started running immediately.

