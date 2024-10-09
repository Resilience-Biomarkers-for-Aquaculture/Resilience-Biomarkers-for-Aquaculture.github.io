---
layout: post
title: Running RNAseq pipeline with reference transcript
tags: RNAseq nextflow Klone referencetranscript
---

## RNASeq workflow 

This workflow is fully documented in a separate post [here](https://resilience-biomarkers-for-aquaculture.github.io/ES-RNAseq-pipeline-construction-post/).

Datasets:    
- [1: Arredondo-Espinoza et al 2023 C. gigas thermotolerance](#dataset-1-arredondo-espinoza-et-al-2023)    


Data is downloaded via nextflow pipeline [fetchNGS](https://nf-co.re/fetchngs/1.12.0/) and see Shelly's example download post [here](https://resilience-biomarkers-for-aquaculture.github.io/a-fetchNGSKlone/). 

Creating a mamba environment on my user for following workflows:

```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh

bash Mambaforge-Linux-x86_64.sh -b -p $HOME/mambaforge

source ~/mambaforge/bin/activate
mamba init 
## close and reopen shell for changes to take place 

mamba create -n nextflow -c bioconda -c conda-forge python=3.8 nf-core nextflow -y -q
```

## Dataset 1: Arredondo-Espinoza et al 2023

Paper found [here](https://www.sciencedirect.com/science/article/pii/S1744117X23000345?via%3Dihub). Shelly downloaded this data on [20240925 with fetchNGS](https://resilience-biomarkers-for-aquaculture.github.io/a-fetchNGSKlone/). 

This study focused on the Pacific oyster, *Crassostrea gigas*, and completed a laboratory-simulated daily oscillatory thermal challenge (26 to 34 °C) for 30 days. 

Path to fastq files: `/mmfs1/gscratch/scrubbed/<ShellynetID>/analyses/20240925/fastq`   
Path to samplesheet created by fetchNGS: `/mmfs1/gscratch/scrubbed/<ShellynetID>/analyses/20240925/samplesheet/samplesheet.csv` 

Path to Emma's analysis: `/mmfs1/gscratch/scrubbed/<EmmanetID>/Cgigas_ArredondoEspinoza2023/`

Testing this pipeline: 

```
# create screen session
screen -S rnaseq_cgigas_AE2023

# request a compute node (mem and time requests can be modified)
salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=40GB --time=12:00:00

# activate mamba environment (created at beginning of this post)
mamba activate nextflow

# update rnaseq workflow in nextflow within mamba environment
nextflow pull nf-core/rnaseq

## attach back to screen 
tmux attach-session -t rnaseq_cgigas_AE2023
```

