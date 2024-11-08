---
layout: post
title: Run RNAseq on perkinsus datasets
tags: nextflow Klone RNAseq
---

This is a meta analysis of 4 RNAseq datasets each comparing different oysters infected with Perkinsus marinus.

I first got all the meta data from the SRA Run Selector (e.g. https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP301630)
I downloaded this csv file by clicking the 'Metadata' button under the download heading.

I then got lists of accessions and concatenated them into one list ([https://raw.githubusercontent.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/refs/heads/main/data/ids.csv?token=GHSAT0AAAAAACZ6576TNOOXXECEX3GPQQ3WZZJH34Q](https://raw.githubusercontent.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/refs/heads/main/data/ids.csv?token=GHSAT0AAAAAACZ6576TNOOXXECEX3GPQQ3WZZJH34Q)) to be fed to the FetchNGS pipeine.

Then I ran FetchNGS to get all the datasets

### FetchNGS

```
# create a screen session
screen -S nextflow

# request a compute node (mem and time requests can be modified)
salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=16GB --time=12:00:00

# load the nextflow environment
mamba activate nextflow

# run nextflow pipeline
nextflow run \
nf-core/fetchngs \
-resume \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20241104_RNAseq/ids.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20241104_RNAseq \
--download_method sratools

```

### genome for C. virginica
[https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=6565](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=6565)

### Next steps
1. Run fastqc on reads to see if any trimming needs to be done
2. Run Emma's RNAseq pipeline
