---
layout: post
title: Re-run FetchNGS for Cvirg_Perkinsus_RNAseq
author: Shelly Wanamaker
tags: nextflow Klone FetchNGS Perkinsus Dermo
---

I realized I made a couple mistakes in getting the data downloaded for the Perkinsus RNAseq reanalysis.

1) I had included the bisulfite data in the meta data linked with Johnson et al. in the ['datasets processed table'](https://resilience-biomarkers-for-aquaculture.github.io/about/)
  - corrected meta data: [SraRunTable (1).csv](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/data/SraRunTable%20(1).csv)
2) I had linked to Roberto's data instead of the Chan et al. data
  - corrected meta data: [SraRunTable (2).csv](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/data/SraRunTable%20(2).csv)

I corrected these links in the [table](https://resilience-biomarkers-for-aquaculture.github.io/about/) and recreated the [ids.csv](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/data/ids.csv) so it now has the correct runs. There are 155 total.
1) 61 from Proestou et al
2) 40 from Johnson et al
3) 10 from Chan et al.
4) 44 from Proestou et al.

I then ran the following:
```
# create a screen session
screen -S nextflow

# request a compute node (mem and time requests can be modified)
salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=16GB --time=12:00:00

# load the nextflow environment
mamba activate nextflow
nextflow run \
nf-core/fetchngs \
-resume \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20241205_FetchNGS/ids.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20241205_FetchNGS \
--download_method sratools
```

output will be `/gscratch/scrubbed/strigg/analyses/20241205_FetchNGS/``
