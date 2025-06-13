---
layout: post
title: Compare fastp parameters with TAG-seq data
author: Shelly Wanamaker
tags: TAG-seq RNA-seq fastp adapter
---

code from 03-27-25

```
# create a screen session
screen -S nextflow

# request a compute node (mem and time requests can be modified)
salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=16GB --time=12:00:00

# load the nextflow environment
apptainer shell --bind /usr/bin/ /gscratch/ /gscratch/srlab/containers/srlab-R4.4-bioinformatics-container-a83e518.sif 





nextflow run \
nf-core/fetchngs \
-resume \
-c /gscratch/srlab/nextflow/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20250327_FetchNGS/ids.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250327_FetchNGS \
--download_method sratools


installed miniforge in nextflow bin

/gscratch/srlab/nextflow/bin/miniforge/bin/mamba create -p miniforge/envs/nextflow -c bioconda -c conda-forge python=3.8 nf-core nextflow -y -q

/gscratch/srlab/nextflow/bin/miniforge/bin/mamba init


#activate srlab nextflow env
/gscratch/srlab/nextflow/bin/miniforge/bin/mamba activate /gscratch/srlab/nextflow/bin/miniforge/envs/nextflow




#variables don't seem set an i don't know where NXF_HOME or NXF_TEMP are. so I'll just run the pipeline and see what happens

NXF_HOME='/gscratch/srlab/programs/nextflow'
```


If you don't have mamba in your path, run
`/gscratch/srlab/nextflow/bin/miniforge/bin/mamba init`

then 
`source ~/.bashrc`

# create a screen session
screen -S nextflow

# request a compute node (mem and time requests can be modified)
salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=16GB --time=12:00:00


nextflow run nf-core/fetchngs -resume -c /gscratch/srlab/nextflow/uw_hyak_srlab.config --input /gscratch/scrubbed/strigg/analyses/20250327_FetchNGS/ids.csv --outdir output --download_method sratools




# load the nextflow environment
#apptainer shell --bind /gscratch/ srlab-R4.4-bioinformatics-container-0751826.sif

# run nextflow pipeline
nextflow run \
nf-core/<pipeline> \
-resume \
-c /gscratch/srlab/nextflow/uw_hyak_srlab.config \
<pipeline specific flags>


nextflow run nf-core/fetchngs -resume -c /gscratch/srlab/nextflow/uw_hyak_srlab.config --input /gscratch/scrubbed/strigg/analyses/20250327_FetchNGS/ids.csv --outdir /gscratch/scrubbed/strigg/analyses/20250327_FetchNGS --download_method sratools