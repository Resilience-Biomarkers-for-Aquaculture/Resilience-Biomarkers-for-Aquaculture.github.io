---
layout: post
title: Run denovotranscript on Klone attempt 2
tags: nextflow Klone denovotranscript
---


The 9/25/24 attempt failed
![](https://raw.githubusercontent.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_denovotranscript/master/analyses/20240925/Screenshot%202024-10-01%142013.png)

I can't find the log because I'm not sure where I ran the code from.

I'm going to try to just process one sample all the way through the pipeline

#create samplesheet with just one sample
head -2 samplesheet.csv > samplesheet_test.csv

```

screen -S nextflow

salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=100GB --time=2-12:00:00

mamba activate nextflow

nextflow run nf-core/denovotranscript \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20240925/samplesheet/samplesheet_test.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20241009_denovo \
--extra_fastp_args='--trim_front1 10 --trim_front2 10' \
--remove_ribo_rna \
--busco_lineage= 'mollusca_odb10' \
-resume \
-with-report nf_report.html \
-with-trace \
-with-timeline nf_timeline.html

```
This failed at Trinity step
- see log file


nextflow run nf-core/denovotranscript \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20240925/samplesheet/samplesheet_test.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20241009_denovo \
--extra_fastp_args='--trim_front1 10 --trim_front2 10' \
--remove_ribo_rna \
--assemblers rnaspades \
--busco_lineage= 'mollusca_odb10' \
-resume \
-with-report nf_report.html \
-with-trace \
-with-timeline nf_timeline.html
