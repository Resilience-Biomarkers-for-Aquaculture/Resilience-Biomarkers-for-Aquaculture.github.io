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
- [nextflow report](https://htmlpreview.github.io/?https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_denovotranscript/blob/main/analyses/20241009/attempt01/nf_report.html)

![](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_denovotranscript/blob/main/analyses/20241009/attempt01/Screenshot%202024-10-09%20100350.png)

I reran the pipeline specifying to only use RNAspades as an assembler:
```
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
```

This completed with errors
- [nextflow report](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cgigas_denovotranscript/20241009/nf_report.html)

There seemed to be an issue with BUSCO and RNAQUAST
```
Command error:
  INFO:    Environment variable SINGULARITYENV_TMPDIR is set, but APPTAINERENV_TMPDIR is preferred
  INFO:    Environment variable SINGULARITYENV_NXF_TASK_WORKDIR is set, but APPTAINERENV_NXF_TASK_WORKDIR is preferred
  INFO:    Environment variable SINGULARITYENV_NXF_DEBUG is set, but APPTAINERENV_NXF_DEBUG is preferred
  usage: busco -i [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] -m [MODE] [OTHER OPTIONS]
  busco: error: argument -l/--lineage_dataset: expected one argument
  ```
I had tried to specify mollusca_odb10 as the lineage dataset to use but it seems to not have liked that. So I tried to remove that parameter and just have it auto-detect which is the [default](https://nf-co.re/denovotranscript/1.0.0/parameters/)

### run on 10-11-2024
```
nextflow run nf-core/denovotranscript \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20240925/samplesheet/samplesheet_test.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20241011_denovo \
--extra_fastp_args='--trim_front1 10 --trim_front2 10' \
--remove_ribo_rna \
--assemblers rnaspades \
-resume \
-with-report nf_report.html \
-with-trace \
-with-timeline nf_timeline.html
```
