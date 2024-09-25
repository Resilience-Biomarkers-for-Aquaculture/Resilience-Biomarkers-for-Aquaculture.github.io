---
layout: post
title: Run denovotranscript on Klone
tags: nextflow Klone denovotranscript
---

This entry describes running denovotranscript on the UW HPC Klone server. I generally followed the [example workflow](https://nf-co.re/denovotranscript/dev/docs/usage/example_workflow) provided in the [denovotranscript usage documentation](https://nf-co.re/denovotranscript/dev/docs/usage).

#### 1. Run the pipeline with --qc_only with default params to check the quality of your reads.

```
# create a screen session
screen -S nextflow

# request a compute node (mem and time requests can be modified)
salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=16GB --time=12:00:00

# load the nextflow environment
mamba activate nextflow

# run nextflow pipeline
nextflow run \
nf-core/denovotranscript \
-resume \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20240925/samplesheet/samplesheet.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20240925_denovo \
--qc_only

```

This completed successfully
![](https://raw.githubusercontent.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_denovotranscript/master/analyses/20240925/Screenshot%202024-09-25%20164327.png)

I looked at the [multiqc](https://htmlpreview.github.io/?https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_denovotranscript/blob/main/analyses/20240925/multiqc_report_qc_only.html) and decided the first 10 bases should be trimmed from all reads. Example below
![](https://raw.githubusercontent.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_denovotranscript/master/analyses/20240925/Screenshot%202024-09-25%20175800.png)

2. Run the pipeline with qc_only and any custom parameters that you have decided to use based on your data. Use resume to avoid unnecessarily rerunning unchanged steps. I next ran:
```bash
nextflow run nf-core/denovotranscript \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20240925/samplesheet/samplesheet.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20240925_denovo \
--qc_only
--extra_fastp_args='--trim_front1 10 --trim_front2 10' \
--remove_ribo_rna \
-resume
```
this failed.
![](https://raw.githubusercontent.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_denovotranscript/master/analyses/20240925/Screenshot%202024-09-25%20182046.png)
I closed the screen with ctrl + A + D and then screen -XS nextflow quit

I decided to just try the whole pipeline and ran:
```bash
screen -S nextflow

salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=100GB --time=2-12:00:00

mamba activate nextflow

nextflow run nf-core/denovotranscript \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20240925/samplesheet/samplesheet.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20240925_denovo \
--extra_fastp_args='--trim_front1 10 --trim_front2 10' \
--remove_ribo_rna \
--busco_lineage= 'mollusca_odb10' \
-resume

```

next time I would like to run the nextflow code with the following parameters
-with-report nf_report
-with-trace
-with-timeline nf_timeline

----
