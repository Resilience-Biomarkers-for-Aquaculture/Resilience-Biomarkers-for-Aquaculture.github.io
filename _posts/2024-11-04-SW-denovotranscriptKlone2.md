---
layout: post
title: Run denovotranscript on Klone attempt 3
tags: nextflow Klone denovotranscript
---

### ran on 11-04-2024
```
nextflow run nf-core/denovotranscript \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
--input /gscratch/scrubbed/strigg/analyses/20240925/samplesheet/samplesheet.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20241104_denovo \
--extra_fastp_args='--trim_front1 10 --trim_front2 10' \
--remove_ribo_rna \
--assemblers rnaspades \
-resume \
-with-report nf_report.html \
-with-trace \
-with-timeline nf_timeline.html
```

This errored out at SORTMERNA step and couldn't figure out why. Here are some screenshots of what was happening:
[![](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_denovotranscript/blob/main/analyses/20241104/Screenshot%202024-11-04%20150221.png)](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_denovotranscript/blob/main/analyses/20241104/Screenshot%202024-11-04%20150221.png)
[![](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_denovotranscript/blob/main/analyses/20241104/Screenshot%202024-11-04%20150236.png)](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_denovotranscript/blob/main/analyses/20241104/Screenshot%202024-11-04%20150236.png)


I reran the pipeline in a new screen session and node allocation.
