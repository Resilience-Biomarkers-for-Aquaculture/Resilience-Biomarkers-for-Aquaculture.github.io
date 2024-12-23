---
layout: post
title: Run Seqera rnaseq pipeline for C. gigas study
tags: nextflow seqera aws rnaseq
---

## 12-23-2024

The goal here was to compare gene counts between `rnaseq` runs using
1. [reference genome and gene annotations from the period of Roberto's study](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000297895.1/) (herein referred to as `2012`)
2. most recent [reference genome and gene annotations](https://ncbi.nlm.nih.gov/datasets/genome/GCA_963853765.1/) (herein referred to as `2023`) for C. gigas.

As previously [noted](https://resilience-biomarkers-for-aquaculture.github.io/SW-CGI_ID_matching_attemp_1/), the
gene count matrix that Roberto provided overwhelmingly used `MSTRG_` (novel or unannotated genes/transcripts identified by
StringTie) and `CGI_` (CpG islands) gene IDs from the annotation
that he used, whereas the avaialbe NCBI annotations for `2023` uses `LOC_` gene IDs
commonly used for annotated genes or predicted genes that do not yet have an official gene symbol or name.
However, we found that the `2012` NCBI genome and annotations did use `LOC_` gene IDs, which facilitates matching by ID
rather than requiring BLAST or other indirect gene ID matching methods. So, for this current task, we used the `2012`
reference genome as as a proxy for Roberto's, given that it represents the state of the art at the time of his
publication.

