---
layout: post
title: PCA for gene count comparison between 2012 and 2023 reference genomes
tags: rnaseq pca 
---

## 12-23-2024

The goal here was to compare gene counts between `rnaseq` runs using
1. [The NCBI reference genome and gene annotations from the period of Roberto's study](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000297895.1/) (herein referred to as `2012`)
2. The most recent [reference genome and gene annotations](https://ncbi.nlm.nih.gov/datasets/genome/GCA_963853765.1/) (herein referred to as `2023`) for C. gigas.

## Motivation
As previously [noted](https://resilience-biomarkers-for-aquaculture.github.io/SW-CGI_ID_matching_attemp_1/), in attempint to compare gene counts between
Roberto's results and `rnaseq` results, we encountered a hurdle:
The gene count matrix that Roberto provided overwhelmingly used gene IDs prefixed with `MSTRG_` (novel or unannotated genes/transcripts identified by
StringTie) and `CGI_` (CpG islands), whereas the avaialbe NCBI annotations, and thus gene count data, for `2023` uses gene IDs prefixed with `LOC_`,
commonly used for annotated genes or predicted genes that do not yet have an official gene symbol or name.
However, we found that the `2012` NCBI genome and annotations did use `LOC_` gene IDs, which facilitates matching by ID
rather than requiring BLAST or other indirect gene ID matching methods.

So, for this current task, we used the `2012`
reference genome (linked above) as as a proxy for Roberto's, given that it represents the state of the art at the time of his
publication.



