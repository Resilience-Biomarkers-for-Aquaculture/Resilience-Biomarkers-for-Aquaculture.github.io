---
layout: post
title: SW-ES Meeting Summaries
tags: meetings
---

Summary of each meeting and priority lists based on discussions

## 2024-10-23 SW-ES

We discussed the output of the traditional vs pseudo-alignment RNAseq workflow with a reference genome (Emma's output). Pseudo-alignment is faster and more efficient than traditional alignment, but they accomplish the same task. I.e. if we need to search for "DNA" in a book chapter, pseudo alignment removes any sentences with "D" before looking for "DNA" because it's impossible to find "DNA" without "D". So this is less computationally intensive. The nextflow RNAseq gives more quality assessments steps for traditional alignment and since we aren't strapped for time, we are ignoring the pseudo alignment path for now.

The reference based RNAseq workflow was successful! ~ 12 hours.

To-Do:  
- Emma: Re-run RNAseq with only traditional alignment and fastp quality arguments     
- Emma: Download metadata from dataset 1 publication to find the heat-stressed vs. ambient information   
- Shelly: Create de novo reference genome (to then be fed into Emma's RNAseq workfow)     
- Shelly: Gather list of disease resistance datasets to use this workflow on  

![](https://github.com/Resilience-Biomarkers-for-Aquaculture/Resilience-Biomarkers-for-Aquaculture.github.io/blob/master/img/IMG_2683.JPG)

## 2024-09-11
Data set: Roberto's thermotolerance RNAseq data

Analyses:
- Original: HiSAT2 -> stringtie
- Our attempt: Salmon
- De novo: RNAspades

Potential intern to help with analysis:
- Funding through Northeast intern through NSF (Emma is familiar with this program) or through MLSC
- Create a job posting for intern to gauge interest

To Do:
- Shelly: get Emma access to Gannet
- Shelly: Run FetchNGS and denovotranscript on Klone
- Shelly: Find data

## 2024-08-29

Do different pipelines (original RNAseq analysis, reference based RNAseq, or denovo genome based RNAseq) or annotations (recent vs. old) lead to more biomarkers? What is gained, what is lost? Are there major differences?

To Do:
- Decide on data set to start with
