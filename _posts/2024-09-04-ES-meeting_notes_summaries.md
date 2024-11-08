---
layout: post
title: SW-ES Meeting Summaries
tags: meetings
---

Summary of each meeting and priority lists based on discussions

## 2024-11-08 SW-ES 

De novo still running into errors. rnaseq pipeline worked well! Moving on to meta analysis of 4 oyster immunity datasets - Shelly ran [fetchNGS](https://resilience-biomarkers-for-aquaculture.github.io/SW-fetchNGS_Cvig_Prkns/) with one csv with all 4 dataset SRA ids. 2 files within this pipeline didn't work and now it's stalled. We checked the current output folder, there were 398 files but should be 362 total.. Where are these extra coming from?

We figured out that we don't have access to the rest of ckpt, just srlab. The config file Shelly made allows for use of other nodes besides srlab. 

We decided to just start fresh because it might be an issue with the -resume flag. There are 221 samples from the 4 datasets with 3 paired-end reads and 1 project with tag seq which is single end reads. Shelly will start running this today to let it download over the weekend. We slack'd Sam and Steven about the ckpt node on Hyak to see if we can get access. Steven and Sam have run jobs on cpkt successfully. 

Emma:  
1. Check Roberto's dataset with meta information on heat treatments and compare to our rnaseq output  
2. Wait for Shelly to finish fetchNGS for 4 disease datasets  

Shelly:  
1. Finish fetchNGS for 4 disease datasets   
2. Continue to troubleshoot de novo pipeline  


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

