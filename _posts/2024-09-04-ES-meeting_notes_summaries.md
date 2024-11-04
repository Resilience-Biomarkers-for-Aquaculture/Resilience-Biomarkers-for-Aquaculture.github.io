---
layout: post
title: SW-ES Meeting Summaries
tags: meetings
---

Summary of each meeting and priority lists based on discussions 

## 2024-10-23 SW-ES

We discussed the output of the tradiitional vs pseudo-alignment rnaseq workflow with a reference genome (Emma's output). Pseudo-alignment is faster and more efficient than traditional alignment, but they accomplish the same task. I.e. if we need to search for "DNA" in a book chapter, pseudo alignment removes any sentences with "D" before looking for "DNA" because it's impossible to find "DNA" without "D". So this is less computationally intensive. The nextflow rnaseq gives more quality assessments steps for tranditional alignment and since we aren't strapped for time, we are ignoring the pseudo alignment path for now. 

The reference based rnaseq workflow was successful! ~ 12 hours. 

To-Do:  
- Emma: Re-run rnaseq with only traditional alignment and fastp quality arguments     
- Emma: Download metadata from dataset 1 publication to find the heat-stressed vs. ambient information   
- Shelly: Create de novo reference genome (to then be fed into Emma's rnaseq workfow)     
- Shelly: Gather list of disease resistance datasets to use this workflow on  

