---
layout: post
title: Compare genes after CGI id conversion attempt 1
author: Shelly Wanamaker
tags: name-conversion identifiers genome-versions
---

In order to get 1:1 matching between genes in Roberto's gene counts matrix and Emma's gene counts matrix, they need to be the same type of identifier. Most genes in Emma's are "LOC" and Roberto's are "CGI"

This is discussed here:[https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_thermo_RNAseq/issues/1](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_thermo_RNAseq/issues/1)

In this attempt, I'm using the genome comparison files from NCBI.

See Jupyter Notebook here: [CGI_id_conversion_using_compare_prev_docs.ipynb](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_thermo_RNAseq/blob/main/analyses/2024-12-05_SW_HeatOysterRNAseq/CGI_id_conversion_using_compare_prev_docs.ipynb)

I first checked how counts compared across samples from using different pipelines, which didn't require matching gene IDs. Overall, Roberto found ~60K genes because his pipeline included novel transcripts and Emma's found ~30K genes because her's didn't include novel transcripts. The distrinbution of counts are nearly the same with Emma's data showing more genes with counts in the 100s range and Roberto's data showing more genes with counts in the 10s range. Below are the density plots showing this. I also made some [histograms](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_thermo_RNAseq/blob/main/analyses/2024-12-05_SW_HeatOysterRNAseq/histogram.jpg) but they aren't normalized to % genes.

 [![](https://raw.githubusercontent.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_thermo_RNAseq/refs/heads/main/analyses/2024-12-05_SW_HeatOysterRNAseq/density_plots.jpg)](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_thermo_RNAseq/blob/main/analyses/2024-12-05_SW_HeatOysterRNAseq/density_plots.jpg)

After the 1:1 matching, I got ~2K CGI genes that had LOC IDs in the new genomes and made some scatter plots. In general genes with less than 10 reads show little correlation, but the correlation improves with more reads.

All the data plotted together:
[![](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_thermo_RNAseq/blob/main/analyses/2024-12-05_SW_HeatOysterRNAseq/corplot_2K.jpg?raw=true)](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_thermo_RNAseq/blob/main/analyses/2024-12-05_SW_HeatOysterRNAseq/corplot_2K.jpg)

Data facetted by sample:
[![](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_thermo_RNAseq/blob/main/analyses/2024-12-05_SW_HeatOysterRNAseq/corplot_2K_eachSample.jpg?raw=true)](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_thermo_RNAseq/blob/main/analyses/2024-12-05_SW_HeatOysterRNAseq/corplot_2K_eachSample.jpg)

There were 45,360 observations for the overlapping genes where both methods show 0 counts (31,540), both methods showed non-zero counts (3959), nf-core detected counts and Roberto's method didn't (9,065), and where Roberto's method detected counts and nf-core didn't (796).

Here is the scatter plot filtered for the observation where both methods showed non-zero:
[![](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_thermo_RNAseq/blob/main/analyses/2024-12-05_SW_HeatOysterRNAseq/scatter_nonzeros.jpg?raw=true)](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_thermo_RNAseq/blob/main/analyses/2024-12-05_SW_HeatOysterRNAseq/scatter_nonzeros.jpg)

Here's the same plot but binned and colored with the number of observations indicated in lighter blue:
[![](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_thermo_RNAseq/blob/main/analyses/2024-12-05_SW_HeatOysterRNAseq/hexbin_nonzeros.jpg?raw=true)](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_thermo_RNAseq/blob/main/analyses/2024-12-05_SW_HeatOysterRNAseq/hexbin_nonzeros.jpg)

For observations where counts were only detected by either nf-core or Roberto's methods, the majority of these had very low counts.
[![](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_thermo_RNAseq/blob/main/analyses/2024-12-05_SW_HeatOysterRNAseq/hist_zeroVSnonzeros.jpg?raw=true)](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_thermo_RNAseq/blob/main/analyses/2024-12-05_SW_HeatOysterRNAseq/hist_zeroVSnonzeros.jpg)

As an additional check, I subsetted out the genes from Roberto's supplementary table 6 which were genes that showed high expression in heat resistant oysters. I subsetted these genes in Emma's gene count matrix, then plotted them to see if they still show the same trend.

[![](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_thermo_RNAseq/blob/main/analyses/2024-12-05_SW_HeatOysterRNAseq/heat_resistant_high_exp_genes_boxplots.jpg?raw=true)](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_thermo_RNAseq/blob/main/analyses/2024-12-05_SW_HeatOysterRNAseq/heat_resistant_high_exp_genes_boxplots.jpg)

It is not obvious to me...

I will attemt a different method of gene ID matching that involves BLAST instead of ID:ID matching.

R analysis here: [https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_thermo_RNAseq/blob/main/analyses/2024-12-05_SW_HeatOysterRNAseq/2024-12-05_GeneCountsCmpr.Rmd](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cgigas_thermo_RNAseq/blob/main/analyses/2024-12-05_SW_HeatOysterRNAseq/2024-12-05_GeneCountsCmpr.Rmd)
