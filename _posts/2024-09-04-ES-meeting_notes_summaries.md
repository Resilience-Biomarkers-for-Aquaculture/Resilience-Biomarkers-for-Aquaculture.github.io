---
layout: post
title: SW-ES Meeting Summaries
tags: meetings
---

Summary of each meeting and priority lists based on discussions

## 2024-11-20 SW-ES 

Steve Yost will be joining as a volunteer to help with bioinformatic workflows and github project management. Shelly and Emma served on a panel for ocean and coastal acidification put on by the New England aquarium and Benchmark Strategies. 

Since last time, Shelly downloaded the 4 disease dataset and Emma began running this set but with the rnaseq pipeline. Currently stuck on issues with the gtf file. We figured out that my issue is probably only the blank gene ID issue instead of the '+'/'.' issue. Shelly ran the below code and figured out that there are only '+' and '-' characters in the strand column (column #7).

```
(base) [strigg@klone-login03 C_virginica]$ zcat GCF_002022765.2_C_virginica-3.0_genomic.gtf.gz | awk -F"\t" '{print  $7}
' | sort | uniq -c
      5
 761431 +
 770111 -
```

To address the blank gene IDs issue, we changed the one gene with several exons to read gene_id="unknown_transcript_1"

```
## check that changes worked 
check for the changes: zcat GCF_002022765.2_C_virginica-3.0_genomic.gtf.gz | awk -F"\t" '{if($9 ~/gene_id ""/ ) gsub(/gene_id ""/,"gene_id \"unknown_transcript_1\"",$9);print $0}' | grep "unknown_t" | less

## save as new gtf file 
zcat GCF_002022765.2_C_virginica-3.0_genomic.gtf.gz | awk -F"\t" '{if($9 ~/gene_id ""/ ) gsub(/gene_id ""/,"gene_id \"unknown_transcript_1\"",$9);print $0}' | gzip > mod_GCF_002022765.2_C_virginica-3.0_genomic.gtf.gz
```

Do we need to run all four together or separately? If we want to see what consistently comes out as biomarkers, I think we need to do this separately? Do we get the same answer when analyzing altogether? Try both. 

Shelly has been trying to get the de novo pipeline working - which we think worked! The de novo assembly worked but the QC steps failed so Shelly will troubleshoot those but we can start using the assembly as input for the rnaseq workflow. 

We found a differential abundance workflow that will take rnaseq output and try quantitative comaprisons (e.g., DESEQ2). [Pipeline here](https://nf-co.re/differentialabundance/1.5.0/). 

Emma:  
1. Run disease dataset with modified GTF file.    
2. Run 4 separately and then altogether    
3. Tagseq vs rnaseq pipeline details - does tagseq need anything else?   
4. Identify counts table to use for meta comparison for Cgigas comparison (gene vs. transcript?)   
5. Attempt differential abundance workflow on counts matrix identified in #4   
6. Move Roberto data output to Gannet  

Shelly:  
1. Troubleshoot denovo QC steps   
2. Run rnaseq with denovo transcriptome   
3. Find counts table from Roberto's data 



## 2024-11-08 SW-ES 

De novo still running into errors. rnaseq pipeline worked well! Moving on to meta analysis of 4 oyster immunity datasets - Shelly ran [fetchNGS](https://resilience-biomarkers-for-aquaculture.github.io/SW-fetchNGS_Cvig_Prkns/) with one csv with all 4 dataset SRA ids. 2 files within this pipeline didn't work and now it's stalled. We checked the current output folder, there were 398 files but should be 362 total.. Where are these extra coming from?

We figured out that we don't have access to the rest of ckpt, just srlab. The config file Shelly made allows for use of other nodes besides srlab. 

We decided to just start fresh because it might be an issue with the -resume flag. There are 221 samples from the 4 datasets with 3 paired-end reads and 1 project with tag seq which is single end reads. Shelly will start running this today to let it download over the weekend. We slack'd Sam and Steven about the ckpt node on Hyak to see if we can get access. Steven and Sam have run jobs on cpkt successfully. 

Emma:  
1. Check Roberto's dataset with meta information on heat treatments and compare to our rnaseq output  
2. Wait for Shelly to finish fetchNGS for 4 disease datasets  (done 11-09-2024)  

Shelly:  
1. Finish fetchNGS for 4 disease datasets   (done 11-09-2024)  
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

