---
layout: post
title: Step-wise differential abundance with Cvirg dataset1
author: Shelly Wanamaker
tags: differentialabundance RNA-seq nf-core 
---

## Purpose

This analysis of dataset 1 ([Proestou et al. 2023](https://doi.org/10.3389/fgene.2023.1054558)) is part of the post-analysis data integration approach, where RNAseq datasets are processed independently and the [differential abundance pipeline](https://nf-co.re/differentialabundance) is run on each to identify DEGs. The lists of DEGs will then be compared across studies. 

This study excluded controls from their main analysis and focused on family differences in dosed animals. Because there isn't currently a way to account for controls in multi-factor model, I wondered if I could first run the differential abundance pipeline on infection vs. control, then on breed to look for family differences. 

## Methods

**Step 1: contrast controls vs. injected**

input data:

- [study1_contrasts.csv](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250730_diffabund/ds_1/study1_contrasts.csv)

```
id,variable,reference,target
trait_dose,condition,control,injected
```

- [study1_samplesheet.csv](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250730_diffabund/ds_1/study1_samplesheet.csv)
- [salmon.merged.gene_counts_length_scaled_dataset1.tsv](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250730_diffabund/ds_1/salmon.merged.gene_counts_length_scaled_dataset1.tsv)

run pipeline: 

```
screen -S dfab_1
salloc -A srlab -p cpu-g2-mem2x -N 1 -c 1 --mem=16GB --time=12:00:00
mamba activate nextflow 

nextflow run nf-core/differentialabundance \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
-resume \
--study_name dataset_1 \
--input study1_samplesheet.csv \
--contrasts study1_contrasts.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250730_diffabund \
--matrix salmon.merged.gene_counts_length_scaled_dataset1.tsv \
--feature_type gene \
--features_gtf_feature_type gene \
--filtering_min_abundance 10 \
--filtering_min_proportion 0.75 \
--gtf /gscratch/srlab/strigg/GENOMES/GCF_002022765.2_C_virginica-3.0_genomic_noEmptyGeneIDs.gtf


#this attempt didn't have the threshold so I reran as the code above
	
nextflow run nf-core/differentialabundance \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
-resume \
--study_name dataset_1 \
--input study1_samplesheet.csv \
--contrasts study1_contrasts.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250730_diffabund/ds_1 \
--matrix salmon.merged.gene_counts_length_scaled_dataset1.tsv \
--feature_type gene \
--features_gtf_feature_type gene \
--gtf /gscratch/srlab/strigg/GENOMES/GCF_002022765.2_C_virginica-3.0_genomic_noEmptyGeneIDs.gtf
```

report: [dataset_1.html](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250730_diffabund/ds_1/report/dataset_1.html)

output: [20250730_diffabund/ds_1](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250730_diffabund/ds_1)

**Step 2: contrast families of injected animals using the DEGs from the step 1 contrast**

input files: 

- [trait_dose.deseq2.results_filtered.tsv](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250730_diffabund/ds_1/tables/differential/trait_dose.deseq2.results_filtered.tsv)
- contrast file for sensitive vs. resistant[study1_contrasts.csv](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250702_diffabund/study1_contrasts.csv)

```
id,variable,reference,target
trait_dose,trait,sensitive,resistant
```

- filtered gene counts matrix for DEGs from step 1 
	- [salmon.merged.gene_counts_length_scaled_dataset1_ctrlFilt.tsv](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250730_diffabund/ds_1/salmon.merged.gene_counts_length_scaled_dataset1_ctrlDEseq2filt.tsv)
	- code to create filtered matrix:

```
#filter the gene counts matrix for DEGs from step 1
awk -F "\t" 'NR==FNR{a[$1]=$1;next}$1 in a{print $0}' tables/differential/trait_dose.deseq2.results_filtered.tsv salmon.merged.gene_counts_length_scaled_dataset1.tsv > salmon.merged.gene_counts_length_scaled_dataset1_ctrlFilt.tsv
```

- samplesheet:
	- wget https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250702_diffabund/study1_samplesheet.csv

Run pipeline on filtered data:

```
nextflow run nf-core/differentialabundance \
-c /gscratch/srlab/strigg/bin/uw_hyak_srlab.config \
-resume --study_name dataset_1 \
--input ctrl_filt/study1_samplesheet.csv \
--contrasts ctrl_filt/study1_contrasts.csv \
--outdir /gscratch/scrubbed/strigg/analyses/20250730_diffabund/ds_1/ctrl_filt \
--matrix ctrl_filt/salmon.merged.gene_counts_length_scaled_dataset1_ctrlFilt.tsv \
--feature_type gene --features_gtf_feature_type gene --filtering_min_abundance 10 \
--filtering_min_proportion 0.75 --gtf /gscratch/srlab/strigg/GENOMES/GCF_002022765.2_C_virginica-3.0_genomic_noEmptyGeneIDs.gtf

```

Copy data to Gannet
```
rsync --progress --verbose --archive 20250730_diffabund shellytrigg@gannet.fish.washington.edu:/volume2/web/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq
```

## Results

- report: [dataset_1.html](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250730_diffabund/ds_1/ctrl_filt/report/dataset_1.html)

- output: [20250730_diffabund/ds_1/ctrl_filt](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250730_diffabund/ds_1/ctrl_filt)

**PCA of DEseq2 DEGs**

- filter counts for DEGs
	- use VST counts 
		- [vst_ctrlDEseq2filt.tsv](https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/Cvirg_RNAseq/20250730_diffabund/ds_1/vst_ctrlDEseq2filt.tsv)
		- code: 
```
awk -F"\t" 'NR==FNR{a[$1]=$1;next}$1 in a{print $0}' ctrl_filt/tables/differential/trait_dose.deseq2.results_filtered.tsv tables/processed_abundance/all.vst.tsv > vst_ctrlDEseq2filt.tsv
```

- plot PCA
	- code here: [https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/scripts/20250806_dataset1_stepwise_diffAbund.ipynb](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/scripts/20250806_dataset1_stepwise_diffAbund.ipynb)
	
The PCA shows a good separation of the sensitive vs. tolerant animals using the 26 significant DEGs identified in step 2, which makes the results believable.

[![](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/analyses/DiffAbundStepwise_pca_ds1.png?raw=true)](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/analyses/DiffAbundStepwise_pca_ds1.png?raw=true)