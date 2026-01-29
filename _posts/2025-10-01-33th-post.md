---
layout: post
title: Series of DESeq2 runs indicates innate DEGs for tolerance
author: Steve Yost
tags: differentialabundance RNA-seq
---

## Abstract

Combining data from two compatible studies for meta-analysis, we ran
a series of exploratory DESeq2 runs, using variations of contrasts and sample sets. We found that differential expression of genes between pre-identified tolerant and sensitive samples was relatively independent of treatment (controls vs treated samples), indicating a possible innate DEG profile for tolerance.

## Discussion

Among the five studies that we're currently working with, [Study 1](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2023.1054558/full) and [Study 5](https://www.sciencedirect.com/science/article/abs/pii/S1050464819311295?via%3Dihub) are the most compatible in terms of experimental design. Both use pre-identified tolerant and sensitive families, have Control and Treatment groups, sample gill tissue, and use similar preparation and analysis methodology. So at this juncture we chose to pursue our meta-analysis route with this more constrained set. Each of the two studies included samples taken at Day 7 after infection, and that seemed to be a point where differential expression was well manifested, so we used only Day 7 samples from each study.

We were working intially with the notion that the presence of Control and Treatment groups was important, allowing us to examine primarly the _differences_ between Controls and treated samples across studies.  
We ran a series of [nf-core/differentialabundance] runs for exploration.  
However, intial contrasts for Treatment showed only one significantly differentially expressed gene! A second run, using only treated samples and contrasting for Condition also showed only one DEG. This led to a hunch that tolerance/sensitivity is shown more in innate differential expression (i.e. common to control and treated samples) between families, rather than in response to treatment.  
Further [differentialabundance]([nf-core/differentialabundance]) runs, first on each of the two studies separately, and then again on the combined set from both studies, gives a promising indication that the differential abundance between sensitive and tolerant groups is exhibited in both control and treated samples.

## Details
All runs were performed on the Seqera platform.

Seqera run `condescending_davinci_5` ([config](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/data/differentialabundance_configs/condescending_davinci_5.config)) combined Study 1 and Study 5 samples and contrasted for Treatment.  
[Sample sheet](https://raw.githubusercontent.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/refs/heads/main/data/differential_abundance_sheets/rnaseq_diffabundance_study1and5_samplesheet_filled.csv)  
[Contrasts file](https://raw.githubusercontent.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/refs/heads/main/data/differential_abundance_sheets/rnaseq_diffabundance_study1and5_D7_contrasts.csv") Contrast on Treatment, block for batch.  
[Matrix](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/raw/refs/heads/main/data/rnaseq_gene_counts/merged_gene_counts_gene_name_study1_study5_common.tsv)  
Results are in s3://steveyost-seqera/Cvirg_Pmarinus/deseq2/study1and5_D7_results
PCA showed separation in Condition in PC2, little separation in Treatment, and some in Batch in PC3 and PC4.  
Only two genes were significantly up-regulated for Treatment; none down-regulated.  
This led to the hunch that tolerance/sensitivity is shown more in innate expression rather than response to treatment (infection).

Seqera run `disturbed_fermat_2` ([config](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/data/differentialabundance_configs/disturbed_fermat_2.config)) combined Study 1 and Study 5 samples, **but only treated samples** (no controls) and contrasted for Condition.  
[Sample sheet](https://raw.githubusercontent.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/refs/heads/main/data/differential_abundance_sheets/rnaseq_diffabundance_study1and5_samplesheet_no_controls.csv)  
[Contrasts file](https://raw.githubusercontent.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/refs/heads/main/data/differential_abundance_sheets/rnaseq_diffabundance_study1and5_D7_condition_contrasts.csv) Contrast on condition, block for batch.   
[Matrix](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/raw/refs/heads/main/data/rnaseq_gene_counts/merged_gene_counts_gene_name_study1_study5_common.tsv)  
Results are in s3://steveyost-seqera/Cvirg_Pmarinus/deseq2/study1and5_D7_no_controls_results
Only one gene was found to be differentially expressed for Condition among these treated samples.

Seqera run `crazy_newton` ([config](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/data/differentialabundance_configs/crazy_newton.config)) used Study 5 only, including Control and Treatment samples, but contrasting only on Condition, with no blocking.  
[Sample sheet](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/data/differential_abundance_sheets/rnaseq_diffabundance_study5_D7_samplesheet.csv)  
[Contrasts file](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/data/differential_abundance_sheets/rnaseq_diffabundance_D7_condition_nobatch_contrasts.csv)  
Results are in s3://steveyost-seqera/Cvirg_Pmarinus/deseq2/study5_D7_condition_results  
PCA showed clear separation between Condition tolerant and sensitive in PC1, some separation of Batch in PC2 and little separation of Treatment.  
1158 genes were significantly up-regulated for Condition, 1099 down-regulated.

Seqera run `mad_mayer_3` ([config](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/data/differentialabundance_configs/mad_mayer_3.config)) is the same as `crazy_newton`, but for Study 1 only.   
PCA showed strong separation in Treatment in PC1 and PC2, and good separation in Condition in PC2 and PC3.  
1448 genes were significantly up-regulated for Condition, 1313 down-regulated.

Finally, Seqera run `Issue_44_study1and5_condition` ([config](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/data/differentialabundance_configs/Issue_44_study1and5_condition.config)) used all samples (Control and Treatment) from Study 1 and Study 5, contrasting on Condition, blocking for Batch.  
[Sample sheet](https://raw.githubusercontent.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/refs/heads/main/data/differential_abundance_sheets/rnaseq_diffabundance_study1and5_samplesheet_filled.csv) ( same as `condescending_davinci_5`)  
[Contrasts file](https://raw.githubusercontent.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/refs/heads/main/data/differential_abundance_sheets/rnaseq_diffabundance_study1and5_D7_condition_contrasts.csv) (Same as `disturbed_fermat_2`)  
Note it used filtering_min_samples: 10
Results are in s3://steveyost-seqera/Cvirg_Pmarinus/deseq2/Issue_44_study1and5_D7_condition_results
PCA showed strong batch separation in PC1 (21.3%) (note this is pre-batch-correction in the pipeline), separation for Treatment in PC1, and separation for Condition in PC3 (5.5%)  
1046 genes were significantly up-regulated for Condition, 905 down-regulated.  
Download the [report HTML page](s3://steveyost-seqera/Cvirg_Pmarinus/deseq2/Issue_44_study1and5_D7_condition_results/report/study.html).