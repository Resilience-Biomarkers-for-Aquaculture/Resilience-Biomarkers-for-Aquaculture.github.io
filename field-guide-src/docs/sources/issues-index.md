# GitHub Issues Index

This index organizes the GitHub issues from the [Cvirg_Pmarinus_RNAseq repository](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues) that were referenced in the project timeline and field guide. Issues are grouped thematically for easier navigation.

All issue links point to: `https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/[number]`

---

## Data Preparation & Integration

### [Issue #3](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/3): Created Merged Metadata

**Timeline:** January 2025  
**Status:** Completed  
**Summary:** Initial effort to create unified metadata across multiple RNA-seq datasets for integrated analysis.

### [Issue #9](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/9): Add Study 5

**Timeline:** January 2025  
**Status:** Completed  
**Summary:** Incorporated an additional dataset (Study 5) into the analysis pipeline to increase sample size and assess reproducibility.

### [Issue #54](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/54): Find More Datasets

**Timeline:** September 2025  
**Status:** Postponed  
**Summary:** Exploration of additional RNA-seq datasets for broader validation. Deferred to focus on existing dataset analysis.

---

## Differential Abundance Analysis

### [Issue #4](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/4): Differential Abundance Initial Approach

**Timeline:** January 2025  
**Status:** Completed  
**Summary:** First attempt at differential abundance analysis. Toy example succeeded, but encountered GTF file issues and weak trait separation in full datasets.

**Key finding:** Mutual information analysis on Perkinsus datasets showed insufficient separation between tolerant/sensitive groups.

### [Issue #12](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/12): Decide Best Methods and Execute

**Timeline:** January 2025  
**Status:** Completed  
**Summary:** Decision point for selecting optimal differential abundance methodologies after initial exploration revealed data integration challenges.

### [Issue #29](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/29): Differential Abundance on All Datasets Together

**Timeline:** April 2025  
**Status:** Deferred  
**Summary:** Ran differential abundance on all datasets together but results were not yet interpreted. Could return to this analysis.

### [Issue #31](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/31): Interpret Combined Dataset Results

**Timeline:** April 2025  
**Status:** Deferred  
**Related:** [Analysis results](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/31)  
**Summary:** Follow-up to #29; interpretation was postponed to focus on per-dataset approaches.

### [Issue #32](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/32): Run Differential Abundance on Datasets Separately

**Timeline:** June 2025  
**Status:** Completed  
**Summary:** Pivot to per-dataset analysis. Steve focused on Study 5, Shelly focused on Study 1. Goal was understanding optimal parameters for differential abundance pipeline.

### [Issue #36](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/36): Compare DEGs Across Datasets

**Timeline:** September 2025  
**Status:** Not completed  
**Summary:** Compare differential abundance results run independently for each dataset. Theme: post-data integration approach. Question remains about reproducibility vs. integrated data analysis, but subsetting approach is uncertain.

### [Issue #46](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/46): Integrate All Data Through Differential Abundance Pipeline

**Timeline:** September 2025  
**Status:** Not completed  
**Summary:** Another attempt at integrated data analysis. Could revisit.

---

## Batch Effects & Normalization

### [Issue #18](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/18): Batch Correction Attempts

**Timeline:** February 2025  
**Status:** Completed  
**Summary:** Attempted batch correction using COMBAT and RemoveBatchEffect methods. Results showed little improvement in trait-based separation.

**Key lesson:** Study-specific effects were stronger than trait effects; batch correction couldn't recover sufficient signal.

### [Issue #34](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/34): Combine Study 4 Injected + Study 5

**Timeline:** July 2025  
**Status:** Completed  
**Summary:** Experimental combination of compatible studies (Study 4 injected group + Study 5) to increase sample size.

**Research question:** Would Study 4 injected samples cluster with resistant or susceptible phenotype from Study 5?

**Learning:** Gained understanding of normalization timing in differentialabundance pipeline (PCAs before analysis, normalization during). Started seeing evidence of innate trait.

**Could return to:** Revisit analysis on 567 significant DEGs (by DESeq) to see if clustering improves compared to top 500 most variable genes.

---

## Technical Issues & Parameters

### [Issue #26](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/26): Parameter Selection & GC Bias

**Timeline:** April 2025  
**Status:** Completed  
**Summary:** Identified that Johnson dataset used TAG-seq (not standard RNA-seq) and discovered GC bias. Determined that initial analysis parameters were inappropriate for TAG-seq data.

### [Issue #28](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/28): Rerun Johnson Data with Different Parameters

**Timeline:** April 2025  
**Status:** Completed  
**Related:** [Notebook post](https://resilience-biomarkers-for-aquaculture.github.io/SW-FastPparams4tagseq/)  
**Summary:** Reprocessed TAG-seq data with appropriate FastP parameters to address issues identified in #26.

---

## GSEA & Pathway Analysis

### [Issue #41](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/41): Stepwise Differential Abundance Approach

**Timeline:** August 2025  
**Status:** Completed  
**Summary:** Developed two-step approach: (1) Controls vs. treated, (2) Resistant vs. sensitive from step 1 genes.

**Implementation:** [`analyses/stepwise_differentialabundance/`](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/tree/main/analyses/stepwise_differentialabundance)

**Results on Dataset 1:** Only 1 significant gene. DESeq2 struggled with small gene set (~50 genes).

**Related:** [Stepwise notebook](https://resilience-biomarkers-for-aquaculture.github.io/SW-diffabund_stepwise_ds1/)

### [Issue #45](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/45): Understand GSEA

**Timeline:** August-September 2025  
**Status:** Not completed  
**Summary:** Goal was to better understand and apply Gene Set Enrichment Analysis (GSEA) for pathway-level interpretation. Deferred.

---

## Classifier Development & Validation

### [Issue #42](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/42): Validate SR320 Classification Results

**Timeline:** August 2025  
**Status:** In progress  
**Summary:** Validation of classification results. Question: Are the ~50 candidate markers convincing?

**Action needed:** Make plots to assess marker quality.

### [Issue #43](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/43): SR320's AI Model

**Timeline:** August 2025  
**Status:** Completed  
**Summary:** Development of machine learning classifier for phenotype prediction.

### [Issue #44](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/44): Combined Datasets 1 & 5

**Timeline:** August 2025  
**Status:** ✅ **Completed - 6-gene classifier success**  
**Summary:** Comparison of integrated data analysis vs. post-data integration approaches using datasets 1 and 5.

**Pipeline:**

- **Step 1:** Rank genes by reproducibility, directionality consistency, and heterogeneity
- **Step 2:** Logistic regression for minimal gene set

**Result:** 6-gene classifier panel with strong separation between tolerant and sensitive phenotypes.

**Key lesson:** Only include training set in test set if exploring within study; for cross-study prediction, keep training and test separate (LOSO validation).

**Related:** [Gene classifier notebook](https://resilience-biomarkers-for-aquaculture.github.io/SY-gene-classifier-panel/)

### [Issue #47](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/47): (Details Unknown)

**Timeline:** September 2025  
**Status:** Not pursued  
**Summary:** No need to revisit.

### [Issue #49](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/49): Plot 6 Genes to Gain Confidence

**Timeline:** September 2025  
**Status:** In progress  
**Summary:** Visualization of 6-gene panel performance to assess how well genes distinguish phenotypes across studies.

**Related:** [Six-gene biomarker exploration](https://resilience-biomarkers-for-aquaculture.github.io/SY-six-gene-biomarker-exploration/)

### [Issue #51](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/51): Replot Heatmap with Improved Clustering

**Timeline:** September 2025  
**Status:** To revisit  
**Summary:** Improve heatmap visualizations with better clustering algorithms and labels for clearer interpretation of gene panel performance.

### [Issue #52](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/52): Coverage Density Plots

**Timeline:** September 2025  
**Status:** Needs notebook entry  
**Summary:** Generate coverage density plots for quality control and validation. Still requires documentation in a notebook post.

### [Issue #53](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/53): Innate vs. Reactive Gene Expression

**Timeline:** September 2025  
**Status:** Completed  
**Summary:** Critical analysis determining whether biomarkers are innate (constitutively different in controls) or reactive (induced by stress).

**Key insight:** Biomarkers may exist in control groups if they represent innate resilience traits. The stepwise approach may remove these.

**Related:** [Innate gene expression notebook](https://resilience-biomarkers-for-aquaculture.github.io/SY-innate-gene-expression/)

---

## Literature Comparison

### [Issue #39](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/39): Compare DE Results from Papers

**Timeline:** July 2025  
**Status:** Remaining to be done  
**Summary:** Systematic comparison of project DEG results with published literature on oyster stress response. Create consolidated list of known DEGs/markers for validation.

---

## Navigation

**Main Guide:**

- [Start Here](../index.md) - Field guide overview
- [Timeline](../timeline.md) - Detailed project history
- [Pipelines](../pipelines/decision-tree.md) - Analysis workflows

**Other Sources:**

- [Notebook Posts Index](posts-index.md) - Analysis notebooks organized by theme
