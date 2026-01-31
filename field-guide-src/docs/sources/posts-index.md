# Notebook Posts Index

This index organizes notebook posts from the [project analysis notebook](https://resilience-biomarkers-for-aquaculture.github.io/notebook/) that document the biomarker discovery journey. Posts are grouped by theme and linked chronologically within each theme.

All notebook posts are available at: `https://resilience-biomarkers-for-aquaculture.github.io/[post-slug]/`

---

## Pipeline Development & Infrastructure

### [2024-09-10: RNAseq Workflow with Reference Genome (ES)](/ES-RNAseq_with_reference_dataset1/)

**Timeline:** September 2024  
**Summary:** Initial RNAseq workflow setup using reference genome for *C. gigas* Dataset #1. Established baseline pipeline approach.

### [2024-09-11: RNAseq Workflow with De Novo Transcriptome](/RNAseq-deNovo/)

**Timeline:** September 2024  
**Summary:** Alternative approach using de novo transcriptome assembly for species without well-annotated reference genomes.

### [2024-09-25: Set Up Nextflow on UW Klone and Run FetchNGS](/a-fetchNGSKlone/)

**Timeline:** September 2024  
**Summary:** Infrastructure setup for running nf-core pipelines on UW Klone HPC cluster. FetchNGS workflow for downloading public data.

### [2024-11-23: AWS Batch Compute Environment on Seqera (SY)](/SY-Sequera_AWS_rnaseq_test/)

**Timeline:** November 2024  
**Summary:** Cloud-based compute infrastructure setup using Seqera platform and AWS Batch for scalable RNA-seq processing.

### [2024-12-03: Repoint nf-core projectDir to srlab Directory (SW)](/SW-nextflow_mods/)

**Timeline:** December 2024  
**Summary:** Infrastructure configuration adjustments for consistent nf-core pipeline execution.

### [2025-01-06: Creating a Nextflow Pipeline for Gene Count Analysis (SY)](/SY-explore-creating-nextflow-pipeline/)

**Timeline:** January 2025  
**Summary:** Exploration of custom Nextflow pipeline development for gene count-based analyses beyond standard nf-core workflows.

---

## Differential Abundance Analysis

### [2024-12-11: Exploring nf-core Differential Abundance Workflow (ES)](/ES-differential_abundance_workflow/)

**Timeline:** December 2024  
**Related:** [Timeline - Dec 4, 2024](../timeline.md#week-of-dec-4)  
**Summary:** Initial exploration of nf-core/differentialabundance pipeline. Encountered compute resource limits when attempting to process 4 combined datasets.

**Key lesson:** Integration of heterogeneous datasets requires careful consideration of batch effects and nuances.

### [2025-07-01: Run nf-core Differentialabundance Pipeline on C. virg Data (SW)](/SW-diff-abund/)

**Timeline:** July 2025  
**Related:** [Issue #34](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/34), [Timeline - July 2025](../timeline.md#july-2025-combining-studies--understanding-normalization)  
**Summary:** Per-dataset differential abundance analysis. Combined Study 4 (injected group) + Study 5 to increase sample size.

**Key learning:** Understanding normalization timing in the pipeline:

- PCAs generated before differential abundance analysis
- Normalization (VST, rlog) happens during DESeq2
- Batch correction applied post-normalization

**Observation:** Started seeing evidence of innate trait in batch-corrected PCAs on top 500 most variable genes.

### [2025-08-05: Creating a GMT File for Use with GSEA (SY)](/SY-augment-diffexp-with-dsea/)

**Timeline:** August 2025  
**Related:** [Issue #41](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/41)  
**Summary:** Integrating Gene Set Enrichment Analysis (GSEA) with differential abundance pipeline. Created GMT file with gene descriptions for pathway-level interpretation.

### [2025-08-26: Step-wise Differential Abundance with Cvirg Dataset 1 (SW)](/SW-diffabund_stepwise_ds1/)

**Timeline:** August 2025  
**Related:** [Issue #41](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/41), [Stepwise Pipeline](../pipelines/decision-tree.md)  
**Summary:** Implementation of two-step filtering approach:

1. Step 1: Controls vs. treated (identify stress-responsive genes)
2. Step 2: Resistant vs. sensitive (from Step 1 genes)

**Results:** Only 1 significant gene identified in Step 2

**Problem:** DESeq2's VST couldn't work properly with highly filtered gene set (~50 genes)

**Key insight:** Raised question about whether stepwise filtering removes innate biomarkers present in control samples.

---

## Data Processing & Quality Control

### [2024-11-04: Run RNAseq on Perkinsus Datasets (SW)](/SW-fetchNGS_Cvig_Prkns/)

**Timeline:** November 2024  
**Summary:** Processing and quality control of *Perkinsus marinus* challenge datasets.

### [2024-12-13: Compare Genes After CGI ID Conversion Attempt 1 (SW)](/SW-CGI_ID_matching_attemp_1/)

**Timeline:** December 2024  
**Summary:** Gene identifier harmonization across datasets with different annotation versions. Critical for cross-dataset comparisons.

### [2024-12-23: PCA for Gene Count Comparison Between NCBI Annotation Releases (SY)](/SY-PCA_rnaseq_ref_genomes_annots/)

**Timeline:** December 2024  
**Summary:** Assessed impact of different *C. gigas* annotation releases on gene count matrices and downstream analysis using PCA.

### [2025-04-15: Reprocess TAG-seq Data with Different FastP Parameters (SW)](/SW-FastPparams4tagseq/)

**Timeline:** April 2025  
**Related:** [Issue #26](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/26), [Issue #28](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/28), [Timeline - April 2025](../timeline.md#tag-seq-parameter-problems)  
**Summary:** **Critical technical discovery:** Johnson dataset used TAG-seq (not standard RNA-seq). Initial analysis with standard RNA-seq parameters was inappropriate.

**Issues:**

- GC bias detected
- Wrong adapter trimming parameters
- Inappropriate quality filtering

**Solution:** Reran with TAG-seq-specific FastP parameters, improving data quality.

---

## Gene Classifier & Biomarker Development

### [2025-09-11: Two-Script Pipeline for Gene-Expression Classifier (SY)](/SY-gene-classifier-panel/)

**Timeline:** August-September 2025  
**Status:** ✅ **Core methodology for 6-gene classifier**  
**Related:** [Issue #44](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/44), [Two-Step Classifier Pipeline](../pipelines/classifier-path.md)  
**Summary:** **Landmark analysis** - Developed validated two-step classifier approach combining datasets 1 & 5.

**Methodology:**

- **Step 1 (R script):** Meta-analysis ranking genes by reproducibility, effect consistency, and heterogeneity across studies
- **Step 2 (Python script):** Logistic regression with LASSO regularization for minimal gene panel; LOSO validation

**Result:** **6-gene classifier** with strong separation between tolerant/sensitive phenotypes

**Key innovation:** Explicitly separates feature discovery from classifier construction to reduce overfitting across batches

**Code:** [Two-step pipeline repository](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/tree/main/analyses/Study1and5ThreeWay/two_step_gene_expression_classifier)

**ChatGPT session:** [Development discussion](https://chatgpt.com/share/68bc4d13-e418-800d-9960-f42e37d9f98b)

### [2025-10-01: Series of DESeq2 Runs Indicates Innate DEGs for Tolerance (SY)](/SY-innate-gene-expression/)

**Timeline:** September 2025  
**Related:** [Issue #53](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/53), [Validation & Pitfalls](../pipelines/validation.md#6-innate-vs-reactive-biomarkers)  
**Summary:** **Critical biological insight** - Discovered that differential expression between tolerant/sensitive samples is relatively independent of treatment status.

**Finding:** Biomarkers appear to be **innate** (constitutively different between families) rather than reactive (induced by stress)

**Implication:** Stepwise filtering approach (control vs. treated first) would remove these innate biomarkers

**Approach:**

- Used Day 7 samples from Studies 1 & 5 (most compatible designs)
- Ran contrasts on: Treatment, Treated-only samples, Control-only samples, Combined samples
- Treatment contrast: only 1 DEG
- Treated-only contrast: only 1 DEG
- Control and treated samples showed similar tolerant vs. sensitive DEGs

**Conclusion:** Tolerance/sensitivity is an innate differential expression profile, common to control and treated samples

### [2025-10-01: Exploring Six-Gene Biomarker Across Studies (SY)](/SY-six-gene-biomarker-exploration/)

**Timeline:** September 2025  
**Related:** [Issue #49](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/49)  
**Summary:** Visualization and validation of the 6-gene classifier panel across multiple studies. Assesses how well genes distinguish phenotypes in independent datasets.

### [2026-01-05: Common Genes per LOSO Fold (SY)](/SY-plot-DEGs-per-fold/)

**Timeline:** January 2026  
**Summary:** Analysis of gene consistency across Leave-One-Study-Out (LOSO) cross-validation folds. Identifies genes that are robustly selected regardless of which study is held out.

---

## Exploratory & Supporting Analyses

### [2024-12-31: Using ChatGPT to Explore Thermal Resilience Based on Gene Counts (SY)](/SY-ChatGPT_abinitio_thermal_tolerance/)

**Timeline:** December 2024  
**Summary:** Exploratory use of AI tools (ChatGPT) for hypothesis generation about thermal tolerance based on gene expression patterns. Documents AI-assisted analysis workflow.

---

## Methylation Sequencing

### [2024-11-27: Run Methylseq on Klone (SW)](/SW-methylseq_Cvirg/)

**Timeline:** November 2024  
**Summary:** Epigenetic analysis using nf-core/methylseq pipeline on *C. virginica* whole-genome bisulfite sequencing (WGBS) data.

### [2025-04-07: Reproduce Oyster WGBS Analysis with nf-core Methylseq (SW)](/SW-CompareSrLab2Methylseq/)

**Timeline:** April 2025  
**Summary:** Validation of methylseq pipeline results by comparing with published oyster WGBS analysis.

### [2025-05-02: Run C. virginica WGBS Data Through nf-core Methylseq (SW)](/SW-methylseq_ceasmallr/)

**Timeline:** May 2025  
**Summary:** Processing additional *C. virginica* WGBS datasets for methylation biomarker discovery.

### [2025-08-06: Run Methylseq on M. chilensis WGBS Data (SW)](/SW-methylseq_Mchilensis/)

**Timeline:** August 2025  
**Summary:** Expansion to Chilean blue mussel (*Mytilus chilensis*) methylation data to assess cross-species epigenetic patterns.

---

## Thematic Groups Summary

### 🔧 **Pipeline Development** (6 posts)

Posts focused on infrastructure, workflow setup, and computational resource configuration.

### 📊 **Differential Abundance** (4 posts)

Posts documenting differential expression analysis attempts, including integrated, per-dataset, and stepwise approaches.

### 🎯 **Classifier & Biomarkers** (4 posts)

Posts detailing the development and validation of the 6-gene classifier, including innate vs. reactive characterization.

### 🔬 **Data QC & Processing** (4 posts)

Posts addressing data quality, annotation harmonization, and technical parameter optimization (especially TAG-seq).

### 🧬 **Methylation** (4 posts)

Posts exploring epigenetic biomarkers using WGBS data.

---

## Chronological Navigation

**2024 Posts:** [September](/archivebydate/) | [November](/archivebydate/) | [December](/archivebydate/)

**2025 Posts:** [January](/archivebydate/) | [April](/archivebydate/) | [May](/archivebydate/) | [July-October](/archivebydate/)

**All Posts:** [Archive by Date](/archivebydate/)

---

## Navigation

**Main Guide:**

- [Start Here](../index.md) - Field guide overview
- [Timeline](../timeline.md) - Detailed project history organized by month
- [Pipelines](../pipelines/decision-tree.md) - Validated analysis workflows

**Other Sources:**

- [GitHub Issues Index](issues-index.md) - Issues organized by theme
