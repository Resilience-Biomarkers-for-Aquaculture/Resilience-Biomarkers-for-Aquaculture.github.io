# Glossary

Key terms and concepts used throughout this field guide, defined in the context of the biomarker discovery project.

---

## A

**Batch Effects**  
Non-biological variation in gene expression data arising from technical differences between experimental batches (e.g., different sequencing runs, library preparation dates, labs). In this project, study-specific effects consistently dominated trait effects, limiting the effectiveness of integrated analysis approaches.

**Biomarker**  
A measurable molecular indicator (in this case, gene expression level) that distinguishes between biological states (tolerant vs. sensitive phenotypes). The project identified a 6-gene biomarker panel.

## C

**Classifier**  
A machine learning model that predicts categorical outcomes (phenotypes) based on input features (gene expression). This project used logistic regression with LASSO regularization to build a minimal 6-gene classifier.

**COMBAT**  
A batch correction method that adjusts for systematic technical variation using empirical Bayes. Attempted in [Issue #18](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/18) with limited success.

**Cross-Validation**  
A technique for assessing how well a model generalizes to independent data by systematically holding out portions of data for testing. See: [LOSO](#loso-leave-one-study-out).

## D

**DEG (Differentially Expressed Gene)**  
A gene whose expression level differs significantly between two or more conditions (e.g., resistant vs. sensitive, control vs. treated).

**DESeq2**  
A widely-used R package for differential gene expression analysis from RNA-seq count data. Used throughout this project for identifying DEGs.

**Differential Abundance**  
Analysis method to identify features (genes, transcripts) that differ in quantity/expression between experimental groups. This project primarily used the [nf-core/differentialabundance](https://nf-co.re/differentialabundance) pipeline.

**Directionality Consistency**  
In the two-step classifier approach, the requirement that a gene's fold change direction (up or down) is the same across multiple datasets. Ensures reproducibility of effect direction.

## F

**FastP**  
A tool for quality control and preprocessing of FASTQ files (raw sequencing reads). Parameter selection is critical for different library types (e.g., TAG-seq vs. standard RNA-seq) - see [Issue #26, #28](timeline.md#april-2025-per-dataset-analysis--tag-seq-issues).

**Fold Change**  
The ratio of gene expression levels between two conditions, often expressed as log2 fold change. Positive values indicate upregulation, negative values indicate downregulation.

## G

**GMT File**  
Gene Matrix Transposed file format used for gene set definitions in GSEA. Contains pathway/gene set names with associated genes. Created in this project for [GSEA integration](timeline.md#gsea-integration).

**GSEA (Gene Set Enrichment Analysis)**  
A computational method that determines whether a defined set of genes (e.g., a pathway) shows statistically significant differences between two biological states. Attempted in [Issue #45](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/45).

## H

**Heterogeneity**  
In the context of gene scoring, the variance in gene expression within phenotype groups. Low heterogeneity (high within-group consistency) is desirable for biomarkers.

## I

**Innate Biomarker**  
A gene that is constitutively different in expression between resistant and sensitive individuals, even in the absence of stress. Contrasts with [reactive biomarkers](#reactive-biomarker). Key insight from [Issue #53](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/53).

**Integrated Data Analysis**  
Approach that pools multiple datasets together before analysis, treating all samples as if from a single study (with batch correction). In this project, integrated analysis **failed** due to strong study-specific effects ([Big Lesson #1](timeline.md#january-2025-differential-abundance-beginnings)).

## L

**LASSO (Least Absolute Shrinkage and Selection Operator)**  
A regularization method for regression that can shrink coefficients to exactly zero, performing automatic feature selection. Used in the two-step classifier to identify the minimal gene set.

**Leakage (Data Leakage)**  
When information from the test set inappropriately influences the training process, leading to overly optimistic performance estimates. See [Validation & Pitfalls](pipelines/validation.md#1-data-leakage).

**LOSO (Leave-One-Study-Out)**  
A cross-validation strategy where each study is held out once as a test set, while all other studies form the training set. Critical for assessing cross-study generalization ([Big Lesson #4](timeline.md#august-2025-gsea--stepwise-approach-development)).

## M

**Meta-Analysis**  
Statistical approach that combines results from multiple independent studies to identify consistent patterns. This project used a post-data integration meta-analysis approach.

**Mutual Information**  
A measure of statistical dependence between variables. Used early in the project to assess gene-phenotype associations, but trait separation was weak.

## N

**nf-core**  
A community effort to collect curated bioinformatics pipelines built using Nextflow. This project used [nf-core/rnaseq](https://nf-co.re/rnaseq) and [nf-core/differentialabundance](https://nf-co.re/differentialabundance).

**Normalization**  
Process of adjusting gene expression data to account for technical variation (library size, sequencing depth, GC content) while preserving biological variation. Critical for cross-sample comparison.

## P

**PCA (Principal Component Analysis)**  
A dimensionality reduction technique that identifies major sources of variation in high-dimensional data (like gene expression). Used throughout the project to visualize sample clustering and assess batch effects.

**Perkinsus marinus**  
Protozoan parasite that causes Dermo disease in oysters, the primary stressor in this project's datasets.

**Phenotype**  
Observable trait or characteristic. In this project: tolerant/resistant vs. sensitive/susceptible oyster phenotypes in response to *P. marinus* infection.

**Post-Data Integration**  
Approach that analyzes each dataset independently, then compares results across datasets to identify reproducible findings. This project pivoted to this approach after integrated analysis failed.

## R

**Reactive Biomarker**  
A gene whose differential expression between resistant and sensitive individuals emerges only after stress exposure (i.e., significant in treated samples but not controls). Contrasts with [innate biomarkers](#innate-biomarker).

**RemoveBatchEffect**  
A function from the limma R package for adjusting gene expression data to remove batch effects. Attempted in [Issue #18](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/18) with limited success.

**Reproducibility**  
In the context of the two-step classifier, the degree to which a gene is identified as differentially expressed across multiple independent datasets. High reproducibility indicates robust, generalizable biomarkers.

**RNA-seq**  
RNA sequencing - a high-throughput method for quantifying gene expression by sequencing cellular RNA. The primary data type in this project.

## S

**Stepwise Differential Abundance**  
A two-step filtering approach: (1) identify stress-responsive genes (control vs. treated), then (2) identify resistance-associated genes from Step 1 genes. Has limitations (removes innate biomarkers, breaks down with small gene sets) - see [Stepwise Pipeline](pipelines/decision-tree.md).

**Study-Specific Effects**  
Systematic differences between datasets due to experimental design, protocols, or biological differences between populations. In this project, these effects were **stronger than trait effects**, driving the pivot to post-data integration.

## T

**TAG-seq (3' Tag RNA-Sequencing)**  
A cost-effective alternative to standard RNA-seq that sequences only the 3' end of transcripts. Requires different data processing parameters than full-length RNA-seq (see [TAG-seq issues](timeline.md#april-2025-per-dataset-analysis--tag-seq-issues)).

**Tolerant/Resistant**  
Phenotype category for oysters that survived *P. marinus* infection or showed low infection intensity and minimal pathology.

**Trait Effects**  
Gene expression differences attributable to the biological phenotype of interest (resistant vs. sensitive). In this project, trait effects were initially weaker than [study-specific effects](#study-specific-effects).

## V

**Validation**  
Process of confirming that a biomarker panel or model performs well on independent data not used during development. See [Validation & Pitfalls](pipelines/validation.md).

**VST (Variance-Stabilizing Transformation)**  
A DESeq2 normalization method that transforms count data to a scale where variance is roughly constant across the range of expression values. Enables visualization and clustering. Requires sufficient genes (>1000) for robust estimation.

---

## Abbreviations

- **DE:** Differential Expression
- **DEG:** Differentially Expressed Gene
- **FDR:** False Discovery Rate (adjusted p-value threshold)
- **GC Bias:** Systematic bias related to GC content in sequencing
- **GTF:** Gene Transfer Format (gene annotation file)
- **HPC:** High-Performance Computing
- **LOSO:** Leave-One-Study-Out
- **PCA:** Principal Component Analysis
- **QC:** Quality Control
- **WGBS:** Whole-Genome Bisulfite Sequencing (for methylation)

---

## Species

- ***Crassostrea virginica*:** Eastern oyster (primary study species)
- ***Crassostrea gigas*:** Pacific oyster (used for reference genome/annotation)
- ***Perkinsus marinus*:** Dermo disease-causing parasite
- ***Mytilus chilensis*:** Chilean blue mussel (methylation studies)

---

**Related Pages:**

- [Start Here](index.md) - Field guide overview
- [Problem Framing](problem-framing.md) - Detailed background on resilience and datasets
- [Timeline](timeline.md) - See concepts in chronological context
