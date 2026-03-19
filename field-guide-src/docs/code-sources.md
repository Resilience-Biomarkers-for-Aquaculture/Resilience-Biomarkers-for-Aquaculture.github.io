# Analysis Code & Source Code

*Entry point for reproducibility — where to find and how to run everything.*

---

## GitHub Repository Structure

The primary analysis code lives in the [Cvirg_Pmarinus_RNAseq repository](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq):

```
Cvirg_Pmarinus_RNAseq/
├── analyses/
│   ├── stepwise_differentialabundance/   # Stepwise DE pipeline (Datasets 1 & 5)
│   └── Study1and5ThreeWay/
│       └── two_step_gene_expression_classifier/  # Two-step classifier (R + Python)
├── data/                                 # Input count matrices and metadata
└── README.md
```

The [main project repository](https://github.com/Resilience-Biomarkers-for-Aquaculture/Resilience-Biomarkers-for-Aquaculture.github.io) (this site) contains the field guide source and notebook posts.

---

## Key Analysis Scripts

### Two-Step Gene Expression Classifier

**Location:** [`analyses/Study1and5ThreeWay/two_step_gene_expression_classifier/`](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/tree/main/analyses/Study1and5ThreeWay/two_step_gene_expression_classifier)

| Script | Language | Purpose |
|---|---|---|
| Step 1 script | R | Meta-analysis gene ranking: reproducibility, directionality consistency, heterogeneity |
| Step 2 script | Python | Logistic regression with LASSO; LOSO cross-validation |

**Related notebook:** [Two-script pipeline for gene-expression classifier](https://resilience-biomarkers-for-aquaculture.github.io/SY-gene-classifier-panel/)

### Stepwise Differential Abundance

**Location:** [`analyses/stepwise_differentialabundance/`](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/tree/main/analyses/stepwise_differentialabundance)

Implements the two-step filtering approach: (1) control vs. treated, then (2) resistant vs. sensitive on the filtered gene set. See [Stepwise Pipeline](pipelines/stepwise.md) for methodology and known limitations.

**Related notebook:** [Step-wise differential abundance on Dataset 1](https://resilience-biomarkers-for-aquaculture.github.io/SW-diffabund_stepwise_ds1/)

---

## How to Run Pipelines from Scratch

### Environment Setup

**R dependencies** (for Step 1 and `nf-core/differentialabundance` post-processing):

```r
install.packages(c("DESeq2", "limma", "ggplot2", "dplyr"))
# or via Bioconductor:
BiocManager::install(c("DESeq2", "limma"))
```

**Python dependencies** (for Step 2 classifier):

```bash
pip install scikit-learn pandas numpy
```

**nf-core pipelines** (for upstream processing):

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash

# Run nf-core/rnaseq (read alignment and quantification)
nextflow run nf-core/rnaseq -profile docker --input samplesheet.csv --genome GRCv ...

# Run nf-core/differentialabundance (per-dataset DE analysis)
nextflow run nf-core/differentialabundance -profile docker ...
```

See [nf-core/rnaseq docs](https://nf-co.re/rnaseq) and [nf-core/differentialabundance docs](https://nf-co.re/differentialabundance) for full parameter documentation.

### Key Parameter Considerations

- **TAG-seq datasets** require different FastP parameters (see [Lesson 6](lessons-learned.md#lesson-6-technology-differences-require-different-parameters))
- **PCAs in `differentialabundance`** are generated before normalization — do not interpret them as post-normalization QC
- **VST normalization** requires > 1000 genes for stable estimation; the stepwise pipeline may violate this

---

## Notebook Posts by Pipeline Phase

For detailed documentation of what was run and why, see notebook posts organized by project phase:

| Phase | Relevant Notebooks |
|---|---|
| Phase 1: Initial Integration | [Differential abundance workflow exploration](https://resilience-biomarkers-for-aquaculture.github.io/2024-12-11-ES-differential_abundance_workflow/) |
| Phase 2: Batch Effects | [PCA and annotation comparison](https://resilience-biomarkers-for-aquaculture.github.io/SY-PCA_rnaseq_ref_genomes_annots/) |
| Phase 3: TAG-seq & Per-Dataset | [Reprocess TAG-seq with FastP params](https://resilience-biomarkers-for-aquaculture.github.io/SW-FastPparams4tagseq/) |
| Phase 4: Normalization | [Differential abundance on C.virg data](https://resilience-biomarkers-for-aquaculture.github.io/SW-diff-abund/) |
| Phase 5: Stepwise DE | [Step-wise differential abundance Dataset 1](https://resilience-biomarkers-for-aquaculture.github.io/SW-diffabund_stepwise_ds1/), [GMT file for GSEA](https://resilience-biomarkers-for-aquaculture.github.io/SY-augment-diffexp-with-dsea/) |
| Phase 6: Classifier & Validation | [Two-script gene classifier pipeline](https://resilience-biomarkers-for-aquaculture.github.io/SY-gene-classifier-panel/), [Innate gene expression](https://resilience-biomarkers-for-aquaculture.github.io/SY-innate-gene-expression/), [Six-gene biomarker exploration](https://resilience-biomarkers-for-aquaculture.github.io/SY-six-gene-biomarker-exploration/), [Common genes per LOSO fold](https://resilience-biomarkers-for-aquaculture.github.io/SY-plot-DEGs-per-fold/) |

All notebook posts are also indexed in [Sources & References](sources/posts-index.md).

---

**Next:** Browse the complete [Sources & References](sources/posts-index.md) index, or see the [Glossary](glossary.md) for term definitions.
