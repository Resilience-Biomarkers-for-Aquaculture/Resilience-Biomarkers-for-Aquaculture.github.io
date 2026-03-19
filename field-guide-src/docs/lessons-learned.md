# Big Lessons Learned

*Distilled, numbered insights for researchers adapting this work to new species or stressors.*

These lessons emerged from 10+ months of iterative analysis. Each one represents a decision point where the team had to abandon a promising approach and rethink the strategy.

---

## Lesson 1: Integrated Analysis Fails with Noisy, Weak-Signal Data

**→ Use post-data integration instead**

When datasets are noisy — too much within- and across-study variation — and the biological signal is not strong enough to overcome that noise, pooling all samples together makes things worse, not better.

**What happened:** Pooling all datasets produced PCAs dominated by study-specific clustering. Trait-based (tolerant vs. sensitive) separation was invisible.

**What to do instead:**
1. Analyze each dataset independently
2. Identify genes that are differentially expressed in multiple datasets
3. Score genes by reproducibility and consistency across studies
4. Build a classifier on the reproducible signal

**See:** [Process Narrative — Phase 1](process-narrative.md#phase-1-dec-2024-jan-2025-integrated-analysis-attempt), [Issue #12](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/12)

---

## Lesson 2: Trait Definitions Must Be Specific and Comparable

**→ Generic "stress vs. control" is insufficient**

Early attempts used broad trait labels (e.g., "stressed" vs. "control") across studies that had different stressors, time points, and experimental designs. The result: study effects dominated trait effects in every analysis.

**What happened:** Batch correction (COMBAT, RemoveBatchEffect) could not recover trait signal because the trait definitions were fundamentally incompatible across studies.

**What to do instead:**
- Define phenotypes precisely (e.g., "tolerant to *P. marinus* injection at Day 7" — not "stressed")
- Compare only studies with harmonized phenotype definitions
- Validate phenotype comparability before combining any data

**See:** [Process Narrative — Phase 2](process-narrative.md#phase-2-feb-2025-confronting-batch-effects), [Research Context](problem-framing.md#defining-the-phenotype)

---

## Lesson 3: Innate Biomarkers Live in Control Groups

**→ Stepwise filtering removes them by design**

The intuitive hypothesis is that resilience biomarkers are *reactive* — genes that change in response to stress differently in tolerant vs. sensitive individuals. But the data showed the opposite: the strongest signals were *innate*, present even in unstressed controls.

**What happened:** The stepwise approach (Step 1: filter for stress-responsive genes; Step 2: compare tolerant vs. sensitive from Step 1 genes) systematically discarded genes that were constitutively different between phenotypes. Only 1 significant gene survived this filter in Dataset 1.

**What to do instead:**
- Test both hypotheses: compare tolerant vs. sensitive in *both* control and treated samples
- If the signal is similar in both conditions → innate biomarker
- Use the two-step classifier approach, which does not filter based on stress response

**See:** [Process Narrative — Phase 5](process-narrative.md#phase-5-aug-2025-stepwise-differential-abundance-discovery-of-innate-signal-problem), [Issue #53](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/53), [Notebook: Innate gene expression](https://resilience-biomarkers-for-aquaculture.github.io/SY-innate-gene-expression/)

---

## Lesson 4: Training-Set Leakage Inflates Accuracy

**→ Prevent leakage at every step: feature selection, normalization, and tuning**

Data leakage — where information from the test set influences the training process — is the most common source of falsely optimistic accuracy in biomarker discovery.

**What happened:** Early evaluations included training samples in the test set during cross-study validation. This inflated apparent accuracy.

**Common leakage points:**
1. **Normalization:** If you normalize train + test together, test data influenced the normalization parameters → leakage
2. **Feature selection:** If you select genes using all data before splitting → leakage
3. **Hyperparameter tuning:** If you tune regularization using the test set → leakage

**What to do instead:**
- Use Leave-One-Study-Out (LOSO) validation: train on all studies except one, test on the held-out study
- Fit normalization parameters on training data only; apply to test
- Select features using training data only
- Never include the training study in the test set for cross-study evaluation

**See:** [Validation & Pitfalls](pipelines/validation.md#1-data-leakage), [Issue #44](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/44)

---

## Lesson 5: Pipeline Internals Matter

**→ PCAs in `nf-core/differentialabundance` are generated *before* normalization**

A subtle but important technical lesson: when using `nf-core/differentialabundance`, the QC PCA plots are generated from raw (or lightly filtered) counts *before* DESeq2 normalization occurs. This means early PCAs do not reflect the data that will actually be used in differential abundance analysis.

**What happened:** Initial PCAs showed poor trait separation, leading to early concern about data quality. Upon understanding the pipeline order, it became clear that these PCAs were not representing the normalized data.

**What to do:**
- Do not draw conclusions about biological signal from pre-normalization PCAs
- Interpret post-analysis PCAs (generated after normalization) for biology
- If you want normalized PCAs, explicitly request them after running DESeq2

**See:** [Research Context — Normalization Timing](problem-framing.md#3-normalization-timing), [Issue #34](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/34)

---

## Lesson 6: Technology Differences Require Different Parameters

**→ TAG-seq ≠ standard RNA-seq**

Not all RNA-seq data is the same. When one study used TAG-seq (3′ tag sequencing) and was processed with standard RNA-seq parameters, the results contained GC bias artifacts that undermined downstream analysis.

**What happened:** The Johnson dataset (Dataset 5) used TAG-seq technology. Standard FastP parameters for full-length RNA-seq were inappropriate: wrong adapter trimming, incorrect quality filters, and GC bias.

**What to do:**
- Identify the library type for every dataset before processing
- Use technology-appropriate parameters (especially for adapter trimming and quality filtering)
- Run QC checks (MultiQC, FastQC) and look for GC bias before proceeding to differential abundance

**See:** [Process Narrative — Phase 3](process-narrative.md#phase-3-apr-jun-2025-shift-to-per-dataset-independent-analysis-tag-seq-discovery), [Issue #26](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/26), [Notebook: TAG-seq FastP params](https://resilience-biomarkers-for-aquaculture.github.io/SW-FastPparams4tagseq/)

---

## Summary Table

| Lesson | Core Failure | Better Approach |
|---|---|---|
| 1. Integrated analysis | Study effects dominate trait signal | Post-data integration: analyze separately, compare |
| 2. Trait definitions | Incomparable phenotype labels across studies | Specific, harmonized phenotype definitions |
| 3. Innate biomarkers | Stepwise filter removes constitutive signals | Test tolerant vs. sensitive in both control and treated |
| 4. Leakage | Training data in test set → inflated accuracy | LOSO; fit normalization and feature selection on training only |
| 5. Pipeline internals | Pre-norm PCAs misinterpreted | Understand and use post-normalization QC outputs |
| 6. Library type | Wrong parameters for TAG-seq | Identify and use technology-appropriate processing |

---

**Next:** See [Methods & Pipelines](pipelines/decision-guide.md) for how these lessons shaped the validated analysis approaches.
