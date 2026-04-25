# Process Narrative

*A chronological account of decisions, pivots, and surprises — the honest story of what happened.*

This narrative documents the evolution of analysis approaches from December 2024 through September 2025, organized by phase. Each phase highlights the key decision or pivot that defined that period.

---

## Phase 1 (Dec 2024 – Jan 2025): Integrated Analysis Attempt

### Why Integration Seemed Obvious

The natural starting point was to pool all available datasets together and run a combined analysis — this approach maximizes sample size and, in theory, enables detection of weak signals.

### Week of Dec 4, 2024

**Attempted:** Running nf-core pipeline on 4 combined datasets

**Outcome:** ❌ Failed - compute resource limits and too many unaccounted nuances

!!! quote "Key Insight"
    "Too many nuances that can't be accounted for when everything is pooled"

**Related:**

- [Notebook: Differential abundance workflow exploration](https://resilience-biomarkers-for-aquaculture.github.io/2024-12-11-ES-differential_abundance_workflow/)

### January 2025: Differential Abundance Beginnings

**Jan 3 activities:**

- [Issue #3](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/3): Created merged metadata
- [Issue #4](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/4): Differential abundance initial approach
  - Toy example: ✅ Success
  - GTF file issues encountered
  - Applied mutual information to Perkinsus datasets

**Outcome:** Merged counts table obtained, but separation by trait was weak

**Jan 16:**

**Attempted:** Subset data to improve separation

**Outcome:** ❌ Didn't help significantly

**Considerations identified:**

- Fixed vs. random effects
- Sample size limitations
- Control group treatment
- Batch corrections needed

**Issues opened:**

- [Issue #9](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/9): Add study 5
- [Issue #12](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/12): Decide best methods and execute

!!! failure "Big Lesson #1: When Integration Fails"
    **Integrated data analysis does not work** when data is noisy (too much within and/or across study variation) and signal is not strong enough.

**Pivot discussed:** Meta-analysis approach (inspired by BMC paper with Erin Witkop, Dina co-author)

---

## Phase 2 (Feb 2025): Confronting Batch Effects

**Key finding:** Study-specific effects are much stronger than trait effects

**Attempted:**

- [Issue #18](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/18): Batch correction
  1. COMBAT
  2. RemoveBatchEffect

**Outcome:** ❌ Little effect on improving trait-based separation

!!! failure "Big Lesson #2: Traits Were Oversimplified"
    Generalized trait definitions across studies were insufficient. Study-specific effects dominated.

**Alternative approach:** Limit variation within study by subsetting for common time points, compare DEGs against control groups

---

## Phase 3 (Apr – Jun 2025): Shift to Per-Dataset Independent Analysis; TAG-seq Discovery

### Strategy Shift: Independent Dataset Analysis

**New approach:** Running DifferentialAbundance on each dataset independently

### TAG-seq Parameter Problems (April 2025)

- [Issue #26](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/26): What flags to use? Discovered GC bias
- **Critical finding:** Johnson dataset used TAG-seq; initial RNA-seq analysis parameters were inappropriate
- [Issue #28](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/28): Rerun Johnson data with different parameters

**Related:**

- [Notebook: Reprocess TAG-seq with FastP params](https://resilience-biomarkers-for-aquaculture.github.io/SW-FastPparams4tagseq/)

### Deferred Work (April 2025)

- [Issue #29 & #31](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/31): Ran differential abundance on all datasets together but didn't interpret results yet

!!! note "Could return to this"
    Results exist but interpretation was postponed to focus on per-dataset approach

### Per-Dataset Deep Dives (June 2025)

[Issue #32](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/32): Attempted differentialabundance on datasets separately

- **Steve:** Focused on study 5
- **Shelly:** Focused on study 1
- **Goal:** Understand parameters and how to best run differential abundance

---

## Phase 4 (Jul 2025): Normalization and Study-Combination Experiments

### Study Combination Experiment

[Issue #34](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/34): Combine study 4 (injected group) + study 5

**Research question:** Will study 4 injected group cluster with resistance or susceptible group from study 5?

**Learning:** Gained deeper understanding of the differentialabundance pipeline

- PCAs are generated **before** any differential abundance analysis happens
- Normalization timing matters for interpretation

**Batch correction attempt:**

- Compared PCAs with and without batch correction on top 500 most variable genes
- **Result:** Starting to see evidence of innate trait

!!! note "Could return to this"
    Analysis exists showing 567 genes with significant differential abundance (DESeq). Could revisit to see if these genes show greater clustering than the top 500 most variable.

**Related:**

- [Notebook: Differential abundance on C.virg data](https://resilience-biomarkers-for-aquaculture.github.io/SW-diff-abund/)

### Outstanding Questions

[Issue #39](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/39): Compare DE results from papers and create list of DEGs/markers

!!! warning "Remaining to be done"
    Systematic comparison with published literature still pending

---

## Phase 5 (Aug 2025): Stepwise Differential Abundance; Discovery of Innate-Signal Problem

### GSEA Integration

**Activity:** Run differentialabundance with GSEA

- Created GMT file with gene descriptions
- [Issue #45](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/45): Understand GSEA (not completed)

**Related:**

- [Notebook: Creating GMT file for GSEA](https://resilience-biomarkers-for-aquaculture.github.io/SY-augment-diffexp-with-dsea/)

### Stepwise Differential Abundance

[Issue #41](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/41): Two-step approach

**Steps:**

1. **Step 1:** Controls vs. treated (identify stress-responsive genes)
2. **Step 2:** Resistant vs. sensitive (from step 1 genes)

**Implementation:** [`analyses/stepwise_differentialabundance/`](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/tree/main/analyses/stepwise_differentialabundance)

**Shelly's results on dataset 1:**

- Only 1 significant gene found
- **Problem:** DESeq isn't ideal for highly pared-down gene sets (VST couldn't work well with only ~50 genes)

!!! danger "Big Lesson #3: Biomarkers in Controls"
    Are we removing biomarkers that are **innate**? If resilience biomarkers are constitutively expressed (present in controls), the stepwise filtering approach removes them!

---

## Phase 6 (Aug – Sep 2025): Two-Step Classifier Development; 6-Gene Panel Validation

### Classifier Development Begins

[Issue #42](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/42): Validate SR320 classification results

**Question:** Are the ~50 markers convincing about the difference between sensitive vs. resistant?

[Issue #43](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/43): SR320's AI model

[Issue #44](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/44): Combined datasets 1 & 5

**Comparison:** Integrated data analysis vs. post-data integration approaches

**Result:** ✅ 6-gene classifier completed!

**Pipeline:**

- **Step 1:** Rank genes based on:
  - Reproducibility across datasets
  - Consistency of directionality in expression differences
  - Heterogeneity assessment
- **Step 2:** Logistic regression to find minimum gene set for good classification

!!! success "Six-Gene Classifier Panel Identified"
    Strong separation between tolerant and sensitive phenotypes achieved

!!! warning "Big Lesson #4: Training Set in Test Set"
    **Only include training set in test set if exploring within study.** If trying to predict phenotypes in other studies, you should **definitely not** include the training data in the test set (avoid overfitting/leakage).

**Related:**

- [Notebook: Two-script pipeline for gene classifier](https://resilience-biomarkers-for-aquaculture.github.io/SY-gene-classifier-panel/)
- [Notebook: Stepwise approach on dataset 1](https://resilience-biomarkers-for-aquaculture.github.io/SW-diffabund_stepwise_ds1/)

### September 2025: Validation & Characterization

#### Cross-Study Comparison

[Issue #36](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/36): Run differentialabundance independently for each dataset and compare DEGs

**Theme:** Post-data integration → do we see more overlap?

!!! note "Could return to this"
    Good question about post-data integration vs. integrated data analysis, but uncertain about subsetting approach mentioned in the issue

**Status:** Not completed

#### 6-Gene Panel Characterization

[Issue #49](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/49): Plot 6 genes to gain confidence

**Goal:** Visualize how well these 6 genes distinguish phenotypes

**Related:**

- [Notebook: Exploring six-gene biomarker across studies](https://resilience-biomarkers-for-aquaculture.github.io/SY-six-gene-biomarker-exploration/)

[Issue #53](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/53): Innate vs. reactive gene expression

**Critical insight:** Biomarkers may be constitutively different in resistant vs. sensitive oysters (innate), not just reactive to stress

**Related:**

- [Notebook: Series of DESeq2 runs indicates innate DEGs](https://resilience-biomarkers-for-aquaculture.github.io/SY-innate-gene-expression/)

---

## Key Takeaways

See [Big Lessons Learned](lessons-learned.md) for a distilled list of insights from this project.

### Successful Methodologies

✅ Per-dataset differential abundance analysis  
✅ Post-data integration (compare results across datasets)  
✅ Two-step classifier pipeline (reproducibility + logistic regression)  
✅ Leave-One-Study-Out (LOSO) validation  
✅ Innate vs. reactive DEG characterization  

### Open Questions

- Systematic comparison with published DEG lists
- Optimal visualization of 6-gene panel across studies
- Additional dataset integration for broader validation

---

**Next:** Read the distilled [Big Lessons Learned](lessons-learned.md), or jump to [Methods & Pipelines](pipelines/decision-guide.md)
