# Validation & Pitfalls

## Overview

Validation is critical for ensuring biomarker panels generalize beyond the discovery cohort. This guide covers validation strategies and common pitfalls encountered in this project, with emphasis on what actually went wrong and how to avoid it.

---

## Validation Strategies

### 1. Leave-One-Study-Out (LOSO) Cross-Validation

**What it is:** Train classifier on all datasets except one, test on the held-out dataset

**Why it's essential:**

- Tests **cross-study generalization**
- Avoids overfitting to study-specific effects
- Simulates real-world scenario: predict phenotypes in new, unseen studies

**Implementation:**

```python
datasets = [dataset1, dataset2, dataset3, dataset4, dataset5]

for test_idx, test_dataset in enumerate(datasets):
    # Training set: all except test_dataset
    train_datasets = [ds for i, ds in enumerate(datasets) if i != test_idx]
    
    X_train = combine_expression_data(train_datasets, gene_panel)
    y_train = combine_labels(train_datasets)
    
    # Fit classifier
    classifier.fit(X_train, y_train)
    
    # Test on held-out study
    X_test = test_dataset.expression[gene_panel]
    y_test = test_dataset.labels
    
    accuracy = classifier.score(X_test, y_test)
    print(f"LOSO Fold {test_idx+1} (test={test_dataset.name}): {accuracy:.3f}")
```

**Interpretation:**

- ✅ **Good:** Consistent accuracy (e.g., 0.80-0.90) across all folds
- ⚠️ **Warning:** One fold significantly lower → investigate study-specific effects
- ❌ **Poor:** High variance in accuracy or < 0.70 average → panel doesn't generalize

!!! warning "Big Lesson #4: Training Set in Test Set"
    **Never include the training study in the test set** when evaluating cross-study performance. This leads to data leakage and falsely optimistic results.
    
    ✅ Correct: LOSO (train on studies 1,2,3; test on study 4)  
    ❌ Wrong: Train on studies 1,2,3,4; test on studies 1,2,3,4

**Related:** [Issue #44](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/44)

### 2. Within-Study Cross-Validation

**What it is:** Hold out samples within a single study (k-fold or train/test split)

**When to use:**

- Initial model development and parameter tuning
- Assessing overfitting within a study
- When you only have one study available (but be cautious about generalizability)

**Limitations:**

- Does NOT test cross-study generalization
- May overestimate performance if study-specific effects are strong
- Use as a first step, but follow with LOSO

### 3. Independent Validation Cohort

**Gold standard:** Test on completely independent dataset never used during development

**Challenges:**

- Requires finding additional suitable datasets
- May have different experimental designs or phenotype definitions
- Labor-intensive to process and integrate

**Status in this project:** [Issue #54](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/54) (find more datasets) - postponed

---

## Common Pitfalls

### 1. Data Leakage

**Definition:** Information from test set "leaks" into training process, causing overly optimistic performance estimates

#### Leakage Scenario 1: Normalization Across Training+Test

❌ **Wrong:**

```python
# Normalize all data together
all_data_normalized = vst_transform(combine(train_data, test_data))

# Then split
X_train = all_data_normalized[train_indices]
X_test = all_data_normalized[test_indices]  # LEAKAGE!
```

✅ **Correct:**

```python
# Fit normalization on training data only
vst_params = fit_vst(train_data)

# Apply to train and test separately
X_train = apply_vst(train_data, vst_params)
X_test = apply_vst(test_data, vst_params)
```

#### Leakage Scenario 2: Feature Selection on All Data

❌ **Wrong:**

```python
# Select features using all data
selected_genes = select_top_degs(all_data, phenotypes)  # LEAKAGE!

# Then split for training/testing
X_train = train_data[selected_genes]
X_test = test_data[selected_genes]
```

✅ **Correct:**

```python
# Select features using training data only
selected_genes = select_top_degs(train_data, train_phenotypes)

# Apply to test
X_train = train_data[selected_genes]
X_test = test_data[selected_genes]
```

**In this project:** LOSO properly avoids leakage by keeping test studies completely separate during training

### 2. Batch Effects

**Definition:** Non-biological variation between studies due to technical differences (sequencing platform, library prep, time, lab)

!!! failure "Big Lessons #1 & #2: When Batch Effects Dominate"
    - Study-specific effects were consistently **stronger than trait effects**
    - Attempted corrections (COMBAT, RemoveBatchEffect) had **little effect** on trait separation
    - **Pivot:** Abandoned integrated analysis in favor of post-data integration

**Observations:**

- PCA plots showed clustering by **study**, not by phenotype
- Trait-based separation was weak even after batch correction
- Integration worked only when combining very similar studies (e.g., study 4 injected + study 5)

#### Batch Correction Attempts

**Tried:**

- [Issue #18](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/18): COMBAT and RemoveBatchEffect

**Outcome:** ❌ Insufficient improvement

**Lesson:** Batch correction is not a panacea. When study effects are too strong, no correction method will recover trait signal.

#### Post-Data Integration Solution

Instead of correcting batches and pooling:

1. Analyze each study independently
2. Identify genes significant in multiple studies
3. Build classifier on reproducible genes

**Advantage:** Reproducibility across studies implicitly handles batch effects

**Related:** [Timeline - February 2025](../timeline.md#february-2025-confronting-batch-effects)

### 3. Normalization Timing & Method

#### Understanding the Pipeline

In `nf-core/differentialabundance`:

1. **Raw counts** are input
2. **PCAs are generated** before differential abundance analysis (on raw or lightly filtered counts)
3. **Normalization** (VST, rlog, TPM) happens during DESeq2 analysis
4. **Batch correction** (if enabled) is applied post-normalization

!!! info "Normalization Insight from July 2025"
    Understanding **when** normalization happens is critical for interpreting QC plots:
    
    - Initial PCAs show raw/lightly-processed data
    - Don't expect trait separation in initial PCAs if normalization hasn't been applied
    - Post-analysis PCAs on normalized data are more informative

**Related:** [Issue #34](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/34), [Timeline - July 2025](../timeline.md#july-2025-combining-studies--understanding-normalization)

#### Normalization Methods

**VST (Variance-Stabilizing Transformation):**

- Recommended for visualization and clustering
- Requires sufficient genes (> 1000) for stable estimation
- Breaks down with small gene sets (< 100)

!!! danger "VST with Small Gene Sets"
    In the stepwise approach, Step 1 filtering left only ~50 genes. DESeq2's VST could not work properly with such a small set.
    
    **Symptom:** Unreliable dispersion estimates, poor model fit  
    **Solution:** Use alternative normalization (TMM, TPM) or abandon stepwise approach

**Related:** [Stepwise Pipeline](decision-tree.md#1-deseq2-breaks-down-with-small-gene-sets)

### 4. Oversimplified Phenotype Definitions

!!! failure "Big Lesson #2: Traits Were Oversimplified"
    Initial approach used generalized trait labels across studies:
    
    - "Stress" vs. "Control"
    - "Resistant" vs. "Sensitive"
    
    **Problem:** Different studies measured different stressors with different designs:
    
    - Study 1: Disease challenge, time point A
    - Study 2: Temperature stress, time point B
    - Study 3: Disease challenge, time point C, different infection route
    
    **Result:** Trait effects were **much weaker** than study-specific effects

**Solution:** Use more specific, harmonized phenotype definitions

- Compare only studies with similar experimental designs
- Within-study trait comparisons first, then assess reproducibility across studies
- Document and respect phenotype heterogeneity

**Related:** [Problem Framing - Defining the Phenotype](../problem-framing.md#defining-the-phenotype)

### 5. Small Sample Sizes

**Challenges:**

- Insufficient power to detect true DEGs
- High variance in estimates
- Overfitting risk (model learns noise, not signal)

**Mitigations:**

- Use regularized models (LASSO, elastic net) to reduce overfitting
- Focus on genes reproducible across studies (increases effective n)
- Report confidence intervals, not just point estimates

**In this project:** Combined datasets (e.g., [Issue #34](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/34) merged study 4 + 5) to increase sample size

### 6. Innate vs. Reactive Biomarkers

!!! success "Big Lesson #3: Biomarkers in Controls"
    Biomarkers may be **constitutively expressed** (innate) in resistant vs. sensitive individuals, even in the absence of stress.
    
    **Implication:** If you filter based on "stress-responsive" (control vs. treated), you will remove innate biomarkers.

**Detection:**

- Compare resistant vs. sensitive in **control samples only**
- If significant → innate biomarker
- If not significant in controls but significant in treated → reactive biomarker

**Pipeline choice:**

- ✅ **Two-step classifier** preserves innate biomarkers
- ❌ **Stepwise approach** removes innate biomarkers by design

**Related:**

- [Issue #53: Innate vs. reactive analysis](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/53)
- [Notebook: Innate gene expression](https://resilience-biomarkers-for-aquaculture.github.io/SY-innate-gene-expression/)

---

## Validation Checklist

Before claiming a validated biomarker panel:

- [ ] **Cross-study validation:** LOSO accuracy > 70% (ideally > 80%)
- [ ] **No data leakage:** Normalization, feature selection, and hyperparameter tuning done independently per fold
- [ ] **Batch effect assessment:** Confirmed genes are reproducible across studies, not driven by batch
- [ ] **Sample size adequacy:** Each phenotype group has n > 10 per study (if possible)
- [ ] **Phenotype definition:** Clear, specific, and documented; comparable across studies used
- [ ] **Innate/reactive characterization:** Assessed whether biomarkers are pre-existing or induced
- [ ] **Biological validation:** Genes have known or plausible functional roles in stress response

---

## Troubleshooting Guide

### Poor LOSO Performance

**Symptom:** Accuracy < 0.70 or high variance across folds

**Possible causes:**

1. **Study-specific effects too strong** → Filter datasets or use more stringent reproducibility criteria in Step 1
2. **Phenotype heterogeneity** → Ensure comparable phenotype definitions
3. **Small sample sizes** → Combine compatible studies or seek additional data
4. **Overfitting** → Increase regularization penalty or reduce feature set

### Genes Don't Validate

**Symptom:** DEGs from discovery cohort not significant in validation cohort

**Possible causes:**

1. **Batch effects** → Apply batch-aware normalization or use post-data integration
2. **False discoveries** → Discovery set had high FDR; use more stringent thresholds
3. **Different experimental conditions** → Validation study not truly comparable
4. **Underpowered validation** → Validation cohort too small to detect effects

### Suspect Data Leakage

**Symptom:** Training accuracy much lower than test accuracy, or perfect test accuracy

**Diagnostic steps:**

1. Re-run analysis with explicit train/test separation at each step
2. Check if normalization used all data before splitting
3. Verify feature selection used only training data
4. Ensure no test samples used in hyperparameter tuning

---

## Related Resources

**Key issues:**

- [#18: Batch correction attempts](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/18)
- [#34: Combining studies](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/34)
- [#44: Classifier development with LOSO](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/44)
- [#53: Innate vs. reactive](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/53)

**Key notebooks:**

- [Gene classifier panel (LOSO validation)](https://resilience-biomarkers-for-aquaculture.github.io/SY-gene-classifier-panel/)
- [Innate gene expression analysis](https://resilience-biomarkers-for-aquaculture.github.io/SY-innate-gene-expression/)
- [Six-gene biomarker exploration](https://resilience-biomarkers-for-aquaculture.github.io/SY-six-gene-biomarker-exploration/)

---

**Next:** Explore detailed project history in the [Timeline](../timeline.md) or browse [Source Indices](../sources/issues-index.md)
