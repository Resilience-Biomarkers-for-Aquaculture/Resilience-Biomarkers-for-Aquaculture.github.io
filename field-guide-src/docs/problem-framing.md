# Problem Framing

## What is "Resilience" in This Project?

In the context of this research, **resilience** refers to the ability of *Crassostrea virginica* (Eastern oyster) to tolerate or resist stress from *Perkinsus marinus* (Dermo disease) infection and other environmental challenges.

### Defining the Phenotype

Throughout this project, phenotype definition evolved significantly:

!!! warning "Big Lesson #2: Traits Were Oversimplified"
    Early attempts used generalized trait definitions across studies (e.g., "stress" vs. "control"). This oversimplification failed because:
    
    - Different studies measured different stressors (disease, temperature, salinity)
    - Trait effects were weaker than study-specific effects
    - Batch correction couldn't compensate for fundamental differences in experimental design

**Current approach:** Focus on specific, comparable phenotypes:

- **Tolerant/Resistant**: Oysters that survived or showed low infection intensity
- **Sensitive/Susceptible**: Oysters that died or showed high infection intensity
- **Control vs. Treated**: Within-study comparisons before examining resilience

## The Datasets

This project integrates multiple RNA-seq datasets from *C. virginica* exposed to *P. marinus*:

### Primary Datasets

- **Dataset 1**: [Description of dataset 1 characteristics]
- **Dataset 4**: Injected group samples
- **Dataset 5**: [Description of dataset 5 characteristics]
- Additional datasets analyzed for specific questions

### Data Characteristics

**Common challenges:**

- **Batch effects**: Study-specific effects consistently stronger than trait effects
- **Sample size limitations**: Variable n across studies
- **Time point variation**: Sampling at different post-exposure times
- **Technology differences**: Including TAG-seq vs. standard RNA-seq ([issue #26, #28](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/26))

## Key Constraints & Design Decisions

### 1. Integration vs. Post-Integration Analysis

!!! example "Big Lesson #1: When Integrated Analysis Fails"
    **Integrated data analysis does not work** when data is noisy (too much within and/or across study variation) and signal is not strong enough.
    
    **Initial approach (Failed):**
    - Pooled all datasets together
    - Attempted batch correction (COMBAT, RemoveBatchEffect)
    - Expected to see trait-based clustering
    
    **Outcome:** Study-specific effects dominated; trait separation was poor
    
    **Pivot:** Adopted post-data integration approach:
    1. Analyze each dataset independently
    2. Compare results across datasets
    3. Identify reproducible signatures

### 2. Fixed vs. Random Effects

Early attempts didn't properly account for:

- Sample size differences across studies
- Control group considerations
- Study as a random effect in mixed models

### 3. Normalization Timing

!!! info "When Normalization Happens"
    In the `nf-core/differentialabundance` pipeline:
    
    - PCAs are generated **before** differential abundance analysis
    - Normalization (VST, TPM) happens during analysis
    - Batch correction (if applied) occurs on normalized counts
    
    Understanding this order is critical for interpreting preliminary QC plots

### 4. Innate vs. Reactive Biomarkers

!!! success "Big Lesson #3: Don't Filter Out Innate Signals"
    Biomarkers may exist in **control groups** if they represent innate resilience traits (genes that are constitutively different in resistant vs. sensitive individuals).
    
    **Implication:** The stepwise approach (filter controls first) may inadvertently remove true biomarkers
    
    **Resolution:** Developed alternative classifier approach that preserves innate signals (see [Two-Step Classifier](pipelines/classifier-path.md))

See [issue #53](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/issues/53) and [notebook post](https://resilience-biomarkers-for-aquaculture.github.io/SY-innate-gene-expression/) for full innate vs. reactive analysis.

## Biological Context

### *Perkinsus marinus* (Dermo)

- Major pathogen in oyster aquaculture
- Causes dermo disease (mortality and reduced growth)
- Variable infection response across oyster populations
- Understanding resistance is key to breeding programs

### Molecular Signatures of Resilience

This project aims to identify gene expression patterns that:

1. **Predict** tolerance/resistance phenotypes
2. Are **reproducible** across independent studies
3. Are **biologically interpretable** (pathway-informed)
4. Could be **assayable** in breeding programs (minimal gene panels)

## Research Questions Evolution

### Initial Questions (December 2024)

- Can we identify shared stress response genes across multiple studies?
- Do RNA-seq datasets cluster by trait when integrated?

### Refined Questions (January-April 2025)

- How do batch effects and study design differences limit integration?
- Can batch correction methods recover trait-based signal?
- Which normalization approaches work best for meta-analysis?

### Current Questions (August-September 2025)

- What is the minimal gene set that distinguishes tolerant from sensitive oysters?
- Are biomarkers innate or reactive (expressed before vs. after stress)?
- Can classifiers trained on one study predict phenotypes in independent studies?

## Success Criteria

A successful biomarker panel should:

✅ Show reproducible differential expression across datasets  
✅ Maintain predictive accuracy in Leave-One-Study-Out (LOSO) validation  
✅ Be small enough for practical assay development (< 10 genes)  
✅ Have biological interpretability (annotatable, pathway-linked)  
✅ Distinguish phenotypes in both control and treated conditions  

---

**Next:** Explore validated analysis approaches in [Pipelines](pipelines/decision-tree.md)
