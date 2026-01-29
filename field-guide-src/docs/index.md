# Start Here

## Welcome to the Resilience Biomarkers Field Guide

This field guide documents the real journey of discovering and validating biomarkers for stress resilience in oyster aquaculture using RNA-seq meta-analysis. Rather than presenting an idealized workflow, this guide shares what actually worked, what didn't, and why—grounded in a year-long research effort analyzing multiple *Crassostrea virginica* (Eastern oyster) and *Perkinsus marinus* (parasite) RNA-seq datasets.

## What This Guide Is

!!! info "A Practical Chronicle"
    This is **not** a theoretical textbook. It's a field guide built from:
    
    - **Real analysis decisions** documented in [50+ notebook entries](sources/posts-index.md)
    - **Actual pivots and failures** tracked across [25+ GitHub issues](sources/issues-index.md)
    - **Hard-won lessons** about batch effects, leakage, and when methods break down
    - **A working 6-gene classifier** that emerged from systematic iteration

## What You'll Find Here

### [Problem Framing](problem-framing.md)
Understand what "resilience" means in this project, the datasets we're working with, and the biological and technical constraints that shaped our approach.

### [Pipelines](pipelines/decision-tree.md)
Two validated analysis paths:

1. **Stepwise Differential Abundance**: Control vs. treated → resistant vs. sensitive
2. **Two-Step Classifier Pipeline**: Reproducibility scoring → logistic regression

Both with explicit guidance on normalization, validation (LOSO), and common pitfalls (batch effects, leakage).

### [Year in Review](timeline.md)
A chronological walk through the project from December 2024 to September 2025, highlighting:

- Major methodological pivots (integrated → post-integrated analysis)
- "Big Lessons" about data integration and trait definition
- When to abandon an approach vs. when to iterate

### [Sources](sources/issues-index.md)
Complete index of the [GitHub issues](sources/issues-index.md) and [notebook posts](sources/posts-index.md) referenced throughout this guide, organized by theme for easy exploration.

### [Glossary](glossary.md)
Key terms and concepts defined in the context of this specific project.

## How to Use This Guide

**If you're new to biomarker discovery:**  
Start with [Problem Framing](problem-framing.md) to understand the research context, then explore the [Timeline](timeline.md) to see how approaches evolved.

**If you're implementing a similar analysis:**  
Jump to [Pipelines](pipelines/decision-tree.md) for validated workflows, then check [Validation & Pitfalls](pipelines/validation.md) for things that will bite you.

**If you're troubleshooting an issue:**  
Use the [Sources](sources/issues-index.md) indices to find relevant discussions—chances are we hit the same problem.

**If you're evaluating this approach:**  
Read the "Big Lessons" in the [Timeline](timeline.md) to understand when and why methods were abandoned or refined.

## Core Finding: The 6-Gene Classifier

After extensive iteration, this project identified a **6-gene panel** that successfully distinguishes tolerant from sensitive oyster phenotypes:

- Emerged from post-data integration approach (combining datasets 1 & 5)
- Validated using Leave-One-Study-Out (LOSO) cross-validation
- Uses logistic regression to minimize feature set while maintaining accuracy
- See [Two-Step Classifier](pipelines/classifier-path.md) for full methodology

## Related Resources

- **[Main Project Website](../)** - Jekyll site with publications, presentations, and updates
- **[Analysis Notebook](../notebook)** - All notebook posts organized by date
- **[GitHub Repository](https://github.com/Resilience-Biomarkers-for-Aquaculture)** - Code, data, and issue tracking
- **[Cvirg_Pmarinus_RNAseq Repo](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq)** - Primary analysis repository with issues

## Attribution & Citation

This field guide is developed by Shelly A. Wanamaker as part of the Resilience Biomarkers for Aquaculture project.

If you use methods or insights from this guide, please cite:

> Wanamaker, S.A. (2025). Resilience Biomarkers Field Guide. Resilience Biomarkers for Aquaculture Project. https://resilience-biomarkers-for-aquaculture.github.io/field-guide/

---

*This guide is actively maintained and reflects work through September 2025. Check the [Timeline](timeline.md) for the most recent updates.*
