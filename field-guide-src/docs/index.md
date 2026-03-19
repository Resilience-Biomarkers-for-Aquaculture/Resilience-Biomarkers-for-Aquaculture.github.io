# Start Here

**Purpose:** This field guide documents considerations for building a biomarker discovery pipeline from RNA-seq meta-analysis, with emphasis on lessons learned, reproducibility, and practical recommendations for researchers working with noisy multi-study omics data.

**Grant mission:** Developing molecular resilience biomarkers to support selective breeding and management of shellfish aquaculture under disease pressure. *(Full narrative: [ProjectSummaryandNarrative.pdf](https://github.com/Resilience-Biomarkers-for-Aquaculture/Resilience-Biomarkers-for-Aquaculture.github.io/blob/master/docs/ProjectSummaryandNarrative.pdf))*

**What was done:** Integrated RNA-seq datasets from *Crassostrea virginica* exposed to *Perkinsus marinus*, conducted differential abundance analyses across multiple independent studies, and developed a two-step classifier pipeline that identifies reproducible gene expression signatures.

**Key results:** Identified a validated 6-gene classifier panel that distinguishes tolerant from sensitive oyster phenotypes, confirmed via Leave-One-Study-Out (LOSO) cross-validation.

**Intended application:** This guide is a reusable template for researchers facing batch effects, weak signals, and overfitting risks in multi-study biomarker discovery.

---

## How to Navigate This Guide

| Section | What You'll Find |
|---|---|
| [1. Research Context & Problem Framing](problem-framing.md) | Background on Dermo disease, datasets, phenotype definitions, and key constraints |
| [2. Process Narrative](process-narrative.md) | Chronological account of decisions, pivots, and surprises — the honest story of what happened |
| [3. Big Lessons Learned](lessons-learned.md) | Distilled, numbered insights for researchers adapting this work |
| [4. Methods & Pipelines](pipelines/decision-guide.md) | Decision guide and validated analysis pipelines |
| [5. Analysis Code & Source Code](code-sources.md) | Where to find and how to run everything |
| [6. Glossary](glossary.md) | Terms defined in the context of this project |
| [7. Sources & References](sources/posts-index.md) | Notebook posts, GitHub issues, and external references |

---

*Developed by Shelly Wanamaker and Steve Yost with AI assistance. Cite as: Wanamaker, S.A. and Yost, S. (2025). Resilience Biomarkers Field Guide. https://resilience-biomarkers-for-aquaculture.github.io/field-guide/*
