# Proposed Field Guide Outline

**Audience:** Researchers with similar omics datasets (potentially different species or stressors) who want to learn from this project's lessons, methods, and pipelines for biomarker discovery.

---

## 0. Start Here *(Landing Page — revised to be succinct)*

**Format:** Abstract-style, ~250–300 words total.

**Content:**

1. **Purpose** — This field guide documents considerations for building a biomarker discovery pipeline from RNA-seq meta-analysis, with an emphasis on lessons learned, reproducibility, and practical recommendations for researchers working with noisy multi-study omics data.

2. **Grant mission** — One to two sentences drawn from the project narrative: developing molecular resilience biomarkers to support selective breeding and management of shellfish aquaculture under disease pressure. *(Full narrative: [ProjectSummaryandNarrative.pdf](https://github.com/Resilience-Biomarkers-for-Aquaculture/Resilience-Biomarkers-for-Aquaculture.github.io/blob/master/docs/ProjectSummaryandNarrative.pdf))*

3. **What was done** — Two to three sentences: integrated RNA-seq datasets from *Crassostrea virginica* exposed to *Perkinsus marinus*, conducted differential abundance analyses, and developed a two-step classifier pipeline.

4. **Key results** — One to two sentences: identified a validated 6-gene classifier panel that distinguishes tolerant from sensitive oyster phenotypes, confirmed via Leave-One-Study-Out (LOSO) cross-validation.

5. **Intended application** — One sentence: this guide is a reusable template for researchers facing batch effects, weak signals, and overfitting risks in multi-study biomarker discovery.

6. **How to navigate this guide** — A short table or bullet list pointing to the main sections below.

---

## 1. Research Context & Problem Framing

*Why this study, why these data, and what "resilience" means here.*

- Background: Dermo disease (*P. marinus*) and its impact on oyster aquaculture
- Biological question: Can gene expression predict tolerance/resistance phenotypes?
- Dataset overview: study descriptions, species, stressors, sample sizes
- Phenotype definition (tolerant vs. sensitive) and how it evolved

---

## 2. Process Narrative *(replaces "Year in Review")*

*A chronological account of decisions, pivots, and surprises — the honest story of what happened.*

- **Phase 1 (Dec 2024 – Jan 2025):** Integrated analysis attempt; why it failed
- **Phase 2 (Feb 2025):** Confronting batch effects; COMBAT and RemoveBatchEffect trials
- **Phase 3 (Apr – Jun 2025):** Shift to per-dataset independent analysis; TAG-seq parameter discovery
- **Phase 4 (Jul 2025):** Normalization and study-combination experiments
- **Phase 5 (Aug 2025):** Stepwise differential abundance; discovery of innate-signal problem
- **Phase 6 (Aug – Sep 2025):** Two-step classifier development; 6-gene panel validation

---

## 3. Big Lessons Learned

*Distilled, numbered insights for researchers adapting this work to new species or stressors.*

1. **Integrated analysis fails with noisy, weak-signal data** → use post-data integration instead
2. **Trait definitions must be specific and comparable** → generic "stress vs. control" is insufficient
3. **Innate biomarkers live in control groups** → stepwise filtering removes them by design
4. **Training-set leakage inflates accuracy** → prevent leakage at every step (feature selection, normalization, tuning)
5. **Pipeline internals matter** → PCAs in `nf-core/differentialabundance` are generated *before* normalization
6. **Technology differences require different parameters** → TAG-seq ≠ standard RNA-seq

---

## 4. Methods & Pipelines

*Practical, decision-first guidance for implementing the analyses.*

### 4a. Pipeline Decision Guide
- Flowchart: Which approach to choose and when *(box diagram/flowchart graphic)*

### 4b. Stepwise Differential Abundance Pipeline
- Step 1: control vs. treated (stress-responsive genes)
- Step 2: resistant vs. sensitive (from Step 1 genes)
- Normalization approach; known failure modes

### 4c. Two-Step Classifier Pipeline *(primary validated approach)*
- Step 1: Reproducibility and directionality scoring across datasets *(box diagram/flowchart graphic)*
- Step 2: Logistic regression to minimize feature set
- How the 6-gene panel was identified

### 4d. Validation & Pitfalls
- Avoiding overfitting; train/test hygiene
- Leave-One-Study-Out (LOSO) protocol
- Detecting and handling batch effects
- Cross-study generalizability considerations

---

## 5. Analysis Code & Source Code

*Entry point for reproducibility — where to find and how to run everything.*

- GitHub repository structure overview
- Key analysis scripts with brief descriptions
- How to run pipelines from scratch (environment setup, inputs, outputs)
- Links to relevant notebook entries (organized by pipeline/phase, not just chronologically)

---

## 6. Glossary

*Terms defined in the context of this project.*

- Domain terms (e.g., Dermo, tolerance phenotype)
- Statistical/computational methods (DESeq2, VST, LOSO, logistic regression)
- Acronyms and abbreviations

---

## 7. Sources & References

*Primary source material underlying the guide.*

- Notebook posts index (analysis logs)
- GitHub issues index (decision trail)
- External publications cited

---

## Navigation Changes Summary

| Current Section | Proposed Section | Change |
|---|---|---|
| Start Here | Start Here | Rewritten as concise abstract |
| Problem Framing | Research Context & Problem Framing | Expanded with dataset details |
| Year in Review | Process Narrative | Renamed; framed as narrative |
| *(none)* | Big Lessons Learned | **New section** |
| Pipelines (3 sub-pages) | Methods & Pipelines (4 sub-pages) | Restructured; decision guide added |
| *(none)* | Analysis Code & Source Code | **New section** |
| Glossary | Glossary | Unchanged |
| Sources (posts + issues) | Sources & References | Moved to end; scoped as references |

---

*This outline is a proposal for review. Full implementation will proceed after approval.*
