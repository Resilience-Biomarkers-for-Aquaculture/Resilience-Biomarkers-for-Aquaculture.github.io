# Methods & Pipelines: Decision Guide

*Practical, decision-first guidance for implementing the analyses.*

This guide helps you choose the right analysis approach for your data. It also documents the two validated pipelines developed in this project and the validation strategies used to confirm results.

---

## Which Pipeline Should You Use?

```mermaid
graph TD
    A[Start: Biomarker Discovery Goal] --> B{Multiple independent<br/>datasets available?}
    B -->|No — single study| C[Use within-study validation only;<br/>generalizability will be limited]
    B -->|Yes| D{Do you believe biomarkers<br/>are reactive to stress?}
    D -->|Unsure — test both| E[Run Two-Step Classifier<br/>on all samples]
    D -->|Yes — reactive only| F{Do you have control<br/>AND treated groups?}
    F -->|No| E
    F -->|Yes| G{Will Step 1 yield<br/>> 100 stress genes?}
    G -->|Unsure / No| E
    G -->|Yes| H[Try Stepwise Differential<br/>Abundance first]
    H --> I{Does Step 2 yield<br/>convincing results?}
    I -->|No| E
    I -->|Yes| J[Validate with LOSO]
    E --> J
```

**When in doubt:** Use the [Two-Step Classifier](classifier-path.md). It is the primary validated approach in this project, does not require control groups, and preserves both innate and reactive biomarkers.

**Use [Stepwise Differential Abundance](stepwise.md) only if** you have a strong prior belief that biomarkers are reactive (not innate), and you have sufficient stress-responsive genes (> 100) surviving Step 1.

---

## Pipeline Overview

| | Stepwise Differential Abundance | Two-Step Classifier |
|---|---|---|
| **Status** | ⚠️ Partially validated | ✅ Validated (6-gene panel) |
| **Requires control groups** | Yes | No |
| **Captures innate biomarkers** | ❌ No | ✅ Yes |
| **Handles small gene sets** | ❌ VST breaks down | ✅ Logistic regression works |
| **Best for** | Reactive biomarkers, large gene sets | Any dataset design |
| **Cross-study validation** | LOSO | LOSO |

---

## Pipeline Details

- **[4b. Stepwise Differential Abundance](stepwise.md)** — Two-step filtering: control vs. treated → resistant vs. sensitive
- **[4c. Two-Step Classifier](classifier-path.md)** — Reproducibility scoring → logistic regression; primary validated approach
- **[4d. Validation & Pitfalls](validation.md)** — Avoiding overfitting, LOSO protocol, batch effects

---

**Next:** Start with the [Two-Step Classifier](classifier-path.md) if you're implementing a new analysis, or read [Validation & Pitfalls](validation.md) to understand cross-study validation requirements.
