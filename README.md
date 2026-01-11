# ProteinProfiler: Systematic Residue-Level Profiling System

**ProteinProfiler** is a modular bioinformatics pipeline designed to generate systematic, residue-level profiles for enzymes to support engineering decisions¹. It integrates evolutionary, structural, and machine learning signals into a single "Source of Truth" CSV².

## 🚀 Features

* **Evolutionary Analysis**: Multi-strategy homology search and alignment with $N_{eff}$ and entropy QC³.
* **Structural Annotation**: Secondary structure (SS3/SS8) and Solvent Accessibility (SASA) via DSSP⁴.
* **PTM Risk Flagging**: Predictive motif scanning for glycosylation liabilities⁵.
* **Deep Learning Mutation Scanning**: Zero-shot variant effect prediction using ESM-2 Log-Likelihood Ratios (LLR)⁶.
* **Evidence-Grounded Prioritization**: Deterministic ranking of "Safe Diversification" vs. "Liability" sites⁷.



## 🏗️ Installation & Setup

1. **System Dependencies**:
   * Linux (Ubuntu/Pop!_OS recommended)
   * MAFFT (`sudo apt install mafft`)
   * DSSP (`sudo apt install dssp`)

2. **Python Environment**:
   ```bash
   python -m venv .venv
   source .venv/bin/activate
   pip install -r requirements.txt