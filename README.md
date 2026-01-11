# ProteinProfiler: Systematic Residue-Level Profiling System

[cite_start]ProteinProfiler is a modular bioinformatics pipeline designed to generate systematic, residue-level profiles for enzymes to support engineering decisions[cite: 3]. [cite_start]It integrates evolutionary, structural, and machine learning signals into a single "Source of Truth"[cite: 10, 160].

## 🚀 Features
- [cite_start]**Evolutionary Analysis**: Multi-strategy homology search (UniProt API) and MAFFT alignment with $N_{eff}$ and entropy QC [cite: 22, 31-36].
- [cite_start]**Structural Annotation**: Secondary structure (SS3/SS8) and Solvent Accessibility (SASA) via DSSP[cite: 50, 56, 62].
- [cite_start]**PTM Risk Flagging**: Predictive motif scanning for glycosylation liabilities [cite: 68-71, 79].
- [cite_start]**Deep Learning Mutation Scanning**: Zero-shot variant effect prediction using ESM-2 (Log-Likelihood Ratios) [cite: 95-103].
- [cite_start]**Evidence-Grounded Prioritization**: Deterministic ranking of "Safe Diversification" vs "Liability" sites [cite: 115, 125-129].

## 🛠️ Installation & Setup
1. **System Dependencies**:
   - Linux (Pop!_OS/Ubuntu recommended)
   - MAFFT (`sudo apt install mafft`)
   - DSSP (`sudo apt install dssp`)

2. **Python Environment**:
   ```bash
   python -m venv .venv
   source .venv/bin/activate
   pip install -r requirements.txt