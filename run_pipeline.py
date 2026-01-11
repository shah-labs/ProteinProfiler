import os
import pandas as pd
from protein_profiler.io.ingestor import ProteinIngestor
from protein_profiler.msa.searcher import UniProtSearcher
from protein_profiler.msa.aligner import MAFFTAligner
from protein_profiler.msa.qc import MSAQualityControl

# Define paths
FASTA = "data/1gfl.fasta"
PDB = "data/1gfl.pdb"

# 1. Ingestion [Stage 1]
ingestor = ProteinIngestor(fasta_path=FASTA, pdb_path=PDB)
df = ingestor.load_data()

# 2. Search [Stage A1]
searcher = UniProtSearcher()
raw_msa_path = searcher.search(ingestor.sequence, job_name="1gfl")

if raw_msa_path:
    # 3. Align [Stage A2]
    aligner = MAFFTAligner()
    aligned_msa_path = aligner.align(raw_msa_path, job_name="1gfl")
    
    # 4. Quality Control [Stage A3]
    if aligned_msa_path:
        qc = MSAQualityControl(aligned_msa_path)
        summary = qc.get_summary()
        print(f"📊 MSA QC Summary: {summary}")
        
        # Calculate evolutionary metrics [cite: 31-36]
        # We must align the MSA results back to the query sequence length
        entropies = qc.calculate_entropy()
        gaps = qc.calculate_gap_fractions()
        
        # Attach to the Unified Residue Profile Table [cite: 140, 160]
        # Note: We take only the first L positions (matching our query)
        df['entropy'] = entropies[:len(df)]
        df['gap_fraction'] = gaps[:len(df)]
        
        print("\n--- Updated Residue Profile Table (First 5 residues) ---")
        print(df.head())
        
        # Save progress
        os.makedirs("results", exist_ok=True)
        df.to_csv("results/residue_profile.csv", index=False)
        print("\n✅ Stage A complete. Profile saved to results/residue_profile.csv")
else:
    print("❌ Pipeline halted: Homology search failed.")