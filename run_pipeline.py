import os
import pandas as pd
from protein_profiler.io.ingestor import ProteinIngestor
from protein_profiler.msa.searcher import UniProtSearcher
from protein_profiler.msa.aligner import MAFFTAligner
from protein_profiler.msa.qc import MSAQualityControl
from protein_profiler.structure.annotator import StructureAnnotator
from protein_profiler.ptm.annotator import PTMAnnotator
from protein_profiler.io.manifest import RunManifest
from protein_profiler.mutscan.scanner import ESMScanner
from protein_profiler.prioritizer.ranker import SitePrioritizer



# --- Configuration & Initialization ---
FASTA = "data/1gfl.fasta"
PDB = "data/1gfl.pdb"
CHAIN = "A"
OUTPUT_DIR = "results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Initialize the paper trail [cite: 20, 178]
manifest = RunManifest(output_path=f"{OUTPUT_DIR}/run.yaml")
manifest.capture_tool_versions()

print("🚀 Starting ProteinProfiler Pipeline...")

# --- 1. Ingestion [User-Facing Input Contract] ---
# Ensures canonical residue indexing (1..L) [cite: 11-18]
ingestor = ProteinIngestor(fasta_path=FASTA, pdb_path=PDB, chain_id=CHAIN)
df = ingestor.load_data()
manifest.add_parameter("ingestion", {"chain": CHAIN, "fasta": FASTA, "pdb": PDB})

# --- 2. Homology Search & MSA [Stage A] ---
searcher = UniProtSearcher()
raw_msa_path = searcher.search(ingestor.sequence, job_name="1gfl")
manifest.add_parameter("msa_search", {"size": 50, "database": "UniProtKB"})

if raw_msa_path:
    # Align sequences using MAFFT [cite: 28-30]
    aligner = MAFFTAligner()
    aligned_msa_path = aligner.align(raw_msa_path, job_name="1gfl")
    
    # --- 3. MSA Quality Control (Non-negotiable) [Stage A3] ---
    if aligned_msa_path:
        qc = MSAQualityControl(aligned_msa_path)
        summary = qc.get_summary()
        print(f"📊 MSA QC Summary: {summary}")
        
        # Calculate evolutionary signals [cite: 31-36, 44-46]
        df['entropy'] = qc.calculate_entropy()[:len(df)]
        df['gap_fraction'] = qc.calculate_gap_fractions()[:len(df)]
        df['neff_pos'] = summary['neff'] # Simplified per-position Neff
    
    # --- 4. Structure & SASA [Stage C & D] ---
    # Annotates secondary structure and solvent exposure [cite: 50-67]
    annotator = StructureAnnotator(PDB, chain_id=CHAIN)
    dssp_results = annotator.run_dssp()

    for res in dssp_results:
        mask = df['pos'] == res['pdb_resnum']
        if mask.any():
            df.loc[mask, 'ss8'] = res['ss8']
            df.loc[mask, 'ss3'] = annotator.map_to_ss3(res['ss8'])
            df.loc[mask, 'rsa'] = res['rsa']
            df.loc[mask, 'surface_flag'] = res['rsa'] > 0.2
    
    # --- 5. PTM Annotation [Stage E] ---
    # Flags predictive liabilities (e.g., Glycosylation) [cite: 68-84]
    ptm_engine = PTMAnnotator(ingestor.sequence)
    ptm_flags = ptm_engine.get_ptm_flags()

    df['ptm_flags'] = None
    df['ptm_sources'] = None
    for ptm in ptm_flags:
        mask = df['pos'] == ptm['pos']
        if mask.any():
            df.loc[mask, 'ptm_flags'] = ptm['ptm_type']
            df.loc[mask, 'ptm_sources'] = ptm['source']

    # --- Finalize Artifacts ---
    # This CSV is the "Single Source of Truth" [cite: 140-141, 160]
    output_csv = f"{OUTPUT_DIR}/residue_profile.csv"
    df.to_csv(output_csv, index=False)
    
    # Save the manifest to close the audit loop [cite: 176]
    manifest.save()
    
    print(f"\n✅ Pipeline Stage A-E Complete.")
    print(f"📄 Unified Table: {output_csv}")
    print(df[['pos', 'aa_wt', 'entropy', 'ss3', 'surface_flag']].head(10))

else:
    print("❌ Pipeline halted: Homology search failed.")

# 7. Mutation Scanning [Stage G]
scanner = ESMScanner()
mut_data = scanner.scan_sequence(ingestor.sequence)
mut_df = pd.DataFrame(mut_data)

# Record in manifest
manifest.add_parameter("mutation_scan", {"model": "esm2_t6_8M", "scoring": "LLR"})

# Save mutation scan artifact 
mut_df.to_parquet(f"{OUTPUT_DIR}/mut_scan.parquet")

# 8. Add "Mutation Tolerance" to the Unified Table [cite: 152-155]
# Calculate average LLR per position as a proxy for tolerance
tolerance = mut_df.groupby('pos')['esm_llr'].mean().reset_index()
df = df.merge(tolerance.rename(columns={'esm_llr': 'mut_tolerance'}), on='pos')

# 9. Site Prioritization [Stage H]
prioritizer = SitePrioritizer(df)
df = prioritizer.apply_rules()

# Save the final "Single Source of Truth" [cite: 140, 160]
df.to_csv(f"{OUTPUT_DIR}/residue_profile.csv", index=False)

# 10. Generate human-readable report [cite: 139]
top_sites = df[df['priority_class'] == 'safe'].sort_values('mut_tolerance', ascending=False)
print("\n🚀 TOP ENGINEERING CANDIDATES (Safe Diversification):")
print(top_sites[['pos', 'aa_wt', 'entropy', 'mut_tolerance', 'priority_class']].head(10))

print(f"\n✨ PROJECT COMPLETE. Final table and manifest are in {OUTPUT_DIR}/")