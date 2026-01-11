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
from protein_profiler.rin.network import RINGAdapter

# --- Configuration ---
FASTA = "data/1gfl.fasta"
PDB = "data/1gfl.pdb"
CHAIN = "A"
OUTPUT_DIR = "results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

manifest = RunManifest(output_path=f"{OUTPUT_DIR}/run.yaml")
manifest.capture_tool_versions()

print("🚀 Starting ProteinProfiler Pipeline...")

# --- 1. Ingestion ---
ingestor = ProteinIngestor(fasta_path=FASTA, pdb_path=PDB, chain_id=CHAIN)
df = ingestor.load_data() # This now populates ingestor.structure

# --- 2 & 3. MSA & QC ---
searcher = UniProtSearcher()
raw_msa_path = searcher.search(ingestor.sequence, job_name="1gfl")
if raw_msa_path:
    aligner = MAFFTAligner()
    aligned_msa_path = aligner.align(raw_msa_path, job_name="1gfl")
    if aligned_msa_path:
        qc = MSAQualityControl(aligned_msa_path)
        df['entropy'] = qc.calculate_entropy()[:len(df)]
        df['gap_fraction'] = qc.calculate_gap_fractions()[:len(df)]

# --- 4. Structure (DSSP) ---
annotator = StructureAnnotator(PDB, chain_id=CHAIN)
dssp_results = annotator.run_dssp()
for res in dssp_results:
    mask = df['pos'] == res['pdb_resnum']
    if mask.any():
        df.loc[mask, 'ss3'] = annotator.map_to_ss3(res['ss8'])
        df.loc[mask, 'surface_flag'] = res['rsa'] > 0.2

# --- 5. PTMs ---
ptm_engine = PTMAnnotator(ingestor.sequence)
for ptm in ptm_engine.get_ptm_flags():
    mask = df['pos'] == ptm['pos']
    df.loc[mask, 'ptm_flags'] = ptm['ptm_type']

# --- 6. Stage F: Interaction Network (RING) ---
print("🕸️  Calculating Residue Interaction Network (RING) features...")
ring_adapter = RINGAdapter(PDB, chain_id=CHAIN)
# The ingestor now has the .structure attribute required here:
ring_features = ring_adapter.extract_features(ingestor.structure)

df['degree'] = 0.0 # Initialize columns
df['w_degree'] = 0.0

for pos, feats in ring_features.items():
    mask = df['pos'] == pos
    if mask.any():
        df.loc[mask, 'degree'] = feats['degree']
        df.loc[mask, 'w_degree'] = feats['w_degree']

# --- 7. Stage G: Mutation Scanning ---
scanner = ESMScanner()
mut_data = scanner.scan_sequence(ingestor.sequence)
mut_df = pd.DataFrame(mut_data)
mut_df = scanner.calculate_consensus(mut_df)
mut_df.to_parquet(f"{OUTPUT_DIR}/mut_scan.parquet")

consensus_summary = mut_df.groupby('pos')['esm_llr'].max().reset_index()
df = df.merge(consensus_summary.rename(columns={'esm_llr': 'top_mut_score'}), on='pos')

# --- 8. Stage H: Prioritization ---
median_degree = df['degree'].median()
df['priority_class'] = 'neutral'
safe_mask = (df['top_mut_score'] > -0.5) & (df['surface_flag'] == True) & (df['degree'] < median_degree)
df.loc[safe_mask, 'priority_class'] = 'safe'

# --- Wrap up ---
df.to_csv(f"{OUTPUT_DIR}/residue_profile.csv", index=False)
manifest.save()

print("\n🚀 TOP ENGINEERING CANDIDATES:")
print(df[df['priority_class'] == 'safe'].sort_values('top_mut_score', ascending=False).head(5))
print(f"\n✨ Done. Results in {OUTPUT_DIR}/")