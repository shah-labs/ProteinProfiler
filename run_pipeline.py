from protein_profiler.io.ingestor import ProteinIngestor

# Setup
ingestor = ProteinIngestor(fasta_path="data/1gfl.fasta", pdb_path="data/1gfl.pdb")

# Run
df = ingestor.load_data()
ingestor.validate()

# Inspect
print(df.head())