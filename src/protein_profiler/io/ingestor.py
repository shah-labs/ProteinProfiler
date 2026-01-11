from Bio.PDB import PDBParser
from Bio import SeqIO
import pandas as pd

class ProteinIngestor:
    def __init__(self, fasta_path, pdb_path, chain_id="A"):
        self.fasta_path = fasta_path
        self.pdb_path = pdb_path
        self.chain_id = chain_id
        self.structure = None  # This is what was missing
        self.sequence = None

    def load_data(self):
        # Load FASTA Sequence
        record = SeqIO.read(self.fasta_path, "fasta")
        self.sequence = str(record.seq)

        # Load PDB Structure using Biopython
        parser = PDBParser(QUIET=True)
        # We store the structure here so RINGAdapter can use it later
        self.structure = parser.get_structure("protein_id", self.pdb_path)
        
        # Create initial DataFrame
        df = pd.DataFrame({
            'pos': list(range(1, len(self.sequence) + 1)),
            'aa_wt': list(self.sequence)
        })
        return df