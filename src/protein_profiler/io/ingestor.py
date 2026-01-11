import pandas as pd
from Bio import SeqIO, Align
from Bio.PDB import PDBParser, PPBuilder

class ProteinIngestor:
    def __init__(self, fasta_path, pdb_path, chain_id="A"):
        self.fasta_path = fasta_path
        self.pdb_path = pdb_path
        self.chain_id = chain_id
        self.sequence = None
        self.residue_table = None

    def load_data(self):
        # 1. Load FASTA sequence (canonical reference)
        fasta_rec = SeqIO.read(self.fasta_path, "fasta")
        self.sequence = str(fasta_rec.seq)

        # 2. Extract sequence from PDB ATOM records [cite: 14]
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", self.pdb_path)
        
        ppb = PPBuilder()
        pdb_seq = ""
        # Get sequence specifically from the requested chain [cite: 18]
        for pp in ppb.build_peptides(structure[0][self.chain_id]):
            pdb_seq += str(pp.get_sequence())

        # 3. Validate consistency using modern PairwiseAligner [cite: 15]
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        score = aligner.score(self.sequence, pdb_seq)

        # 4. Initialize the Unified Residue Profile Table [cite: 10, 140-143]
        self.residue_table = pd.DataFrame({
            'pos': range(1, len(self.sequence) + 1),
            'aa_wt': list(self.sequence),
            'chain': self.chain_id
        })
        
        print(f"✅ Loaded {len(self.sequence)}aa. PDB Alignment Score: {score}")
        return self.residue_table

    def validate(self):
        if self.sequence is None:
            raise ValueError("Data not loaded. Call load_data() first.")
        return True