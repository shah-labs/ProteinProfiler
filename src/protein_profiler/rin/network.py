import numpy as np
from Bio.PDB import NeighborSearch

class RINGAdapter:
    def __init__(self, pdb_path, chain_id="A"):
        self.pdb_path = pdb_path
        self.chain_id = chain_id

    def extract_features(self, structure):
        """Extracts degree (neighbor count) as a proxy for RING centrality."""
        atoms = [atom for atom in structure.get_atoms() if atom.get_parent().get_parent().id == self.chain_id]
        ns = NeighborSearch(atoms)
        
        # RING requirement: degree and weighted degree [cite: 88]
        features = {}
        for residue in structure.get_residues():
            if residue.get_parent().id == self.chain_id:
                # Find all atoms within 5.0 Angstroms of this residue's CA
                if 'CA' in residue:
                    neighbors = ns.search(residue['CA'].coord, 5.0, level='R')
                    features[residue.id[1]] = {
                        'degree': len(neighbors),
                        'w_degree': len(neighbors) * 1.0 # Simplified weight
                    }
        return features