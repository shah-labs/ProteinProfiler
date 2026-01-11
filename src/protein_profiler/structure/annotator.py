import subprocess
import os
from Bio.PDB.DSSP import DSSP
from Bio.PDB import PDBParser

class StructureAnnotator:
    def __init__(self, pdb_path, chain_id="A"):
        self.pdb_path = pdb_path
        self.chain_id = chain_id

    def run_dssp(self):
        """Runs DSSP and extracts structural features per residue."""
        print(f"🏗️  Running DSSP on {self.pdb_path}...")
        
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", self.pdb_path)
        model = structure[0]
        
        # Biopython wrapper for DSSP
        # Note: Ensure 'dssp' or 'mkdssp' is in your PATH
        dssp = DSSP(model, self.pdb_path, dssp='dssp') 
        
        results = []
        for key in dssp.keys():
            if key[0] == self.chain_id:
                res_data = dssp[key]
                results.append({
                    'pdb_resnum': key[1][1],
                    'ss8': res_data[2],    # 8-state secondary structure
                    'rsa': res_data[3],    # Relative Solvent Accessibility
                    'phi': res_data[4],
                    'psi': res_data[5]
                })
        
        return results

    @staticmethod
    def map_to_ss3(ss8_char):
        """Maps DSSP 8-state to 3-state (H, E, C)[cite: 57]."""
        if ss8_char in ['H', 'G', 'I']: return 'H' # Helix
        if ss8_char in ['E', 'B']: return 'E'      # Sheet
        return 'C'                                 # Coil/Loop