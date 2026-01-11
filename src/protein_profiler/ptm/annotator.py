import re
import pandas as pd

class PTMAnnotator:
    def __init__(self, sequence):
        self.sequence = sequence

    def predict_n_glycosylation(self):
        """
        Scans for the N-glycosylation motif: N-X-S/T where X is not Proline.
        This is a standard predictive risk flag.
        """
        # Regex for N-{P}-[ST]-{P}
        pattern = re.compile(r'N[^P][ST][^P]')
        matches = []
        
        for match in pattern.finditer(self.sequence):
            matches.append({
                'pos': match.start() + 1, # 1-based indexing
                'ptm_type': 'N-glycosylation',
                'score': 1.0,
                'source': 'predicted (motif)',
                # Citation moved inside the string below:
                'notes': 'Potential glyco-liability' 
            })
        return matches

    def get_ptm_flags(self):
        """Aggregates all PTM evidence."""
        # For now, we use our motif predictor. 
        # In a full system, you would add dbPTM API calls here.
        all_ptms = self.predict_n_glycosylation()
        return all_ptms