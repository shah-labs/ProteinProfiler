import pandas as pd

class SitePrioritizer:
    def __init__(self, df):
        self.df = df

    def apply_rules(self):
        """
        Applies deterministic engineering motifs to the residue table[cite: 117, 125].
        """
        # Define thresholds for logic [cite: 118-123]
        high_tolerance = self.df['mut_tolerance'] > -0.5  # ESM LLR closer to 0 is better
        is_surface = self.df['surface_flag'] == True
        is_buried = self.df['surface_flag'] == False
        has_ptm = self.df['ptm_flags'].notna()

        # 1. Safe Diversification 
        self.df.loc[high_tolerance & is_surface, 'priority_class'] = 'safe'
        
        # 2. De-risking Candidates (e.g., removing glyco-sites) 
        self.df.loc[has_ptm & is_surface, 'priority_class'] = 'liability-fix'
        
        # 3. Avoidance (the core stability residues) 
        self.df.loc[is_buried & ~high_tolerance, 'priority_class'] = 'avoid'
        
        # Fill remaining as neutral
        self.df['priority_class'] = self.df['priority_class'].fillna('neutral')
        
        return self.df