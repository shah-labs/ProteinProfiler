import torch
from transformers import EsmTokenizer, EsmForMaskedLM

class ESMScanner:
    def __init__(self, model_name="facebook/esm2_t6_8M_UR50D"):
        print(f"🤖 Loading PLM: {model_name}...")
        self.tokenizer = EsmTokenizer.from_pretrained(model_name)
        self.model = EsmForMaskedLM.from_pretrained(model_name)
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model.to(self.device)
        self.model.eval()

    def scan_sequence(self, sequence):
        """
        Calculates LLR for all possible mutations.
        LLR = log(P(mutant) / P(wildtype)) [cite: 99-100]
        """
        inputs = self.tokenizer(sequence, return_tensors="pt").to(self.device)
        
        with torch.no_grad():
            logits = self.model(**inputs).logits 
        
        # Remove start/stop tokens and get log probabilities
        log_probs = torch.log_softmax(logits, dim=-1)[0, 1:-1, :]
        
        # Define the 20 standard amino acids
        aa_list = list("ACDEFGHIKLMNPQRSTVWY")
        # Ensure tokens is a list of IDs
        tokens = self.tokenizer.convert_tokens_to_ids(aa_list)
        
        results = []
        for i, wt_aa in enumerate(sequence):
            wt_token = self.tokenizer.convert_tokens_to_ids(wt_aa)
            wt_log_prob = log_probs[i, wt_token].item()
            
            # Use zip to iterate through both lists simultaneously
            for mut_aa, mut_token in zip(aa_list, tokens):
                mut_log_prob = log_probs[i, mut_token].item()
                llr = mut_log_prob - wt_log_prob
                
                results.append({
                    'pos': i + 1,
                    'wt': wt_aa,
                    'mut': mut_aa,
                    'esm_llr': llr
                })
        return results
    
    def calculate_consensus(self, mut_df):
        """
        Calculates consensus rank and disagreement across models [cite: 110-112].
        """
        # Convert LLR to ranks (lower LLR = more deleterious = higher rank)
        mut_df['rank_esm'] = mut_df.groupby('pos')['esm_llr'].rank(ascending=True)
        
        # In a real G3 setup, you'd add SaProt or ProteinMPNN ranks here [cite: 104, 106]
        mut_df['consensus_rank'] = mut_df['rank_esm'] 
        mut_df['support_count'] = 1 # Currently 1 model, framework ready for more
        return mut_df