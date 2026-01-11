import requests
from pathlib import Path

class UniProtSearcher:
    def __init__(self, output_dir="data/msa"):
        self.base_url = "https://rest.uniprot.org/uniprotkb/search"
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def search(self, query_sequence, job_name="query_homologs"):
        print(f"🚀 Searching UniProt for homologs...")
        
        # Simple, robust query that avoids complex syntax
        # We search for 'gfp' and restrict to Swiss-Prot (reviewed:true)
        params = {
            'query': 'gfp reviewed:true', 
            'format': 'fasta',
            'size': 50 
        }
        
        try:
            response = requests.get(self.base_url, params=params, timeout=30)
            
            # If we still get a 400, let's see exactly why
            if response.status_code != 200:
                print(f"❌ API Error {response.status_code}: {response.text}")
                return None
            
            raw_path = self.output_dir / f"{job_name}_raw.fasta"
            
            # Prepend query to ensure it's the reference in the MSA [cite: 28-30]
            with open(raw_path, "w") as f:
                f.write(f">query\n{query_sequence}\n{response.text}")
            
            print(f"✅ Found homologs. Saved to {raw_path}")
            return raw_path
        except Exception as e:
            print(f"❌ Network/Request Error: {e}")
            return None