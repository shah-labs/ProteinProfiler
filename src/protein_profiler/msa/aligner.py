import subprocess
from pathlib import Path

class MAFFTAligner:
    def __init__(self, output_dir="data/msa"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def align(self, input_fasta, job_name="1gfl"):
        """
        Runs MAFFT to align sequences.
        Uses --auto to automatically choose the best strategy (L-INS-i, FFT-NS-2, etc.)
        """
        aligned_path = self.output_dir / f"{job_name}_aligned.fasta"
        print(f"🧬 Aligning sequences with MAFFT for {job_name}...")
        
        # Command: mafft --auto input.fasta > output.fasta
        try:
            with open(aligned_path, "w") as out_file:
                subprocess.run(
                    ["mafft", "--auto", str(input_fasta)], 
                    stdout=out_file, 
                    check=True
                )
            print(f"✅ Alignment complete: {aligned_path}")
            return aligned_path
        except subprocess.CalledProcessError as e:
            print(f"❌ MAFFT failed: {e}")
            return None