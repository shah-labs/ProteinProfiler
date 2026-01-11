import numpy as np
from Bio import AlignIO
from collections import Counter
import math

class MSAQualityControl:
    def __init__(self, aligned_fasta):
        self.alignment = AlignIO.read(aligned_fasta, "fasta")
        self.num_sequences = len(self.alignment)
        self.alignment_len = self.alignment.get_alignment_length()

    def calculate_gap_fractions(self):
        """Calculates percentage of gaps per position[cite: 33]."""
        gap_fractions = []
        for i in range(self.alignment_len):
            column = self.alignment[:, i]
            gap_count = column.count("-")
            gap_fractions.append(gap_count / self.num_sequences)
        return np.array(gap_fractions)

    def calculate_entropy(self):
        """Calculates Shannon Entropy (conservation) per position[cite: 36]."""
        entropies = []
        for i in range(self.alignment_len):
            column = self.alignment[:, i].replace("-", "")
            if not column:
                entropies.append(0.0)
                continue
            counts = Counter(column)
            total = len(column)
            entropy = -sum((count/total) * math.log2(count/total) for count in counts.values())
            entropies.append(entropy)
        return np.array(entropies)

    def calculate_neff(self, identity_threshold=0.8):
        """
        Calculates Neff (Effective Number of Sequences)[cite: 32].
        Down-weights sequences with >80% identity to reduce bias[cite: 39].
        """
        # Convert alignment to a matrix for faster identity calculation
        align_array = np.array([list(rec) for rec in self.alignment])
        n_seq = self.num_sequences
        weights = np.ones(n_seq)

        # Simplified weight calculation: 
        # Sequences with many neighbors above threshold get lower weight
        for i in range(n_seq):
            matches = np.sum(align_array == align_array[i], axis=1) / self.alignment_len
            neighbors = np.sum(matches >= identity_threshold)
            weights[i] = 1.0 / neighbors
            
        return np.sum(weights)

    def get_summary(self):
        """Generates the msa_qc.json data[cite: 43]."""
        return {
            "n_sequences": self.num_sequences,
            "neff": self.calculate_neff(),
            "alignment_length": self.alignment_len,
            "avg_gap_fraction": np.mean(self.calculate_gap_fractions())
        }