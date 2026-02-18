from Bio import SeqIO
import time
import sys

class ProteinAnalyzer:
    STANDARD = set("ACDEFGHIKLMNPQRSTVWY")
    HYDROPHOBIC = set("AVILMFWYPGC")
    HYDROPHILIC_NEUTRAL = set("STNQH")
    HYDROPHILIC_POSITIVE = set("KR")
    HYDROPHILIC_NEGATIVE = set("DE")

    def __init__(self, fasta_path):
        self.fasta_path = fasta_path
        self.stats = {}

    def analyze(self):
        for rec in SeqIO.parse(self.fasta_path, "fasta"):
            seq = str(rec.seq).upper()
            if any(ch not in self.STANDARD for ch in seq):
                continue

            c = {"hydrophobic": 0, "hydrophilic_neutral": 0,
                 "hydrophilic_positive": 0, "hydrophilic_negative": 0}

            for aa in seq:
                if aa in self.HYDROPHOBIC: c["hydrophobic"] += 1
                elif aa in self.HYDROPHILIC_NEUTRAL: c["hydrophilic_neutral"] += 1
                elif aa in self.HYDROPHILIC_POSITIVE: c["hydrophilic_positive"] += 1
                elif aa in self.HYDROPHILIC_NEGATIVE: c["hydrophilic_negative"] += 1

            self.stats[rec.id] = c

    def save(self, out_path):
        with open(out_path, "w", encoding="utf-8") as f:
            for rid, c in self.stats.items():
                f.write(f"{rid}\t{c['hydrophobic']}\t{c['hydrophilic_neutral']}\t"
                        f"{c['hydrophilic_positive']}\t{c['hydrophilic_negative']}\n")


# --- запуск ---
t0 = time.time()
fasta_file = sys.argv[1]
an = ProteinAnalyzer(fasta_file)
an.analyze()
an.save("stats.tsv")
print("Time (s):", time.time() - t0)
