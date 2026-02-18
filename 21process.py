from Bio import SeqIO
import multiprocessing as mp
import time
import sys

STANDARD = set("ACDEFGHIKLMNPQRSTVWY")
HYDROPHOBIC = set("AVILMFWYPGC")
HYDROPHILIC_NEUTRAL = set("STNQH")
HYDROPHILIC_POSITIVE = set("KR")
HYDROPHILIC_NEGATIVE = set("DE")


def worker(task):
    rid, seq = task
    seq = seq.upper()

    # фильтрация записи целиком
    if any(ch not in STANDARD for ch in seq):
        return None

    c = {"hydrophobic": 0, "hydrophilic_neutral": 0,
         "hydrophilic_positive": 0, "hydrophilic_negative": 0}

    for aa in seq:
        if aa in HYDROPHOBIC:
            c["hydrophobic"] += 1
        elif aa in HYDROPHILIC_NEUTRAL:
            c["hydrophilic_neutral"] += 1
        elif aa in HYDROPHILIC_POSITIVE:
            c["hydrophilic_positive"] += 1
        elif aa in HYDROPHILIC_NEGATIVE:
            c["hydrophilic_negative"] += 1

    return rid, c


class ProteinAnalyzer:
    def __init__(self, fasta_path: str, processes: int | None = None):
        self.fasta_path = fasta_path
        self.processes = processes or mp.cpu_count()
        self.stats = {}

    def _tasks(self):
        for rec in SeqIO.parse(self.fasta_path, "fasta"):
            yield rec.id, str(rec.seq)  # только (id, строка), чтобы легко сериализовалось

    def analyze(self):
        t0 = time.time()
        with mp.Pool(processes=self.processes) as pool:
            for res in pool.imap_unordered(worker, self._tasks(), chunksize=100):
                if res is None:
                    continue
                rid, c = res
                self.stats[rid] = c
        return time.time() - t0

    def save(self, out_path: str):
        with open(out_path, "w", encoding="utf-8") as f:
            for rid, c in self.stats.items():
                f.write(
                    f"{rid}\t{c['hydrophobic']}\t{c['hydrophilic_neutral']}\t"
                    f"{c['hydrophilic_positive']}\t{c['hydrophilic_negative']}\n"
                )


if __name__ == "__main__":
    fasta_file = sys.argv[1]
    an = ProteinAnalyzer(fasta_file, processes=None)  # None -> все ядра
    elapsed = an.analyze()
    an.save("stats.tsv")
    print("Time (s):", elapsed)
    print("Records saved:", len(an.stats))
