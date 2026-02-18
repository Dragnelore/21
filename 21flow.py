from Bio import SeqIO
from queue import Queue
from threading import Thread
import time
import sys

class ProteinAnalyzer:
    STANDARD = set("ACDEFGHIKLMNPQRSTVWY")

    HYDROPHOBIC = set("AVILMFWYPGC")
    HYDROPHILIC_NEUTRAL = set("STNQH")
    HYDROPHILIC_POSITIVE = set("KR")
    HYDROPHILIC_NEGATIVE = set("DE")

    def __init__(self, fasta_path: str):
        self.fasta_path = fasta_path
        self.stats = {}  # {rec_id: {class: count}}

        self.q = Queue(maxsize=1000)  # синхронизирующая очередь
        self.SENTINEL = None

    def _producer(self): # если очередь полная то он ждет maxsize
        for rec in SeqIO.parse(self.fasta_path, "fasta"):
            seq = str(rec.seq).upper() #фильтруем записи чтобы не попали маленькие буквы
            self.q.put((rec.id, seq)) # передаем в очередь айди и саму запись
        self.q.put(self.SENTINEL) # передаем None и останавливамся

    def _consumer(self): # если очередь пустая то он ждет
        while True:
            item = self.q.get() # получаем что передал нам продюссер
            if item is self.SENTINEL: #если None то заканчиваем задачу и потом программу
                self.q.task_done()
                break # все 

            rid, seq = item #пара которую нам передал продюесссер распоковываем

            # фильтрация записи если есть не-стандартные символы пропуск
            if any(ch not in self.STANDARD for ch in seq):
                self.q.task_done()
                continue

            c = {"hydrophobic": 0, "hydrophilic_neutral": 0,
                 "hydrophilic_positive": 0, "hydrophilic_negative": 0}

            for aa in seq:
                if aa in self.HYDROPHOBIC:
                    c["hydrophobic"] += 1
                elif aa in self.HYDROPHILIC_NEUTRAL:
                    c["hydrophilic_neutral"] += 1
                elif aa in self.HYDROPHILIC_POSITIVE:
                    c["hydrophilic_positive"] += 1
                elif aa in self.HYDROPHILIC_NEGATIVE:
                    c["hydrophilic_negative"] += 1

            self.stats[rid] = c # заполняем наш список прочтениями который мы совершили
            self.q.task_done()

    def analyze(self):
        t_prod = Thread(target=self._producer) # создаем поток продюсер
        t_cons = Thread(target=self._consumer) # создаем поток консумер

        t0 = time.time()
        t_prod.start() # пускаем потоки
        t_cons.start()

        t_prod.join() # ждем пока потоки доработают
        self.q.join() # дождаться обработки всех элементов
        t_cons.join()

        return time.time() - t0

    def save(self, out_path: str):
        with open(out_path, "w", encoding="utf-8") as f:
            for rid, c in self.stats.items():
                f.write(
                    f"{rid}\t{c['hydrophobic']}\t{c['hydrophilic_neutral']}\t"
                    f"{c['hydrophilic_positive']}\t{c['hydrophilic_negative']}\n"
                )


fasta_file = sys.argv[1]
an = ProteinAnalyzer(fasta_file)
elapsed = an.analyze()
an.save("stats.tsv")
print("Time (s):", elapsed)
print("Records saved:", len(an.stats))
