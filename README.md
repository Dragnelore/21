Что делает программа

Считывает белковые последовательности из FASTA (Biopython SeqIO.parse)

Для каждой записи считает количество аминокислот по 4 классам:

hydrophobic

hydrophilic_neutral

hydrophilic_positive

hydrophilic_negative

Если в последовательности есть символы не из 20 стандартных аминокислот — запись фильтруется (пропускается)

Сохраняет результат в TSV (таблица, разделитель \t)

Измеряет время работы через time

Требования

Python 3.x

Biopython

Установка Biopython
Вариант 1: conda (рекомендуется)
conda install -c conda-forge biopython


Проверка:

python -c "from Bio import SeqIO; print('OK')"

Вариант 2: pip
pip install biopython


Проверка:

python -c "from Bio import SeqIO; print('OK')"

Подготовка файлов

Сохраните ваш код в файл, например:

protein_analyzer.py

Подготовьте FASTA-файл, например:

insulin.fasta

Формат FASTA:

>id1
MALWMRLLPLLALLALWGPDPAAA...
>id2
MALWTRLLPLLALLALWGPDPAAA...

Запуск (CLI через sys.argv)

Программа принимает:

sys.argv[1] — входной FASTA

sys.argv[2] (опционально) — выходной файл (TSV). Если не указан, используется stats.tsv.

Запуск с одним аргументом (только вход)
python protein_analyzer.py insulin.fasta


Результат сохранится в:

stats.tsv

Запуск с входом и выходом
python protein_analyzer.py insulin.fasta results.tsv

Выходной файл (TSV)

Файл содержит по одной строке на каждую FASTA-запись, которая прошла фильтрацию:

Формат:

<record_id>    <hydrophobic>    <hydrophilic_neutral>    <hydrophilic_positive>    <hydrophilic_negative>


Пример:

sp|P01308|INS_HUMAN    20    8    2    4
sp|P01315|INS_PIG      19    9    2    4

Примечания

Если в последовательности встречаются символы типа X, B, Z, U, * и т.п. — запись пропускается (это и есть “фильтрация записи”).

В терминал дополнительно выводится:

время выполнения (сек)
