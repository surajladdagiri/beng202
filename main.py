import csv
from code import Score

# Read fluorescent aptamers from CSV
fluorescent_seqs = {}
with open('data/fluorescent_aptamers.csv', mode='r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        fluorescent_seqs[row['name']] = row['sequence']

mangoI = fluorescent_seqs.get('mangoI')

# Read genes from CSV
genes = {}
with open('data/genes.csv', mode='r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        genes[row['name']] = row['sequence']

if not fluorescent_seqs:
    raise RuntimeError("Please pass in some fluorescent sequences!")

if not genes:
    raise RuntimeError("Please pass in at least one gene sequence!")

for g_name, g_seq in genes.items():
    print(f"Running for gene {g_name}")
    best_score = float('-inf')
    best_f = ""
    for i in range(len(fluorescent_seqs)):
        curr_f = fluorescent_seqs[i]
        sub_str, score = Score.score(g_seq, curr_f)
        if best_score < score:
            best_score, best_f = score, curr_f
        print(f"For fluorophore sequence {i + 1}:")
        print(f"Optimal substring in gene: {sub_str}\tScore: {score}")

    print(f"Best fluorophore: {best_f}\nScore: {best_score}")
    print("\n\n")
