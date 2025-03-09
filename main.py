import csv
from source import Score

# Read fluorescent aptamers from CSV
fluorescent_seqs = {}
with open('data/fluorescent_aptamers.csv', mode='r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        fluorescent_seqs[row['name']] = row['sequence']

# Read genes from CSV
genes = {}
with open('data/genes.csv', mode='r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        genes[row['name']] = row['sequence']

# Ensure mangoI and epsH are available
if "mangoI" not in fluorescent_seqs:
    raise RuntimeError("mangoI sequence not found in fluorescent aptamers!")

if "epsH" not in genes:
    raise RuntimeError("epsH sequence not found in genes!")

# Compare mangoI against epsH
genes = {name: seq for name, seq in genes.items() if name == "epsH"}
fluorescent_seqs = {name: seq for name, seq in fluorescent_seqs.items() if name == "mangoI"}

for g_name, g_seq in genes.items():
    print(f"Running for gene {g_name}")
    best_score = float('-inf')
    best_f = ""
    for i in range(len(fluorescent_seqs)):
        curr_f = fluorescent_seqs[i]
        sub_str, score = Score.Score(g_seq, curr_f)
        if best_score < score:
            best_score, best_f = score, curr_f
        print(f"For fluorophore sequence {i + 1}:")
        print(f"Optimal substring in gene: {sub_str}\tScore: {score}")

    print(f"Best fluorophore: {best_f}\nScore: {best_score}")
    print("\n\n")
