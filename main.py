import csv
import time
import matplotlib.pyplot as plt
from source import Score


# convert a DNA sequence to an RNA sequence
def convert_DNA_to_RNA(seq: str) -> str:
    convert = {'T': 'U', 'A': 'A', 'C': 'C', 'G': 'G', 'U': 'U'}
    return_str = ""
    for n in seq:
        if n not in convert:
            raise RuntimeError
        return_str += convert[n]

    return return_str


def plot(xs: list[float], ys: dict[str, list[float]], x_label: str, y_label: str, title: str):
    colors = []
    for i in range(len(xs)):
        x_vals = xs[i]
        y_label, y_vals = list(ys.items())[i]
        plt.plot(x_vals, y_vals, linestyle='-', color=colors[i % len(colors)], label=y_label)

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.legend()
    plt.show()


# Read fluorescent aptamers from CSV
fluorescent_seqs = {}
with open('data/fluorescent_aptamers.csv', mode='r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        fluorescent_seqs[row['name']] = convert_DNA_to_RNA(row['sequence'])

# Read genes from CSV
genes = {}
with open('data/genes.csv', mode='r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        genes[row['name']] = convert_DNA_to_RNA(row['sequence'])

top_subs = 10
a = 0.5

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
        f_name, f_seq = list(fluorescent_seqs.items())[i]
        sub_str, score = Score.Score(g_seq, f_seq, top_subs, a)
        if best_score < score:
            best_score, best_f = score, f_name
        print(f"For fluorophore sequence {i + 1}:")
        print(f"Optimal mRNA substring in gene: {sub_str}\tScore: {score}")

    print(f"Best fluorophore: {best_f}\nScore: {best_score}")
    print("\n\n")
