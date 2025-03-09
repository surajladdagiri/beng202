import csv
import time
import matplotlib.pyplot as plt

from source import Score
from source import RandomData


# convert a DNA sequence to an RNA sequence
def convert_DNA_to_RNA(seq: str) -> str:
    """
    Converts a DNA sequence to an RNA sequence
    :param seq: DNA sequence string
    :return: RNA sequence equivalent for seq
    """
    convert = {'T': 'U', 'A': 'A', 'C': 'C', 'G': 'G', 'U': 'U'}
    return_str = ""
    for n in seq:
        if n not in convert:
            raise RuntimeError
        return_str += convert[n]

    return return_str


def plot(xs: list[list[float]], ys: list[tuple[str, list[float]]], x_label: str, y_label: str, title: str) -> None:
    """
    Plots data using matplotlib
    :param xs: set of sets of x values to plot
    :param ys: corresponding y values to plot along with a title to use for the legend
    :param x_label: x-axis label for the plot
    :param y_label: y-axis label for the plot
    :param title: title for the plot
    :return: None
    """
    colors = ['red', 'blue', 'green', 'yellow', 'orange', 'purple', 'brown', 'pink', 'gray', 'black']
    for i in range(len(xs)):
        x_vals = xs[i]
        y_label, y_vals = ys[i]
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

def_k = 10
def_alpha = 0.5


def subs_and_scores(gs: [str, str], fs: dict[str, str], top_n: int, alpha: float) -> dict[str, tuple[str, float]]:
    return_dict = {}
    for g_name, g_seq in gs.items():
        print(f"Running for gene {g_name}")
        best_score = float('-inf')
        best_f = ""
        for i in range(len(fs)):
            f_name, f_seq = list(fs.items())[i]
            sub_str, score = Score.Score(g_seq, f_seq, top_n, alpha)
            if best_score < score:
                best_score, best_f = score, f_name
            print(f"For fluorophore sequence {i + 1}:")
            print(f"Optimal mRNA substring in gene: {sub_str}\tScore: {score}")

        print(f"Best fluorophore: {best_f}\nScore: {best_score}")
        print("\n\n")
        return_dict[g_name] = (best_f, best_score)

    return return_dict


# Testing k vs. score
score_save = {}
ks = [1, 5, 10, 15, 25, 50]
for k in ks:
    f_test = {"mangoI": fluorescent_seqs["mangoI"]}
    g_test = {"epsC": genes['epsC']}
    f, score = subs_and_scores(g_test, f_test, k, def_alpha)['epsC']
    score_save[k] = score

print(score_save)


# Testing alpha vs. score
score_save = {}
alphs = [0, 0.001, 0.01, 0.1, 0.2, 0.3, 0.5, 0.75, 0.9, 1]
for alph in alphs:
    f_test = {"mangoI": fluorescent_seqs["mangoI"]}
    g_test = {"epsC": genes['epsC']}
    f, score = subs_and_scores(g_test, f_test, def_k, alph)['epsC']
    score_save[alph] = score

print(score_save)


# Testing length of g vs. runtime
time_save = {}
g_lens = []
for g_len in g_lens:
    random_g = RandomData.RandomData(g_len)
    f_test = {"mangoI": fluorescent_seqs["mangoI"]}
    g_test = {"randomGene": random_g}
    start_time = time.time()
    subs_and_scores(g_test, f_test, def_k, def_alpha)
    end_time = time.time()
    time_save[g_len] = end_time - start_time

print(time_save)


# Testing k vs. runtime
time_save = {}
ks = [1, 5, 10, 15, 25, 50, 100]
for k in ks:
    random_g = RandomData.RandomData(100000)
    f_test = {"mangoI": fluorescent_seqs["mangoI"]}
    g_test = {"randomGene": random_g}
    start_time = time.time()
    subs_and_scores(g_test, f_test, k, def_alpha)
    end_time = time.time()
    time_save[k] = end_time - start_time

print(time_save)
