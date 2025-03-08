#  main file to run our code
import Score

fluorescent_seqs = []
genes = []

if fluorescent_seqs:
    k = len(fluorescent_seqs[0])
else:
    raise RuntimeError("Please pass in some fluorescent sequences!")

for i in range(len(fluorescent_seqs)):
    if fluorescent_seqs[i] != k:
        raise RuntimeError("All fluorescent sequences must be the same length!")

if not genes:
    raise RuntimeError("Please pass in some gene sequences!")

best_score = float('-inf')
best_gs = []
for i, g in enumerate(genes, start=1):
    curr_score = Score.Score(g, fluorescent_seqs)
    print(f"Score for gene {i}: {curr_score}")
    if curr_score > best_score:
        best_score = curr_score
        best_gs = [i]
    elif curr_score == best_score:
        best_gs.append(i)

print("Best overall score: {}".format(best_score))
print("\nGenes with that score: ")
for g_i in best_gs:
    print(g_i)
