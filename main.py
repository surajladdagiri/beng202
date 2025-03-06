#  main file to run our code
import Score

fluorescent_seqs = []

if fluorescent_seqs:
    k = len(fluorescent_seqs[0])
else:
    raise RuntimeError("Please pass in some input sequences!")

for i in range(len(fluorescent_seqs)):
    if fluorescent_seqs[i] != k:
        raise RuntimeError("All input sequences must be the same length!")


def kmer_set(k: int) -> list[str]:
    """
    Generates all possible RNA nucleotide kmers of a given length k
    :param k: length of kmers to generate
    :return: list of all possible RNA kmers
    """
    nucs = ['A', 'U', 'C', 'G']
    kmers = []
    if k == 0:
        return kmers
    else:
        sub_kmers = kmer_set(k - 1)
        for kmer in sub_kmers:
            new_kmers = [n + kmer for n in nucs]
            kmers += new_kmers
        return kmers


best_kmers = []
best_score = float('inf')

for kmer in kmer_set(k):
    curr_score = Score.Score(kmer, fluorescent_seqs)
    if curr_score > best_score:
        best_score = curr_score
        best_kmers = [kmer]
    elif curr_score == best_score:
        best_kmers.append(kmer)

print("Best overall score: {}".format(best_score))
print("\nKmers with that score: ")
for kmer in best_kmers:
    print(kmer)
