from source import RNAFold


def StructSim(S: list[str], f: str) -> dict[str, float]:
    """
    Function that calculates the structural similarity between an RNA sequence and each RNA sequence in a set of RNA
    sequences
    :param S: a set of RNA sequences of length k
    :param f: an RNA sequence of length k
    :return: a dictionary of containing the sequences in S and their similarity scores
    """

    scores = {}
    for s in S:
        if len(s) != len(f):
            raise RuntimeError
        RNAs = RNAFold.rna_folding(s)
        RNAf = RNAFold.rna_folding(f)
        intersect = []
        for pair in RNAs:
            if pair in RNAf:
                intersect.append(pair)
        j_index = len(intersect) / (len(RNAs) + len(RNAf) - len(intersect))
        scores[s] = j_index

    return scores
