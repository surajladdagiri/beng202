from source import RNAFold
from source import RNAFoldEnhanced


def StructSim(S: list[str], f: str, enhanced_fold: bool) -> dict[str, float]:
    """
    Function that calculates the structural similarity between an RNA sequence and each RNA sequence in a set of RNA
    sequences
    :param S: a set of RNA sequences of length k
    :param f: an RNA sequence of length k
    :param enhanced_fold: whether to use the enhanced folding algorithm or not
    :return: a dictionary of containing the sequences in S and their similarity scores
    """

    scores = {}
    for s in S:
        if len(s) != len(f):
            raise RuntimeError

        if enhanced_fold:
            RNAs = RNAFoldEnhanced.rna_folding_enhanced(s)
            RNAf = RNAFoldEnhanced.rna_folding_enhanced(f)
        else:
            RNAs = RNAFold.rna_folding(s)
            RNAf = RNAFold.rna_folding(f)
        intersect = []
        for pair in RNAs:
            if pair in RNAf:
                intersect.append(pair)
        j_index = len(intersect) / (len(RNAs) + len(RNAf) - len(intersect))
        scores[s] = j_index

    return scores
