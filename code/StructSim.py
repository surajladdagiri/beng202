from RNAFold import RNAFold


def StructSim(s: str, f: str) -> float:
    """
    Function that calculates the structural similarity between each ith element pairs of two list of RNA sequences
    :param s: an RNA sequence of length k
    :param f: an RNA sequence of length k
    :return: the max Jaccard index between each secondary structure of s and each secondary structure of f
    """

    if len(s) != len(f):
        raise RuntimeError
    RNAs = RNAFold(s)
    RNAf = RNAFold(f)
    j_index = 0
    for ss_s in RNAs:
        for ss_f in RNAf:
            intersect = []
            for pair in ss_s:
                if pair in ss_f:
                    intersect.append(pair)
            curr_j = len(intersect) / (len(ss_s) + len(ss_f) - len(intersect))
            j_index = max(j_index, curr_j)

    return j_index
