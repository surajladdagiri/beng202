from RNAFold import RNAFold


def StructSim(S: list[str], F: list[str]) -> float:
    """
    Function that calculates the structural similarity between each ith element pairs of two list of RNA sequences
    :param S: a set of RNA sequences of a length k
    :param F: a set of RNA sequences of the same length k
    :return: the sum of the max Jaccard index between each secondary structures of s_i and each secondary structure of
    the corresponding f_i in F
    """

    total_index = 0
    if len(S) != len(F):
        raise RuntimeError
    for i in range(len(S)):
        RNAsi = RNAFold(S[i])
        RNAfi = RNAFold(F[i])
        j_index = 0
        for ss_s in RNAsi:
            for ss_f in RNAfi:
                intersect = []
                for pair in ss_s:
                    if pair in ss_f:
                        intersect.append(pair)
                curr_j = len(intersect) / (len(ss_s) + len(ss_f) - len(intersect))
                j_index = max(j_index, curr_j)
        total_index += j_index

    return total_index
