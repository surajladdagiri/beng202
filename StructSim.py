from RNAFold import RNAFold


def StructSim(s: str, F: list[str]) -> float:
    """
    Function that calculates the structural similarity between an RNA sequence and a list of RNA sequences
    :param s: an RNA sequence
    :param F: a set of RNA sequences
    :return: the sum of the max Jaccard index between each secondary structures of s and each secondary structure of f
    in F
    """

    total_index = 0
    RNAss = RNAFold(s)
    for f in F:
        j_index = 0
        RNAfs = RNAFold(f)
        for RNAf in RNAfs:
            for RNAs in RNAss:
                intersect = []
                for pair in RNAs:
                    if pair in RNAf:
                        intersect.append(pair)
                curr_j = len(intersect) / (len(RNAs) + len(RNAf) - len(intersect))
                j_index = max(j_index, curr_j)
        total_index += j_index

    return total_index
