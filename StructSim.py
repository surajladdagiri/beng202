from RNAFold import RNAFold


def StructSim(s: str, F: list[str]) -> float:
    """
    Function that calculates the structural similarity between an RNA sequence and a list of RNA sequences
    :param s: an RNA sequence
    :param F: a set of RNA sequences
    :return: the sum of the Jaccard indices between s and each RNA sequence in F
    """
    RNAF = [RNAFold(f) for f in F]
    RNAs = RNAFold(s)

    total_index = 0
    for RNAf in RNAF:
        intersect = []
        for pair in RNAs:
            if pair in RNAf:
                intersect.append(pair)
        j_index = len(intersect) / (len(RNAs) + len(RNAf) - len(intersect))
        total_index += j_index

    return total_index
