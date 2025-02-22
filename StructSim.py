import RNAFold


def StructSim(s: str, F: list[str]) -> float:
    RNAF = [RNAFold.RNAFold(f) for f in F]
    RNAs = RNAFold.RNAFold(s)

    total_index = 0
    for RNAf in RNAF:
        intersect = []
        for pair in RNAs:
            if pair in RNAf:
                intersect.append(pair)
        j_index = len(intersect) / (len(RNAs) + len(RNAf) - len(intersect))
        total_index += j_index

    return total_index
