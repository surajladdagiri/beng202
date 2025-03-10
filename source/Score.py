from source import StructSim, SeqSim


def Score(g: str, f: str, n: int, a: float = 0.5, enhanced: bool = False) -> tuple[str, float]:
    seq_scores = SeqSim.SeqSim(g, f, n)
    struc_score = StructSim.StructSim(list(seq_scores.keys()), f, enhanced)
    total_score = {sub: (a * seq_scores[sub] + (1-a) * struc_score[sub]) for sub in struc_score.keys()}

    return sorted(total_score.items(), key=lambda x: x[1], reverse=True)[0]
