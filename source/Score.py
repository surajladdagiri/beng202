import StructSim
import SeqSim


def Score(g: str, f: str, n: int, a: int = 0.5) -> tuple[str, float]:
    seq_scores = SeqSim.SeqSim(g, f, n)
    struc_score = StructSim.StructSim(list(seq_scores.keys()), f)
    total_score = {sub: (a * seq_scores[sub] + (1-a) * struc_score[sub]) for sub in struc_score.keys()}

    return sorted(total_score.items(), key=lambda x: x[1], reverse=True)[0]
