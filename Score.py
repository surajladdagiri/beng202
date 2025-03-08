import StructSim
import SeqSim


def Score(g: str, f: str, a: int = 1, b: int = 1) -> float:
    seq_score, g_sub = SeqSim.SeqSim(g, f)
    struc_score = StructSim.StructSim(g_sub, f)
    return a * seq_score + b * struc_score
