import StructSim
import SeqSim


def Score(g: str, f: str, a: int = 0.5) -> float:
    seq_score, g_sub = SeqSim.SeqSim(g, f)
    struc_score = StructSim.StructSim(g_sub, f)
    return a * seq_score + (1-a) * struc_score
