import StructSim
import SeqSim


def Score(g: str, F: list[str], a: int = 1, b: int = 1) -> float:
    seq_score, seq_align = SeqSim.SeqSim(g, F)
    struc_score = StructSim.StructSim(seq_align, F)
    return a * seq_score + b * struc_score
