import StructSim
import SeqSim


def Score(s: str, F: list[str], a: int = 1, b: int = 1) -> float:
    return a * SeqSim.SeqSim(s, F) + b * StructSim.StructSim(s, F)
