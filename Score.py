from typing import List, Dict, Iterable
import StructSim
import SeqSim


def Score(s: str, F: List[str], a = 1, b = 1):
    return a * SeqSim.SeqSim(s, F) + b * StructSim.StructSim(s, F)