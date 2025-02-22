from typing import List, Dict, Iterable
import RNAFold

def StructSim(s: str, F: List[str]):
    RNAf = [RNAFold.RNAFold(f) for f in F]
    
    pass