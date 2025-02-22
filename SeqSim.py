from typing import List, Dict, Iterable, Tuple
import numpy as np

def local_alignment(match_reward: int, mismatch_penalty: int, indel_penalty: int, s: str, t: str):
    """
    Compute the local alignment of two strings based on match reward, mismatch penalty, and indel penalty.
    """
    score = np.zeros((len(s)+1, len(t)+1))
    
    max_score = 0
    max_pos = (0,0)
    
    for i in range(1,len(s)+1):
        for j in range(1,len(t)+1):
            if s[i-1] == t[j-1]:
                match = match_reward
            else:
                match = -mismatch_penalty
                
            score[i,j] = max(score[i-1, j]-indel_penalty, score[i,j-1]-indel_penalty, score[i-1,j-1]+match, 0)
            if score[i,j] > max_score:
                max_score = score[i,j]
                max_pos = (i,j)
     
    i = max_pos[0]
    j = max_pos[1]
    return score[i,j]


def SeqSim(s: str, F: List[str], match_reward: int = 1, mismatch_penalty: int = 1, indel_penalty: int = 1) -> int:
    score = 0
    for f in F:
        score += local_alignment(match_reward, mismatch_penalty, indel_penalty, s, f)
    return score