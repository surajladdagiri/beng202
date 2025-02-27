from typing import List, Dict, Iterable, Tuple
import numpy as np


def fitting_alignment(s: str, t: str, match_reward: int, mismatch_penalty: int, indel_penalty: int) -> Tuple[int, str, str]:
    """
    Compute the fitting alignment of two strings given a match reward, mismatch penalty, and indel penalty.
    String t is fit to s (all of t is aligned to s)
    """
    cache = {}

    def count_and_path(i, j):
        if (i, j) in cache:
            return cache[(i, j)]
        max_count = float('-inf')
        max_first = ""
        max_second = ""
        if i == 0 and j == 0:
            max_count = 0
        else:
            if i == len(s) and j == len(t):
                for k in range(len(s) + 1):
                    if k == len(s) and j == len(t):
                        continue
                    count, first, second = count_and_path(k, j)
                    if count >= max_count:
                        max_count = count
                        max_first = first
                        max_second = second
            elif j == 0:
                max_count = 0
            if i != 0 and j != 0:
                diag_count, diag_first, diag_second = count_and_path(i - 1, j - 1)
                match_score = match_reward if s[i - 1] == t[j - 1] else 0 - mismatch_penalty
                diag_count += match_score
                if diag_count > max_count:
                    max_count = diag_count
                    max_first = diag_first + s[i - 1]
                    max_second = diag_second + t[j - 1]
            if i != 0:
                down_count, down_first, down_second = count_and_path(i - 1, j)
                down_count -= indel_penalty
                if down_count > max_count:
                    max_count = down_count
                    max_first = down_first + s[i - 1]
                    max_second = down_second + '-'
            if j != 0:
                right_count, right_first, right_second = count_and_path(i, j - 1)
                right_count -= indel_penalty
                if right_count >= max_count:
                    max_count = right_count
                    max_first = right_first + '-'
                    max_second = right_second + t[j - 1]

        info = (max_count, max_first, max_second)
        cache[(i, j)] = info
        return info
    info = count_and_path(len(s), len(t))
    return info


def SeqSim(s: str, F: List[str], match_reward: int = 1, mismatch_penalty: int = 1, indel_penalty: int = 1) -> int:
    score = 0
    for f in F:
        curr_score, f_align, s_align = fitting_alignment(f, s, match_reward, mismatch_penalty, indel_penalty)
        score += curr_score
    return score
