def fitting_alignment(s: str, t: str, match_reward: float, mismatch_penalty: float, indel_penalty: float) -> tuple[float, str, str]:
    """
    A function that finds the best fitting alignment between two strings
    :param s: the string that can be partially matched
    :param t: the string that must be fully matched
    :param match_reward: the reward for a character match
    :param mismatch_penalty: the penalty for a character mismatch
    :param indel_penalty: the penalty for an insertion or deletion in the alignment
    :return: the total alignment score, the alignment for string s, the alignment for string t
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
                for k in range(len(s)):
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


def hamming(s: str, t: str) -> int:
    """
    A function that calculates the Hamming distance between two strings
    :param s: a string
    :param t: a string
    :return: the Hamming distance between s and t
    """
    count = 0
    for i in range(len(s)):
        if s[i] != t[i]:
            count += 1
    return count


def SeqSimWithHamming(s: str, t: str, n: int, match_reward: float, mismatch_penalty: float) -> dict[str, float]:
    subs_and_scores = {}
    for i in range(len(s) - len(t) + 1):
        sub_str = s[i:i + len(t)]
        ham = hamming(sub_str, t)
        score = match_reward*(len(t) - ham) - mismatch_penalty*ham
        subs_and_scores.update({sub_str: score})
    subs_and_scores = dict(sorted(subs_and_scores.items(), key=lambda x: x[1], reverse=True)[:n])
    return subs_and_scores


def SeqSim(s: str, f: str, n: int, match_reward: float = 1, mismatch_penalty: float = 1, indel_penalty: float = float('-inf')) -> dict[str, float]:
    """
    A function that calculates the sequence similarity score between all substrings of s and a string f and returns the
    top n substrings along with the scores
    :param s: a string of interest
    :param f: a string to compare to s
    :param n: the number of top substrings of s to return
    :param match_reward: the reward for a match in the fitting alignment
    :param mismatch_penalty: the penalty for having a mismatch in the fitting alignment
    :param indel_penalty: the penalty for having an insertion or deletion in the fitting alignment
    :return: a dictionary of the top n substrings with the best fitting alignment scores along with their scores
    """

    return SeqSimWithHamming(s, f, n, match_reward, mismatch_penalty)
