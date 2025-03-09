def dot_bracket_notation(pairs, s):
    """
    Converts a list of base pairs into dot-bracket notation.

    Args:
        pairs (List[Tuple[int, int]]): A list of tuples where each tuple (i, j) represents a pair of indices in the RNA sequence that form a base pair.
        s (str): The original RNA sequence.

    Returns:
        str: The dot-bracket notation of the RNA secondary structure.
    """
    n = len(s)
    dot_bracket = ['.'] * n

    for i, j in pairs:
        dot_bracket[i] = '('
        dot_bracket[j] = ')'

    return ''.join(dot_bracket)