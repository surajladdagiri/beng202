import RNA  # pip install viennarna


def RNAFold(s: str) -> list[tuple[int, int]]:
    """
    A function to that returns the base pair interactions in an RNA secondary structure given the RNA sequence
    :param s: an RNA sequence consisting of RNA nucleotides {A, C, G, U}
    :return: the RNA secondary structure of s as tuples of interacting base pairs (index starts at 1)
    """
    # Obtain RNA secondary structure in dot-bracket notation
    dot_bracket, mfe = RNA.fold(s)

    # Stack to track base pair positions
    stack = []
    base_pairs = []

    # Iterate through structure to find matching parentheses
    for i, char in enumerate(dot_bracket, start=1):
        if char == '(':  # Opening bracket, push position onto stack
            stack.append(i)
        elif char == ')':  # Closing bracket, pop from stack and record pair
            if stack:
                j = stack.pop()
                base_pairs.append((j, i))

    return base_pairs
