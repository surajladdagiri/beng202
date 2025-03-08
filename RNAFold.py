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

def rna_folding(s: str):
    n = len(s)
    dp = [[0] * n for _ in range(n)]
    backtrack = [[None] * n for _ in range(n)]  # Store pairing decisions
    pairs = {('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C')}

    # Filling the DP table
    for length in range(2, n+1):  # Consider substrings of increasing lengths
        for i in range(n - length + 1):
            j = i + length - 1

            # Case 1: Skip the current nucleotide
            dp[i][j] = dp[i+1][j]
            backtrack[i][j] = (i+1, j)  # Move to the next nucleotide
            
            # Case 2: Try to pair i with a valid k (where i < k ≤ j)
            for k in range(i+2, j+1):  # Ensure k - i > 1
                if (s[i], s[k]) in pairs:
                    left_part = dp[i+1][k-1] if k-1 >= i+1 else 0
                    right_part = dp[k+1][j] if k+1 <= j else 0
                    val = 1 + left_part + right_part

                    if val > dp[i][j]:
                        dp[i][j] = val
                        backtrack[i][j] = (i+1, k-1, k+1, j, i, k)  # Store the optimal pairing

    # Reconstruct the actual pairs
    def traceback(i, j):
        if i >= j:
            return []
        if backtrack[i][j] == (i+1, j):
            return traceback(i+1, j)  # Move to the next nucleotide
        else:
            i1, j1, i2, j2, pi, pk = backtrack[i][j]
            return [(pi, pk)] + traceback(i1, j1) + traceback(i2, j2)

    return traceback(0, n-1)

# Example usage
rna_string = "AUGCGAU"
print(rna_folding(rna_string))  # Output: List of (i, j) pairs


### CHECK OUT Nussinov Algorithm for RNA structure prediction