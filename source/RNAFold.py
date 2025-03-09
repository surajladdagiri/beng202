def rna_folding(s: str):
    """
    Computes the optimal RNA secondary structure using DP.

    This is an implementation of the Nussinov Algorithm.

    Input:
    s (str): A string representing the RNA sequence consisting of characters {A,U,C,G}

    Output:
    List[Tuple[int, int]]: A list of tuples where each tuple (i, j) represents a pair of indices in the RNA sequence that form a base pair in the optimal secondary structure.
    """
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

    # Utils function to reconstruct the actual pairs
    def traceback(i, j):
        if i >= j:
            return []
        if backtrack[i][j] == (i+1, j):
            return traceback(i+1, j)  # Move to the next nucleotide
        else:
            i1, j1, i2, j2, pi, pk = backtrack[i][j]
            return [(pi, pk)] + traceback(i1, j1) + traceback(i2, j2)

    return traceback(0, n-1)

# #Test
# print(rna_folding("UAGCAUCGAUCGAUCUACGU"))