import random


def RandomData(k: int) -> str:
    """
    Generate a random DNA sequence of length k
    :param k: the length of the DNA sequence
    :return: a random DNA sequence of length k
    """
    return ''.join(random.choices('ACGU', k=k))
