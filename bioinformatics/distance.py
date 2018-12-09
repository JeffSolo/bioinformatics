from typing import List


def hamming_distance(p: str, q: str) -> int:
    """ Compute hamming distance between two strings

    Parameters
    ----------
    p : str
        first string to be compared
    q : str
        second string to be compared

    Returns
    -------
    int
        Distance between the strings
    """
    if len(p) != len(q):
        raise ValueError("Hamming Distance requires strings to be of equal length")
    distance = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            distance += 1

    return distance


def get_neighborhood(pattern: str, alphabet: str, max_distance: int) -> List[str]:
    """ Get all patterns within hamming distance *max_distance* using characters in *alphabet*
        e.g. get_neighbors('AB', 'ABC', 1) -> ['AA', 'AB', 'BB', 'CB', 'AC']

    Parameters
    ----------
    pattern: str
        Pattern we want to get neighbors of
    alphabet : str
        All allowable characters that can be found in pattern
    max_distance: int
        Maximum hamming distance of neighbors

    Returns
    -------
    list
        List of all patterns within specified distance
    """
    single_nucs = tuple(alphabet)
    if max_distance == 0:
        return [pattern]
    if len(pattern) == 1:
        return single_nucs
    neighborhood = []
    suffix = pattern[1:]
    suffix_neighbors = get_neighborhood(suffix, alphabet, max_distance)
    for neighbor in suffix_neighbors:
        if hamming_distance(suffix, neighbor) < max_distance:  # less than since we'll add another nucleotide
            for nuc in single_nucs:
                neighborhood.append(nuc + neighbor)
        else:
            neighborhood.append(pattern[0] + neighbor)
    return neighborhood
