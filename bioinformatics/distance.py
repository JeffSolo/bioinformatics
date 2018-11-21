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
