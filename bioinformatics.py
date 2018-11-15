from collections import Counter
from itertools import product
from typing import List, Tuple


DNA_COMPLEMENTS = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}


NUCLEOTIDE_INT_MAP = {
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3
}


def read_genome(file_path: str, has_header=False, has_footer=False) -> Tuple[str, str, str]:
    """ Read in genome file, and keep header/footer separate if necessary.

    Parameters
    ----------
    file_path
        File we want to read genome from
    has_header
        Whether or not file has a header, assumes 1 line
    has_footer
        Whether or not file has a footer, assumes 1 line
    Returns
    -------
    tuple
        index 0 is the genome, index 1 is the header, index 2 is the footer
    """
    with open(file_path, 'r') as infile:
        header = ''
        footer = ''
        body = infile.read().splitlines()

        if has_header:
            header = body.pop(0)
        if has_footer:
            footer = body.pop(-1)

        return ''.join(body), header, footer


def pattern_count(genome: str, pattern: str) -> int:
    """ Count number of times a pattern occurs in text

    Parameters
    ----------
    genome
        The genome sequence we want to search through
    pattern
        The pattern we want to search for in the genome

    Returns
    -------
    int
        Number of times pattern occurred in genome
    """
    count = 0
    for i in range(len(genome)):
        if genome[i:i+len(pattern)] == pattern:
            count += 1
    return count


def pattern_to_number(pattern: str) -> int:
    """Convert nucleotide pattern into base 4 representation where:
        A = 0, C = 1, G = 2, T = 3

    Parameters
    ----------
    pattern
        Pattern we want to convert into integer representation

    Returns
    -------
    int
        Integer representation of given pattern
    """
    base = 4  # nucleotide map is base 4
    mapped = pattern
    for nuc, num in NUCLEOTIDE_INT_MAP.items():
        mapped = mapped.replace(nuc, str(num))
    return int(mapped, base)


def number_to_pattern(number: int, k: int) -> str:
    """ Convert base 4 representation of nucleotide pattern back into the pattern
    A = 0, C = 1, G = 2, T = 3

    Parameters
    ----------
    number
        The base 4 representation of a nucleotide sequence
    k
        The length of the nucleotide sequence

    Returns
    -------
    str
        nucleotide sequence
    """
    base = 4
    new_str = ''
    while number > 0:
        new_str += str(number % base)
        number = number // base
    for nuc, num in NUCLEOTIDE_INT_MAP.items():
        new_str = new_str.replace(str(num), nuc)
    while len(new_str) < k:
        new_str = new_str + 'A'
    return new_str[::-1]


def reverse_complement(pattern: str) -> str:
    """ Get the complement of a nucleotide sequence and reverse it
    e.g. ACTG -> TGAC

    Parameters
    ----------
    pattern
        Nucleotide sequence we want the reverse complement of

    Returns
    -------
    str
        The reverse complement of a nucleotide sequence
    """
    complement = ''
    for nuc in pattern:
        complement += DNA_COMPLEMENTS[nuc]
    return complement[::-1]


def pattern_match(genome: str, pattern: str) -> List[int]:
    """ Get indices for start location of all matching patterns in genome

    Parameters
    ----------
    genome
        The genome sequence we want to search through
    pattern
        The pattern whose indices we wish to find

    Returns
    -------
    list
        The index of starting position for all occurrences of *pattern* in *genome*
    """
    indices = []
    for i in range(len(genome)):
        if genome[i:i+len(pattern)] == pattern:
            indices.append(i)
    return indices


def compute_frequencies(genome: str, k: int) -> List[int]:
    """ Make frequency array for each possible pattern of length *k*

    Parameters
    ----------
    genome
        Nucleotide sequence to create frequency array for
    k
        length of kmers

    Returns
    -------
    List(int)
        Frequencies of each kmer alphabetically indexed ('AA', 'AC', 'AG'...)
    """
    freq = Counter()
    for i in range(len(genome) - k + 1):
        pattern = genome[i: i + k]
        freq.update([str(pattern)])
    order = [freq[''.join(list(i))] for i in product('ACGT', repeat=k)]
    return order


def get_kmer_counts(genome: str, k: int) -> Counter:
    """ Count how many times each kmer appears in *genome*

    Parameters
    ----------
    genome
       The genome sequence get counts for
    k
        Length of kmers to get counts of

    Returns
    -------
    Counter
        Counter of how many times each pattern appeared
    """
    frequency = Counter()
    for i in range(len(genome) - k + 1):
        pattern = genome[i: i + k]
        frequency.update([str(pattern)])

    return frequency


def get_frequent(frequency: Counter, min_frequency: int) -> List[str]:
    """ Get all patterns in Counter that occur more than specified minimum

    Parameters
    ----------
    frequency
        Counter of frequencies
    min_frequency
        Minimum number of times an item needs to be in counter to be selected

    Returns
    -------
    List(str)
        Items that occurred at least the minimum number of times
    """
    frequent = []  # list of all strings occurring >= *min_frequency* times

    # since freq.most_common is sorted max->min, grab until count < min_freq
    for counts in frequency.most_common():
        if counts[1] < min_frequency:
            break
        frequent.append(counts[0])
    return frequent


def frequent_words(genome: str, k: int, min_frequency: int) -> List[str]:
    """ Get all kmers in genome that appear a minimum number of times

    Parameters
    ----------
    genome
        The genome sequence to get counts for
    k
        Length of kmers to get counts of
    min_frequency
        Minimum number of times a kmer needs to appear

    Returns
    -------
    List(str)
        List of kmers that occurred at least the min number of times
    """
    frequency = get_kmer_counts(genome, k)

    freq_str = []  # list of all strings occurring >= *min_frequency* times
    for counts in frequency.most_common():
        if counts[1] < min_frequency:
            break
        freq_str.append(counts[0])

    return list(set(freq_str))


def most_frequent_words(genome: str, k: int) -> List[str]:
    """ Get only the most frequently occurring kmers in genome

    Parameters
    ----------
    genome
        The genome sequence to get counts for
    k
        Length of kmers to get counts of

    Returns
    -------
    List(str)
        Most frequently occurring kmers
    """
    frequency = get_kmer_counts(genome, k)
    max_freq = frequency.most_common(1)[0][1]

    return get_frequent(frequency, max_freq)


def find_clumps(genome, k, L, t) -> List[str]:
    """ Find regularly occurring kmers in each clump of length *L*
    genome: str
        The genome to search
    k: int
        Length of each kmer to check
    L: int
        Length of a clump to search in
    t: int
        Minimum number of times a kmer must appear to be considered

    Returns
    -------
    list(str)
        All kmers that appear at least *t* times in a clump of size *L*
    """

    # initialize counter with first clump
    counter = get_kmer_counts(genome[0: L], k)
    # define beginning and end k-mer for first clump
    first_kmer = genome[0: k]
    frequent = get_frequent(counter, t)

    for i in range(1, len(genome) - L + 1):  # start at index 1 since we initialized with the first clump already
        # grab our next clump
        clump = genome[i: i + L]

        # remove first k-mer from previous clump
        counter.subtract([first_kmer])

        first_kmer = clump[0: k]
        last_kmer = clump[L - k: L]

        # add last k-mer from current clump
        counter.update([last_kmer])

        if counter[last_kmer] >= t:
            frequent.append(last_kmer)

    return list(set(frequent))
