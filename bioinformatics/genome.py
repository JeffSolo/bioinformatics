from collections import Counter
from itertools import chain, product
from typing import List
from .distance import hamming_distance


class Genome:
    """ A genome consisting of a string of nucleotides

    Parameters
    ----------
    sequence : str default ''
        Sequence of nucleotides, e.g. 'ACGTGACTAAGATCGGG'

    Attributes
    ----------
    dna_complement (class attribute) : dict
        Complement of each nucleobase, i.e. {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    _nucleotide_int_map (class attribute): dict
        Map of nucleotides to ints in alphabetical order i.e. {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    sequence : str
        Sequence of nucleotides that make up a genome, e.g. 'ACGTAATCGTTCGCCC'
    """

    dna_complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }

    _nucleotide_int_map = {
        'A': 0,
        'C': 1,
        'G': 2,
        'T': 3
    }

    def __init__(self, sequence: str = ''):
        self.sequence = sequence

    def read_genome(self, file_path: str, skip_header_rows=0, skip_footer_rows=0) -> str:
        """ Read in genome from file path (combines all lines into 1 string)

        Parameters
        ----------
        file_path : str
            File we want to read genome from
        skip_header_rows : int, optional default 0
            Do not read in the first n lines of file
        skip_footer_rows : int, optional default 0
            Do not read in the last n lines of file

        Returns
        -------
        str
            returns genome sequence
        """
        with open(file_path, 'r') as infile:
            body = infile.read().splitlines()

            for i in range(skip_header_rows):
                body.pop(0)
            for j in range(skip_footer_rows):
                body.pop(-1)

            self.sequence = ''.join(body)
            return self.sequence

    def reverse_complement(self) -> str:
        """ Get the complement of our genome sequence and reverse it
        e.g. ACTG -> CAGT

        Returns
        -------
        str
            The reverse complement of genome sequence
        """

        return self._get_reverse_complement(self.sequence)

    @classmethod
    def _get_reverse_complement(cls, pattern: str) -> str:
        """ Get reverse complement of a nucleotide sequence
        e.g. ACTG -> CAGT

        Parameters
        ----------
        pattern : str
            Nucleotide sequence

        Returns
        -------
        str
            Reverse complement of *pattern*
        """
        return ''.join(reversed([cls.dna_complement[nuc] for nuc in pattern]))

    def pattern_count(self, pattern: str, max_distance=0) -> int:
        """ Count number of times a pattern occurs in genome

        Parameters
        ----------
        pattern : str
            The pattern we want to count occurrences of
        max_distance : int, optional default 0
             Maximum allowable hamming distance from *pattern* to count as a match

        Returns
        -------
        int
            Number of times pattern occurred
        """
        count = 0
        for i in range(len(self.sequence) - len(pattern) + 1):
            sub_sequence = self.sequence[i:i+len(pattern)]
            if sub_sequence == pattern or hamming_distance(sub_sequence, pattern) <= max_distance:
                count += 1

        return count

    def pattern_match_index(self, pattern: str, max_distance=0) -> List[int]:
        """ Get indices for start location of all matching patterns in genome

        Parameters
        ----------
        pattern : str
            The pattern whose indices we wish to find
        max_distance : int, optional default 0
            Maximum allowable hamming distance from *pattern* to count as a match

        Returns
        -------
        list
            The index of starting position for all occurrences of *pattern* in *genome*
        """
        indices = []
        for i in range(len(self.sequence) - len(pattern) + 1):
            sub_sequence = self.sequence[i:i + len(pattern)]
            if sub_sequence == pattern or hamming_distance(sub_sequence, pattern) <= max_distance:
                indices.append(i)
        return indices

    def compute_all_frequencies_alphabetically(self, k: int, max_distance=0) -> List[int]:
        """ Make frequency array for each possible pattern of length *k*, alphabetically indexed

        Parameters
        ----------
        k : int
            length of kmers
        max_distance : int, optional default 0
            Maximum hamming distance from a pattern to count as a match

        Returns
        -------
        list
            Frequencies of each kmer alphabetically indexed.
            Index 0 is frequency of AA, index 1 is frequency of AC, index 2 is frequency of AG...
        """
        freq = self.get_kmer_counts(self.sequence, k, max_distance)
        order = [freq[''.join(list(i))] for i in product('ACGT', repeat=k)]
        return order

    @classmethod
    def get_kmer_counts(cls, sub_sequence: str, k: int, count_reverse_complement=False, max_distance=0) -> Counter:
        """ Count how many times each kmer appears in sub_sequence

        Parameters
        ----------
        sub_sequence : str
            Part or all of a genome sequence to get kmer counts for
        k : int
            Length of kmers to get counts of
        count_reverse_complement : bool, optional default False
            Whether we also want to add occurrences of the (approximate) reverse complement to our frequency
        max_distance : int, optional default 0
            Maximum hamming distance from a pattern to count as a match

        Returns
        -------
        Counter
            Counter of how many times each kmer occurred
        """
        def get_frequencies(pattern):
            freq = Counter()

            for i in range(len(sub_sequence) - k + 1):
                neighbors = cls.get_neighbors(pattern[i: i + k], max_distance)
                for neighbor in neighbors:
                    freq.update([str(neighbor)])

            return freq

        frequency = get_frequencies(sub_sequence)

        if count_reverse_complement:
            rc = cls._get_reverse_complement(sub_sequence)
            rc_frequency = get_frequencies(rc)
            frequency.update(rc_frequency)

        return frequency

    def frequent_kmers(self, k: int, min_frequency=0, count_reverse_complement=False, max_distance=0) -> List[str]:
        """ Get all kmers in genome that appear a minimum number of times

        Parameters
        ----------
        k : int
            Length of kmers to get counts of
        min_frequency : int
            Minimum number of times a kmer needs to appear
        count_reverse_complement : bool, optional default False
            Whether we also want to add occurrences of the (approximate) reverse complement to our frequency
        max_distance : int, optional default 0
            Maximum allowable hamming distance to be considered a matching kmer

        Returns
        -------
        list
            List of kmers that occurred at least the min number of times
        """
        frequency = self.get_kmer_counts(self.sequence, k, count_reverse_complement, max_distance)
        return self.get_frequent_kmer(frequency, min_frequency)

    def most_frequent_kmer(self, k: int, max_distance=0, count_reverse_complement=False) -> List[str]:
        """ Get only the most frequently occurring kmers in genome

        Parameters
        ----------
        k : int
            Length of kmers to get counts of
        count_reverse_complement : bool, optional default False
            Whether we also want to add occurrences of the (approximate) reverse complement to our frequency
        max_distance : int, optional default 0
            Maximum allowable hamming distance to be considered a matching kmer

        Returns
        -------
        list
            Most frequently occurring kmers
        """
        frequency = self.get_kmer_counts(self.sequence, k, count_reverse_complement, max_distance)
        max_freq = frequency.most_common(1)[0][1]

        return self.get_frequent_kmer(frequency, max_freq)

    def find_clumps(self, k: int, L: int, t: int) -> List[str]:
        """ Find regularly occurring kmers in each clump of length *L*
        k : int
            Length of each kmer to check
        L : int
            Length of a clump to search in
        t : int
            Minimum number of times a kmer must appear to be considered

        Returns
        -------
        list
            All kmers that appear at least *t* times in a clump of size *L*
        """

        # initialize counter with first clump
        counter = self.get_kmer_counts(self.sequence[0: L], k)
        # define beginning and end k-mer for first clump
        first_kmer = self.sequence[0: k]
        frequent = self.get_frequent_kmer(counter, t)

        for i in range(1, len(self.sequence) - L + 1):  # start at index 1 since we initialized with the first clump
            # grab our next clump
            clump = self.sequence[i: i + L]

            # remove first k-mer from previous clump
            counter.subtract([first_kmer])

            first_kmer = clump[0: k]
            last_kmer = clump[L - k: L]

            # add last k-mer from current clump
            counter.update([last_kmer])

            if counter[last_kmer] >= t:
                frequent.append(last_kmer)

        return list(set(frequent))

    def minimum_skew(self) -> List[int]:
        """ Find the indices with the minimum skew
        When we encounter C decrease by 1, G increase by 1
        Helps find origin of replication via mutations due to how DNA replicates

        Returns
        -------
        list
            indices where the skew is lowest
        """
        skew = 0
        steps = [0]
        for nucleotide in self.sequence:
            if nucleotide == 'C':
                skew -= 1
            elif nucleotide == 'G':
                skew += 1
            steps.append(skew)

        min_skew = min(steps)
        return [i for i, val in enumerate(steps) if val == min_skew]

    @classmethod
    def median_string(cls, dna_strands: List[str], k: int) -> List[str]:
        """ Find kmers minimizing hamming distance amongst all dna strands

        Parameters
        ----------
        dna_strands : list
            DNA strands
        k : int
            length of kmers

        Returns
        -------
        list
            kmers with minimum distance
        """
        min_distance = k * len(dna_strands)
        median = []
        for pattern in [''.join(i) for i in product('ACGT', repeat=k)]:
            distance = cls.distance_between_pattern_and_strands(dna_strands, pattern)
            if distance < min_distance:
                min_distance = distance
                median = [pattern]
            elif min_distance == distance:
                median.append(pattern)
        return median

    @classmethod
    def motif_enumeration(cls, dna_strands: List[str], k: int, max_distance=0) -> List[str]:
        """ Brute force method for finding motifs that occur in dna strands

        Parameters
        ----------
        dna_strands : list(str)
            DNA strands
        k : int
            Length of kmer
        max_distance : int, optional default 0
            Maximum hamming allowable distance for a motif

        Returns
        -------
        List
            All (k, d) motifs in dna strand
        """
        neighbors_list = []
        for strand in dna_strands:
            temp = []
            for i in range(len(strand) - k + 1):
                kmer = strand[i: i + k]
                temp.append(set(cls.get_neighbors(kmer, max_distance)))
            neighbors_list.append(list(chain.from_iterable(temp)))

        patterns = set(neighbors_list[0])
        for neighbors in neighbors_list[1:]:
            patterns = patterns & set(neighbors)

        return list(set(patterns))

    @staticmethod
    def distance_between_pattern_and_strands(dna_strand: List[str], pattern: str) -> int:
        """ Sum the hamming distance between a pattern and each dna strand

        Parameters
        ----------
        dna_strand : list
            List of dna strands
        pattern : str
            Pattern to check distance against

        Returns
        -------
        int
            Total distance between pattern and strands
        """
        k = len(pattern)
        total_distance = 0

        for strand in dna_strand:
            distance = k
            for i in range(len(strand) - k + 1):
                kmer = strand[i: i+k]
                if distance > hamming_distance(pattern, kmer):
                    distance = hamming_distance(pattern, kmer)
            total_distance += distance
        return total_distance

    @staticmethod
    def get_frequent_kmer(frequency: Counter, min_frequency: int) -> List[str]:
        """ Get all items in Counter that occur more than specified minimum

        Parameters
        ----------
        frequency : Counter object
            Counter of frequencies
        min_frequency : int
            Minimum number of times an item needs to be in counter to be selected

        Returns
        -------
        list
            Items that occurred at least the minimum number of times
        """
        frequent = []  # list of all strings occurring >= *min_frequency* times

        # since freq.most_common is sorted max->min, grab until count < min_freq
        for counts in frequency.most_common():
            if counts[1] < min_frequency:
                break
            frequent.append(counts[0])
        return frequent

    @classmethod
    def get_neighbors(cls, pattern: str, max_distance: int) -> List[str]:
        """ Get all nearby patterns within hamming distance *max_distance*

        Parameters
        ----------
        pattern:
            Pattern we want to get neighbors of

        max_distance:
            Maximum hamming distance of neighbors

        Returns
        -------
        list
            List of all patterns within specified distance
        """
        single_nucs = list(cls._nucleotide_int_map.keys())
        if max_distance == 0:
            return [pattern]
        if len(pattern) == 1:
            return single_nucs
        neighborhood = []
        suffix = pattern[1:]
        suffix_neighbors = cls.get_neighbors(suffix, max_distance)
        for neighbor in suffix_neighbors:
            if hamming_distance(suffix, neighbor) < max_distance:  # less than since we'll add another nucleotide
                for nuc in single_nucs:
                    neighborhood.append(nuc + neighbor)
            else:
                neighborhood.append(pattern[0] + neighbor)
        return neighborhood

    @classmethod
    def pattern_to_number(cls, pattern: str, nucleotide_map: dict=None) -> int:
        """" Convert nucleotide pattern into base 4 representation where:

        Parameters
        ----------
        pattern: str
            Pattern we want to convert into number
        nucleotide_map : dict, optional default {'A': 0, 'C': 1, 'G': 2, 'T': 3}
            What nucleotide we want to map each number to, should be unique and < 4

        Returns
        -------
        int
            Base 4 representation of genome
        """
        if nucleotide_map is None:
            nucleotide_map = cls._nucleotide_int_map

        base = 4  # nucleotide map is base 4
        mapped = pattern
        for nuc, num in nucleotide_map.items():
            mapped = mapped.replace(nuc, str(num))
        return int(mapped, base)

    @classmethod
    def number_to_pattern(cls, number: int, k: int, nucleotide_map: dict=None) -> str:
        """ Convert base 4 representation of nucleotide pattern back into the pattern
        Parameters
        ----------
        number : int
            The base 4 representation of a nucleotide sequence
        k : int
            The length of the nucleotide sequence
        nucleotide_map : dict, optional default {'A': 0, 'C': 1, 'G': 2, 'T': 3}
            What nucleotide we want to map each number to, should be unique and < 4

        Returns
        -------
        str
            nucleotide sequence
        """
        if nucleotide_map is None:
            nucleotide_map = cls._nucleotide_int_map

        base = 4
        new_str = ''
        while number > 0:
            new_str += str(number % base)
            number = number // base
        for nuc, num in nucleotide_map.items():
            new_str = new_str.replace(str(num), nuc)
        while len(new_str) < k:
            new_str = new_str + 'A'
        return new_str[::-1]
