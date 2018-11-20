from collections import Counter
from itertools import product
from typing import List


class Genome:
    """ A genome consisting of a string of nucleotides

    Parameters
    ----------
    sequence : str default None
        Sequence of nucleotides, e.g. 'ACGTGACTAAGATCGGG'

    Attributes
    ----------
    __dna_complement (class attribute) : dict
        Complement of each nucleobase, i.e. {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    __nucleotide_int_map (class attribute): dict
        Map of nucleotides to ints in alphabetical order i.e. {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    sequence : str
        Sequence of nucleotides that make up a genome, e.g. 'ACGTAATCGTTCGCCC'
    """

    __dna_complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }

    __nucleotide_int_map = {
        'A': 0,
        'C': 1,
        'G': 2,
        'T': 3
    }

    def __init__(self, sequence: str = None):
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

    def pattern_count(self, pattern: str) -> int:
        """ Count number of times a pattern occurs in genome

        Parameters
        ----------
        pattern : str
            The pattern we want to count occurrences of

        Returns
        -------
        int
            Number of times pattern occurred
        """
        count = 0
        for i in range(len(self.sequence)):
            if self.sequence[i:i+len(pattern)] == pattern:
                count += 1
        return count

    def reverse_complement(self) -> str:
        """ Get the complement of a nucleotide sequence and reverse it
        e.g. ACTG -> CAGT

        Returns
        -------
        str
            The reverse complement of a nucleotide sequence
        """
        complement = ''
        for nuc in self.sequence:
            complement += self.__dna_complement[nuc]
        return complement[::-1]

    def pattern_match_index(self, pattern: str) -> List[int]:
        """ Get indices for start location of all matching patterns in genome

        Parameters
        ----------
        pattern : str
            The pattern whose indices we wish to find

        Returns
        -------
        list
            The index of starting position for all occurrences of *pattern* in *genome*
        """
        indices = []
        for i in range(len(self.sequence)):
            if self.sequence[i:i + len(pattern)] == pattern:
                indices.append(i)
        return indices

    # TODO return dictionary instead. Only reason it currently doesn't is because of course required output format
    def compute_frequencies(self, k: int) -> List[int]:
        """ Make frequency array for each possible pattern of length *k*, alphabetically indexed

        Parameters
        ----------
        k : int
            length of kmers

        Returns
        -------
        list
            Frequencies of each kmer alphabetically indexed.
            Index 0 is frequency of AA, 1 is frequency of AC, 2 is frequency of AG...
        """
        freq = Counter()
        for i in range(len(self.sequence) - k + 1):
            pattern = self.sequence[i: i + k]
            freq.update([str(pattern)])
        order = [freq[''.join(list(i))] for i in product('ACGT', repeat=k)]
        return order

    @staticmethod
    def get_kmer_counts(sub_sequence: str, k: int) -> Counter:
        """ Count how many times each kmer appears in sub_sequence

        Parameters
        ----------
        sub_sequence : str
            Part or all of a genome sequence to get kmer counts for
        k : int
            Length of kmers to get counts of

        Returns
        -------
        Counter
            Counter of how many times each kmer occurred
        """
        frequency = Counter()
        for i in range(len(sub_sequence) - k + 1):
            pattern = sub_sequence[i: i + k]
            frequency.update([str(pattern)])

        return frequency

    def frequent_kmers(self, k: int, min_frequency: int) -> List[str]:
        """ Get all kmers in genome that appear a minimum number of times

        Parameters
        ----------
        k : int
            Length of kmers to get counts of
        min_frequency : int
            Minimum number of times a kmer needs to appear

        Returns
        -------
        list
            List of kmers that occurred at least the min number of times
        """
        frequency = self.get_kmer_counts(self.sequence, k)

        return self.get_frequent_kmer(frequency, min_frequency)

    def most_frequent_kmer(self, k: int) -> List[str]:
        """ Get only the most frequently occurring kmers in genome

        Parameters
        ----------
        k : int
            Length of kmers to get counts of

        Returns
        -------
        list
            Most frequently occurring kmers
        """
        frequency = self.get_kmer_counts(self.sequence, k)
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
    def pattern_to_number(cls, pattern: str, nucleotide_map: dict=None) -> int:
        """" Convert nucleotide pattern into base 4 representation where:

        Parameters
        ----------
        pattern: str
            Pattern we want to convert into number
        nucleotide_map : dict, default = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
            What nucleotide we want to map each number to, should be unique and < 4

        Returns
        -------
        int
            Base 4 representation of genome
        """
        if nucleotide_map is None:
            nucleotide_map = cls.__nucleotide_int_map

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
        nucleotide_map : dict, default = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
            What nucleotide we want to map each number to, should be unique and < 4

        Returns
        -------
        str
            nucleotide sequence
        """
        if nucleotide_map is None:
            nucleotide_map = cls.__nucleotide_int_map

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
