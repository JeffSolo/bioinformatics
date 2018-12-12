from itertools import product, chain
from typing import Dict, List, Union
from .distance import hamming_distance
from .dna import DNA


class Motifs(DNA):
    """ Motif finding

    Parameters
    ----------
    dna_strands : list or str
        Set of DNA strands
    separator : str, optional default ' '
        What separate is used between strands if dna_strands is given as a single string
        If dna_strands is not a string, this parameter is not utilized

    Attributes
    ----------
    strands : list
        All given DNA strands
    """

    def __init__(self, dna_strands: Union[List[str], str], separator=' '):
        self.strands = dna_strands.split(separator) if isinstance(dna_strands, str) else dna_strands

    def median_string(self, k: int) -> List[str]:
        """ Find kmers minimizing hamming distance amongst all dna strands

        Parameters
        ----------
        k : int
            length of kmers

        Returns
        -------
        list
            kmers with minimum distance
        """
        min_distance = k * len(self.strands)
        median = []
        for pattern in [''.join(i) for i in product('ACGT', repeat=k)]:
            distance = self.distance_between_pattern_and_strands(pattern)
            if distance < min_distance:
                min_distance = distance
                median = [pattern]
            elif min_distance == distance:
                median.append(pattern)
        return median

    def motif_enumeration(self, k: int, max_distance=0) -> List[str]:
        """ Brute force method for finding motifs that occur in dna strands

        Parameters
        ----------
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
        for strand in self.strands:
            temp = []
            for i in range(len(strand) - k + 1):
                kmer = strand[i: i + k]
                temp.append(set(self._get_neighbors(kmer, max_distance)))
            neighbors_list.append(list(chain.from_iterable(temp)))

        patterns = set(neighbors_list[0])
        for neighbors in neighbors_list[1:]:
            patterns = patterns & set(neighbors)

        return list(set(patterns))

    def distance_between_pattern_and_strands(self, pattern: str) -> int:
        """ Sum the hamming distance between a pattern and each dna strand

        Parameters
        ----------
        pattern : str
            Pattern to check distance against

        Returns
        -------
        int
            Total distance between pattern and strands
        """
        k = len(pattern)
        total_distance = 0

        for strand in self.strands:
            distance = k
            for i in range(len(strand) - k + 1):
                kmer = strand[i: i+k]
                if distance > hamming_distance(pattern, kmer):
                    distance = hamming_distance(pattern, kmer)
            total_distance += distance
        return total_distance

    @staticmethod
    def most_probable_kmer(sequence: str, probability_profile: Dict[str, List[float]], k: int) -> list:
        """ Find the most probable kmers occurring in the sequence

        Parameters
        ----------
        sequence :
            String of nucleotides
        probability_profile : dict
            4xM probability matrix
                keys are the nucleotides
                values are list corresponds to probability of nucleotide occurring in it's index in kmer
        k : int
            length of kmer, should equal to number of columns in *probability_profile*

        Returns
        -------
        list
            Most probable kmers
        """
        most_probable = []
        max_prob = 0
        for i in range(len(sequence) - k + 1):
            prob = 1
            kmer = sequence[i: i + k]
            for index, char in enumerate(kmer):
                prob *= probability_profile[char][index]
            if prob == max_prob:
                most_probable.append(kmer)
            elif prob > max_prob:
                most_probable = [kmer]
                max_prob = prob
        return most_probable
