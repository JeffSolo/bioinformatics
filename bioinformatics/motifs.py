from collections import Counter
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
                kmer = strand[i: i + k]
                if distance > hamming_distance(pattern, kmer):
                    distance = hamming_distance(pattern, kmer)
            total_distance += distance
        return total_distance

    def distance_between_patterns_and_strands(self, patterns: Union[list, str]) -> int:
        """ Sum the hamming distance between pattern(s) and each dna strand

        Parameters
        ----------
        patterns : list or str
            Patterns to check distance against,
            If single pattern, compare to all strands, otherwise compare pattern to strand - must be same length

        Returns
        -------
        int
            Total distance between pattern and strands
        """
        if isinstance(patterns, str):
            patterns = [patterns]
        k = len(patterns[0])
        total_distance = 0
        index = 0
        #('x', self.strands, patterns)
        for strand in self.strands:
            distance = k
            for i in range(len(strand) - k + 1):
                kmer = strand[i: i + k]
                if distance > hamming_distance(patterns[index], kmer):
                    distance = hamming_distance(patterns[index], kmer)
            total_distance += distance
            if len(patterns) > 0:
                index += 1
        return total_distance

    @staticmethod
    def _most_probable_strand(probability_profile: Dict[str, List[float]]) -> str:
        """ Return most probable string based on probability profile

        Parameters
        ----------
        probability_profile : dict
            Keys are the nucleotides
            Values are lists of floats whose indices correspond to probability of appearing in that location

        Returns
        -------
        Most probable string, if there are multiple, just return one

        """
        nuc = [i for i in probability_profile.keys()]
        probs = list(zip(*probability_profile.values()))

        most_probable = ''
        for i in probs:
            max_prob = max(i)
            max_index = i.index(max_prob)
            most_probable += nuc[max_index]
        return most_probable

    @staticmethod
    def _most_probable_kmer(sequence: str, probability_profile: Dict[str, List[float]], k: int) -> list:
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

    @staticmethod
    def _make_profile(kmers: List[str], pseudocount=False) -> dict:
        """ Create a probability profile for kmers

        Parameters
        ----------
        kmers : list
            kmers we want to profile
        pseudocount : bool
            Whether we want to increment all profile counts by 1 (to avoid some probabilities being 0)

        Returns
        -------
        dict
            Probability profile
        """
        profile = {'A': [], 'C': [], 'G': [], 'T': []}
        transpose = [''.join(i) for i in zip(*kmers)]
        for i in transpose:
            for key in profile.keys():
                if pseudocount:
                    profile[key].append(Counter(i+'ACGT')[key] / (len(kmers) + 4))
                else:
                    profile[key].append(Counter(i)[key] / len(kmers))
        return profile

    def greedy_motif_search(self, k, use_pseudocount=False) -> List[str]:
        """ Find motif in dna strands. WARNING: not well defined which get returned if multiple are equally probable

        Parameters
        ----------
        k : length of motif
        use_pseudocount : bool
            Whether we want to increment all profile counts by 1 (to avoid some probabilities being 0)

        Returns
        -------
        list
            list of strings making up motif
        """
        best_motifs = [strand[:k] for strand in self.strands]
        best_distance = len(self.strands) * k
        profile = {}
        num_kmers = len(self.strands[0]) - k + 1

        for i in range(num_kmers):
            motifs = [self.strands[0][i: i + k]]
            for strand in self.strands[1:]:
                profile = self._make_profile(motifs, pseudocount=use_pseudocount)
                motifs.append(self._most_probable_kmer(strand, profile, k)[0])  # if tied, get first
            distance = self.distance_between_pattern_and_strands(self._most_probable_strand(profile))
            if distance < best_distance:
                best_motifs = motifs
                best_distance = distance
        return best_motifs
