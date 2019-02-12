from collections import Counter
from itertools import product, chain
from random import randint, sample
from typing import Dict, List, Union

from numpy.random import choice

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
        All given DNA strands, should be of the same length
    """

    def __init__(self, dna_strands: Union[List[str], str], separator=' '):
        # TODO check if strands are the same length
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
    def _build_kmer_probabilities(sequence: str, probability_profile: Dict[str, List[float]], k: int) -> list:
        """" Calculate probabilities for each kmer in a sequence based on it's probability profile

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
            List of probabilities - index of probability corresponds to start index of kmer
        """
        probabilities = []
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i: i + k]
            prob = 1
            for index, char in enumerate(kmer):
                prob *= probability_profile[char][index]
            probabilities.append(prob)
        return probabilities

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
        k : int
            length of motif
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
        num_kmers = len(self.strands[0]) - k + 1 # assume all strands are the same length

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

    def randomized_motif_search(self, k: int, iterations=1000):
        """ Randomly select kmers from each strand, create profile, and calculate the best score over many iterations

        Parameters
        ----------
        k : int
            length of motif
        iterations : int optional default 1000
            number of times to run algorithm

        Returns
        -------
        list
            list of strings making up best motif
        """
        max_distance = len(self.strands) * k
        num_possible_kmers = len(self.strands[0]) - k + 1  # assume all strands are the same length, TODO fix or check
        best_overall_motifs = []
        best_overall_distance = max_distance

        for _ in range(iterations):
            kmer_starts = sample(range(num_possible_kmers), len(self.strands))
            best_iter_motifs = [strand[kmer_starts[i]: kmer_starts[i] + k] for i, strand in enumerate(self.strands)]
            best_iter_distance = max_distance + 1
            motifs = best_iter_motifs
            distance = max_distance
            profile = self._make_profile(motifs, pseudocount=True)

            while distance < best_iter_distance:
                best_iter_distance = distance
                best_iter_motifs = motifs

                # use first most probable kmer encountered
                motifs = [self._most_probable_kmer(strand, profile, k)[0] for strand in self.strands]
                profile = self._make_profile(motifs, pseudocount=True)

                distance = hamming_distance(motifs, self._most_probable_strand(profile))

            if best_iter_distance < best_overall_distance:
                best_overall_distance = best_iter_distance
                best_overall_motifs = best_iter_motifs

        return best_overall_motifs

    #TODO consider refactor - too much duplicate code from randomized_motif_search
    def gibbs_sampler(self, k: int, restarts=20, iterations=1000):
        """ Randomly select kmers from each strand, create profile, and calculate the best score,
         only changing 1 kmer between iterations

        Parameters
        ----------
        k : int
            length of motif
        restarts : int optional default 20
            number of times we want to try  with different initial random kmers
        iterations : int optional default 1000
            number of times to run algorithm for each random set of kmers


        Returns
        -------
        list
            list of strings making up best motif
        """
        num_strands = len(self.strands)
        max_distance = num_strands * k
        num_possible_kmers = len(self.strands[0]) - k + 1  # assume all strands are the same length, TODO fix or check
        best_overall_motifs = []
        best_overall_distance = max_distance

        for _ in range(restarts):
            kmer_starts = sample(range(num_possible_kmers), num_strands)

            best_iter_motifs = [strand[kmer_starts[i]: kmer_starts[i] + k] for i, strand in enumerate(self.strands)]
            best_iter_distance = max_distance + 1
            motifs = best_iter_motifs

            for _ in range(iterations):
                exclude_index = randint(0, num_strands - 1)
                motifs.pop(exclude_index)
                deleted_strand = self.strands[exclude_index]
                profile = self._make_profile(motifs, pseudocount=True)
                unweighted_probs = self._build_kmer_probabilities(deleted_strand, profile, k)
                total_prob = sum(unweighted_probs)
                probs = [i / total_prob for i in unweighted_probs]

                random_index = choice(range(num_possible_kmers), p=probs)

                new_random_kmer = deleted_strand[random_index: random_index + k]

                motifs.insert(exclude_index, new_random_kmer)

                distance = hamming_distance(motifs, self._most_probable_strand(profile))

                if distance < best_iter_distance:
                    best_iter_distance = distance
                    best_iter_motifs = motifs

            if best_iter_distance < best_overall_distance:
                best_overall_distance = best_iter_distance
                best_overall_motifs = best_iter_motifs

        return best_overall_motifs
