class DNA:
    """ DNA information

    Attributes
    ----------
    nucleobases : dict
        Key is first letter of nucleobase, value is the full name
    dna_complement (class attribute) : dict
        Complement of each nucleobase, i.e. {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    _nucleotide_int_map (class attribute): dict
        Map of nucleotides to ints in alphabetical order i.e. {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    """
    dna_complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }

    nucleobases = {
        'A': 'Adenine',
        'C': 'Cytosine',
        'G': 'Guanine',
        'T': 'Thymine'
    }

    _nucleotide_int_map = {
        'A': 0,
        'C': 1,
        'G': 2,
        'T': 3
    }

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

    @classmethod
    def get_reverse_complement(cls, pattern: str) -> str:
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