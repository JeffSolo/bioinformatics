from os.path import exists
from typing import Tuple
from errno import EEXIST

def parse_genome_file(file_path: str, has_header=False, has_footer=False, join_character='') -> Tuple[str, str, str]:
    """ Read in genome file, and keep header/footer separate if necessary.

    Parameters
    ----------
    file_path
        File we want to read genome from
    has_header
        Whether or not file has a header, assumes 1 line
    has_footer
        Whether or not file has a footer, assumes 1 line
    join_character
        What character to separate each line with when turning into a single string
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

        return join_character.join(body), header, footer


# submissions are expected to only have space separators
def print_formatted_output(answer):
    if isinstance(answer, list):
        if isinstance(answer[0], str):
            print(' '.join(answer))
        else:
            print(' '.join(str(i) for i in answer))
    else:
        print(answer)


def save_to_file(file_path: str, output: str, overwrite=False):
    """ Save data to file in 'output' directory

    Parameters
    ----------
    file_path : str
        Where and what to name the file
    output : str
        Information we want to save in file
    overwrite : bool, default False
    """
    if not overwrite and exists(file_path):
        raise IOError(EEXIST, f'File {file_path} already exists')

    with open(f'{file_path}', 'w') as outfile:
        outfile.write(output)


# sometimes multiple parameters are given in the header and footer, separated by a space
def parse_parameters(parameter_string):
    return parameter_string.split(' ')
