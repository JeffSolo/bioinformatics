from datetime import datetime
import bioinformatics as bio



# submissions are expected to only have space separators
def print_formatted_output(answer):
    if isinstance(answer, list):
        if isinstance(answer[0], str):
            print(' '.join(answer))
        else:
            print(' '.join(str(i) for i in answer))
    else:
        print(answer)


# sometimes multiple parameters are given in the header and footer, separated by a space
def parse_parameters(parameter_string):
    return parameter_string.split(' ')

# pattern count
genome, _, parameter = bio.read_genome('./datasets/week_1/dataset_2_7.txt', has_footer=True)
print_formatted_output(bio.pattern_count(genome, parameter))

# frequent words
genome, _, parameter = bio.read_genome('./datasets/week_1/dataset_2_10.txt', has_footer=True)
print_formatted_output(bio.most_frequent_words(genome, int(parameter)))

# reverse complement
genome, _, _ = bio.read_genome('./datasets/week_1/dataset_3_2.txt')
print_formatted_output(bio.reverse_complement(genome))

# pattern match
genome, parameter, _ = bio.read_genome('./datasets/week_1/dataset_3_5.txt', has_header=True)
print_formatted_output(bio.pattern_match(genome, parameter))

# clump finding
genome, _, parameter = bio.read_genome('./datasets/week_1/dataset_4_5.txt', has_footer=True)
params = parse_parameters(parameter)
print_formatted_output(bio.find_clumps(genome, int(params[0]), int(params[1]), int(params[2])))

# computing frequencies
genome, _, parameter = bio.read_genome('./datasets/week_1/dataset_2994_5.txt', has_footer=True)
print_formatted_output(bio.compute_frequencies(genome, int(parameter)))

# pattern to number
genome, _, _ = bio.read_genome('./datasets/week_1/dataset_3010_2.txt')
print_formatted_output(bio.pattern_to_number(genome))

# number to pattern
number, _, parameter = bio.read_genome('./datasets/week_1/dataset_3010_5.txt', has_footer=True)
print_formatted_output(bio.number_to_pattern(int(number), int(parameter)))

# clump finding for E. Coli
with open('./datasets/E_coli.txt', 'r') as infile:
    gen = infile.read()

print(datetime.now())
with open('./output/e_coli_clumps.txt', 'w') as outfile:
    outfile.write(str(bio.find_clumps(gen, 9, 500, 3)))
print(datetime.now())