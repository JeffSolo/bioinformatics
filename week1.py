from datetime import datetime
from bioinformatics.genome import Genome
from bioinformatics.course_helper import parse_genome_file, parse_parameters, print_formatted_output, save_to_file

# pattern count
sequence, _, pattern = parse_genome_file('./datasets/week_1/dataset_2_7.txt', has_footer=True)
genome = Genome(sequence)
print_formatted_output(genome.pattern_count(pattern))

# frequent words
sequence, _, k = parse_genome_file('./datasets/week_1/dataset_2_10.txt', has_footer=True)
genome = Genome(sequence)
print_formatted_output(genome.most_frequent_kmer(int(k)))

# reverse complement
sequence, _, _ = parse_genome_file('./datasets/week_1/dataset_3_2.txt')
genome = Genome(sequence)
print_formatted_output(genome.reverse_complement())

# pattern match
sequence, pattern, _ = parse_genome_file('./datasets/week_1/dataset_3_5.txt', has_header=True)
genome = Genome(sequence)
print_formatted_output(genome.pattern_match_index(pattern))

# clump finding
sequence, _, parameters = parse_genome_file('./datasets/week_1/dataset_4_5.txt', has_footer=True)
k, L, t = parse_parameters(parameters)
genome = Genome(sequence)
print_formatted_output(genome.find_clumps(int(k), int(L), int(t)))

# computing frequencies
sequence, _, k = parse_genome_file('./datasets/week_1/dataset_2994_5.txt', has_footer=True)
genome = Genome(sequence)
print_formatted_output(genome.compute_frequencies(int(k)))

# pattern to number
sequence, _, _ = parse_genome_file('./datasets/week_1/dataset_3010_2.txt')
print_formatted_output(Genome.pattern_to_number(sequence))

# number to pattern
number, _, l = parse_genome_file('./datasets/week_1/dataset_3010_5.txt', has_footer=True)
print_formatted_output(Genome.number_to_pattern(int(number), int(l)))

# clump finding for E. Coli
print(datetime.now())
genome = Genome()
genome.read_genome('./datasets/E_coli.txt')
clumps = genome.find_clumps(9, 500, 3)
print(datetime.now())

save_to_file('./output/e_coli_clumps.txt', clumps)
