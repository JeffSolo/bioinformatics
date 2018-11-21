from bioinformatics.genome import Genome
from bioinformatics.course_helper import parse_genome_file, parse_parameters, print_formatted_output, save_to_file
from bioinformatics.distance import hamming_distance


# minimum skew
sequence, _, _ = parse_genome_file('./datasets/week_2/dataset_7_6.txt')
genome = Genome(sequence)
print_formatted_output(genome.minimum_skew())

# hamming distance
_, str_1, str_2 = parse_genome_file('./datasets/week_2/dataset_9_3.txt', has_header=True, has_footer=True)
print_formatted_output(hamming_distance(str_1, str_2))

# approximate pattern match
sequence, pattern, max_dist = parse_genome_file('./datasets/week_2/dataset_9_4.txt', has_header=True, has_footer=True)
genome = Genome(sequence)
print_formatted_output(genome.pattern_match_index(pattern, int(max_dist)))

# approximate pattern count
sequence, pattern, max_dist = parse_genome_file('./datasets/week_2/dataset_9_6.txt', has_header=True, has_footer=True)
genome = Genome(sequence)
print_formatted_output(genome.pattern_count(pattern, int(max_dist)))
