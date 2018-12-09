from bioinformatics.genome import Genome
from bioinformatics.course_helper import parse_genome_file, parse_parameters, print_formatted_output, save_to_file

# motif enumeration
motif_str, parameters, _ = parse_genome_file('./datasets/week_3/dataset_156_8.txt', has_header=True, join_character=' ')
k, max_dist = parse_parameters(parameters)
motifs = parse_parameters(motif_str)
print_formatted_output(Genome.motif_enumeration(motifs, int(k), int(max_dist)))

# distance between pattern and motifs
motif_str, pattern, _ = parse_genome_file('./datasets/week_3/dataset_5164_1.txt', has_header=True)
motifs = parse_parameters(motif_str)
print_formatted_output(Genome.distance_between_pattern_and_strands(motifs, pattern))

# median string
motif_str, k, _ = parse_genome_file('./datasets/week_3/dataset_158_9.txt', has_header=True, join_character=' ')
motifs = parse_parameters(motif_str)
print_formatted_output(Genome.median_string(motifs, int(k)))
