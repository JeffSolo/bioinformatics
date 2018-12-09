from bioinformatics.motifs import Motifs
from bioinformatics.course_helper import parse_genome_file, parse_parameters, print_formatted_output, save_to_file

# motif enumeration
dna, parameters, _ = parse_genome_file('./datasets/week_3/dataset_156_8.txt', has_header=True, join_character=' ')
k, max_dist = parse_parameters(parameters)
motifs = Motifs(dna)
print_formatted_output(motifs.motif_enumeration(int(k), int(max_dist)))

# distance between pattern and motifs
dna, pattern, _ = parse_genome_file('./datasets/week_3/dataset_5164_1.txt', has_header=True)
motifs = Motifs(dna)
print_formatted_output(motifs.distance_between_pattern_and_strands(pattern))

# median string
dna, k, _ = parse_genome_file('./datasets/week_3/dataset_158_9.txt', has_header=True, join_character=' ')
motifs = Motifs(dna)
print_formatted_output(motifs.median_string(int(k)))
