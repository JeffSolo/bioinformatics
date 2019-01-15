from bioinformatics.motifs import Motifs
from bioinformatics.course_helper import parse_genome_file, parse_parameters, print_formatted_output, save_to_file

# motif enumeration
dna, parameters, _ = parse_genome_file('./datasets/week_4/dataset_161_5.txt', has_header=True, join_character=' ')
k, _ = parse_parameters(parameters)
motifs = Motifs(dna)
print_formatted_output(motifs.randomized_motif_search(int(k), iterations=1000), joiner='\n')
