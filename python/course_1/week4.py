from python.bioinformatics import Motifs
from python.bioinformatics import parse_genome_file, parse_parameters, print_formatted_output

# motif enumeration
dna, parameters, _ = parse_genome_file('./datasets/week_4/dataset_161_5.txt', has_header=True, join_character=' ')
k, _ = parse_parameters(parameters)
motifs = Motifs(dna)
print_formatted_output(motifs.randomized_motif_search(int(k), iterations=1000), joiner='\n')

# gibbs sampling
dna, parameters, _ = parse_genome_file('./datasets/week_4/dataset_163_4.txt', has_header=True, join_character=' ')
k, _, n = parse_parameters(parameters)
motifs = Motifs(dna)
print_formatted_output(motifs.gibbs_sampler(int(k), restarts=20, iterations=int(n)), joiner='\n')