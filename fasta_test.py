from utils.input_processor import read_score_matrix, read_sequence_database, make_k_partition, read_input_sequence
from fasta.fasta import fasta_search


alphabet, score_matrix = read_score_matrix("scoring_matrix.txt")
sequences = read_sequence_database("database.txt")
k_sequences = make_k_partition(sequences, 10)
input_sequence = read_input_sequence("input_sequence.txt")
fasta_search(alphabet, score_matrix, sequences, k_sequences, 10, 10, input_sequence)