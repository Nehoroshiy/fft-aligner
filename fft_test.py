__author__ = 'Const'

from bioalgo.fft_align import take_data, fft_align
from utils.seq import print_aligned_seqences


# Simple nice test that show how algorithm works
fseq, sseq, alphabet, score_matrix = take_data(
    seq_input_first_file_path='seq_input_first',
    seq_input_second_file_path='seq_input_second',
    matrix_filepath='score_matrix_dna.txt')

"""# More complex test that show memory bottleneck in our algorithm
fseq, sseq, alphabet, score_matrix = take_data(
    seq_input_first_file_path='27_million_first_seq',
    seq_input_second_file_path='big_seq_input_second',
    matrix_filepath='score_matrix_dna.txt')"""


aligned, shift = fft_align((fseq, sseq, alphabet, score_matrix), amount=8, test_flag=True)
print 'Best local alignment + shift step'
print aligned, shift
print '*'*80

print 'best aligned_sequences:'
print_aligned_seqences(fseq, sseq, alphabet, aligned, shift)

