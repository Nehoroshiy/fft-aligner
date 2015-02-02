__author__ = 'Const'

from bioalgo.fft_align import take_data, fft_align
from utils.seq import print_aligned_seqences



fseq, sseq, alphabet, score_matrix = take_data(
    seq_input_first_file_path='seq_input_first',
    seq_input_second_file_path='seq_input_second',
    matrix_filepath='score_matrix_dna.txt')


aligned, shift = fft_align((fseq, sseq, alphabet, score_matrix), amount=5)
print 'local alignment + shift step'
print aligned, shift
print '*'*80

print 'aligned_sequences:'
print_aligned_seqences(fseq, sseq, alphabet, aligned, shift)

