from bioalgo import needleman_wunsch, smith_waterman, hirschberg
from utils import vizualize

similarity_matrix = {
    'A': {'A':  2, 'G': -1, 'C': -1, 'T': -1},
    'G': {'A': -1, 'G':  2, 'C': -1, 'T': -1},
    'C': {'A': -1, 'G': -1, 'C':  2, 'T': -1},
    'T': {'A': -1, 'G': -1, 'C': -1, 'T':  2},
}

#   Needleman Wunsch
seq1, seq2 = 'GGATCGA', 'GAATTCAGTTA'

alignment_seq1, alignment_seq2, (way_mtx, way) = needleman_wunsch.sequence_alignment(
    seq1, seq2, similarity_matrix, gap_penalty=-5, ret_way=True)

vizualize.viz_way2D(seq1, seq2, alignment_seq1, alignment_seq2, way_mtx, way)

#   Smith Waterman
seq1, seq2 = 'AGCACACA', 'ACACACTA'

alignment_seq1, alignment_seq2, (way_mtx, way) = smith_waterman.sequence_alignment(
    seq1, seq2, similarity_matrix, gap_penalty=-1, ret_way=True)

vizualize.viz_way2D(seq1, seq2, alignment_seq1, alignment_seq2, way_mtx, way)

#   Hirschberg
seq1, seq2 = 'AGTACGCA', 'TATGC'

alignment_seq1, alignment_seq2 = hirschberg.sequence_alignment(seq1, seq2, similarity_matrix, gap_penalty=-2)
print alignment_seq1
print alignment_seq2

#Todo vizualize.viz_way2D(seq1, seq2, alignment_seq1, alignment_seq2, way_mtx, way)