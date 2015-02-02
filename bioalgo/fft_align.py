__author__ = 'Const'

import numpy as np
from numpy.fft import fft, ifft
from utils.input_processor import read_nucleotide_alphabet_and_matrix
from smith_waterman import mtx_construct, way_back


def take_data(seq_input_first_filepath, seq_input_second_filepath, matrix_filepath):
    alphabet, score_matrix = read_nucleotide_alphabet_and_matrix(matrix_filepath)

    # This is data preparing part.
    fseq = np.fromfile(seq_input_first_filepath,  dtype=np.int8, count=-1)
    sseq = np.fromfile(seq_input_second_filepath, dtype=np.int8, count=-1)
    N = fseq.size
    M = len(sseq)
    K = len(alphabet) - 1

    for i in xrange(K):
        print ord(alphabet[i])
        fseq[fseq == ord(alphabet[i])] = i
    for i in xrange(K):
        sseq[sseq == ord(alphabet[i])] = i
    print 'fseq:', fseq
    print 'sseq:', sseq
    return fseq, sseq, alphabet, score_matrix



def fft_align((fseq, sseq, alphabet, score_matrix), amount=1):
    N = fseq.size
    M = len(sseq)
    K = len(alphabet) - 1

    fdata = np.zeros(N)
    correlation = np.zeros(N)
    sdata = np.zeros(N)
    for i in xrange(K):
        scoring_row = score_matrix[i]
        for j in xrange(K):
            fdata[fseq == j] = scoring_row[j]
        sdata_view = sdata[:M]
        sdata_view[sseq == i] = 1
        correlation += ifft(np.conj(fft(fdata)) * fft(sdata)).real

    delta_args = np.argsort(correlation)[::-1][:amount]
    print delta_args
    for delta_arg in delta_args:
        if delta_arg != 0:
            delta_args[(delta_args != delta_arg) & (np.abs(delta_args - delta_arg) < M/2)] = 0
    print delta_args

    delta_args = delta_args[delta_args != 0]
    max_score = -100000000

    for delta_arg in delta_args:
        left_bound = min(N - delta_arg - M, 0)
        right_bound = max(N - delta_arg + M, N)

        search_region = fseq[left_bound: right_bound]
        way_mtx = mtx_construct(search_region, sseq, score_matrix, -1)
        if way_mtx[-1].max() > max_score:
            alignment, shift = way_back(search_region, sseq, way_mtx, score_matrix, -1)
    return alignment, shift


    """fresult = ''.join(alphabet[charcode] for charcode in fseq)
    sbuff = np.zeros(N, dtype=np.int8)
    for d in deltas:
        sbuff[:] = -1
        sbuff[N - d:N - d + M] = sseq[:]
        sresult = ''.join(alphabet[charcode] for charcode in sbuff)
        print fresult
        print sresult
        print 'correlation:', correlation[d]
        print 'shift is ', N - d
        print '*'*80"""