__author__ = 'Const'

import numpy as np
np.set_printoptions(suppress=True, linewidth=200)
from numpy.fft import fft, ifft, rfft, irfft
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
        fseq[fseq == ord(alphabet[i])] = i
    for i in xrange(K):
        sseq[sseq == ord(alphabet[i])] = i
    return fseq, sseq, alphabet, score_matrix


def fft_align((fseq, sseq, alphabet, score_matrix), amount=1):
    N = fseq.size
    M = len(sseq)
    K = len(alphabet) - 1

    N2 = power_2_bound(N)
    fdata = np.zeros(N2, np.float)
    correlation = np.zeros(N2, np.float)
    sdata = np.zeros(N2, np.float)
    for i in xrange(K):
        scoring_row = score_matrix[i]
        for j in xrange(K):
            fdata[fseq == j] = scoring_row[j]
        sdata_view = sdata[:M]
        sdata_view[sseq == i] = 1
        correlation += ifft_fast(np.conj(fft(fdata)) * fft_fast(sdata)).real

    delta_args = np.argsort(correlation)[::-1][:amount]
    print 'initial best shifts:'
    print N2 - delta_args
    for delta_arg in delta_args:
        if delta_arg != 0:
            delta_args[(delta_args != delta_arg) & (np.abs(delta_args - delta_arg) < M)] = 0


    delta_args = delta_args[delta_args != 0]
    delta_args = N2 - delta_args
    print 'filtered best shifts:'
    print delta_args
    max_score = -100000000

    for delta_arg in delta_args:
        left_bound = max(delta_arg - M, 0)
        right_bound = min(delta_arg + M, N)

        search_region = fseq[left_bound: right_bound]
        way_mtx = mtx_construct(search_region, sseq, score_matrix, -1)
        way_max = way_mtx[-1].max()
        if way_max > max_score:
            alignment, shift = way_back(search_region, sseq, way_mtx, score_matrix, -1)
            max_score = way_max
            max_shift = left_bound + shift
    return alignment, max_shift


# returns smallest power of 2 that is greater or equal than v
def power_2_bound(v):
    v -= 1
    v |= v >> 1
    v |= v >> 2
    v |= v >> 4
    v |= v >> 8
    v |= v >> 16
    return v + 1


def fft_fast(a):
    N = a.shape[0]

    if N & (N - 1) != 0:
        raise ValueError("size of x must be a power of 2")

    N_min = min(N, 32)

    n = np.arange(N_min)
    k = n[:, None]
    M = np.exp(-2j * np.pi * n * k / N_min)
    X = np.dot(M, a.reshape([N_min, -1]))

    n = np.arange(N_min)
    while X.shape[0] < N:
        X_even = X[:, :X.shape[1] / 2]
        X_odd = X[:, X.shape[1] / 2:]
        factor = np.exp(-1j * np.pi * np.arange(X.shape[0])
                        / X.shape[0])[:, None]
        X = np.vstack([X_even + factor * X_odd,
                       X_even - factor * X_odd])
    return X.ravel()

def ifft_fast(y):
    N = y.shape[0]

    if N & (N - 1) != 0:
        raise ValueError("size of x must be a power of 2")

    N_min = min(N, 32)

    n = np.arange(N_min)
    k = n[:, None]
    #T = np.exp(-2j * np.pi * n * k / N_min)
    #M = np.conj(T / N_min)
    M = np.exp((2j * np.pi * n) * k / N_min) / N_min
    X = np.dot(M, y.reshape([N_min, -1]))

    n = np.arange(N_min)
    while X.shape[0] < N:
        X_even = X[:, :X.shape[1] / 2]
        X_odd = X[:, X.shape[1] / 2:]
        factor = np.exp(1j * np.pi * np.arange(X.shape[0])
                        / X.shape[0])[:, None]
        X = np.vstack([((X_even + factor * X_odd) / 2),
                       (X_even - factor * X_odd) / 2])
    return X.ravel()


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