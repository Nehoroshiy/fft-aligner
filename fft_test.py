__author__ = 'Const'

from bioalgo.fft_align import take_data, fft_align, sequencer_test
from bioalgo.smith_waterman import mtx_construct, construct_way_mtx, way_back
import numpy as np
import cProfile

fseq, sseq, alphabet, score_matrix = take_data(
    seq_input_first_filepath='seq_input_first',
    seq_input_second_filepath='seq_input_second',
    matrix_filepath='score_matrix.txt')



sequencer_test((fseq, sseq, alphabet, score_matrix), amount=5)





"""
#mtx = mtx_construct(fseq, sseq, score_matrix)

#aligned, shift = way_back(fseq, sseq, mtx, score_matrix)
cProfile.run("fft_align((fseq, sseq, alphabet, score_matrix), amount=5)", sort=True)
aligned, shift = fft_align((fseq, sseq, alphabet, score_matrix), amount=5)

print aligned, shift


N = fseq.size
M = sseq.size
NMOD = N + aligned.shape[1] - aligned[1].sum()
MMOD = M + aligned.shape[1] - aligned[0].sum()
fslen = aligned[0].sum()
fbuf = np.zeros(NMOD, dtype=np.int8)
sbuf = np.zeros(NMOD, dtype=np.int8)
fbuf[:] = -1
sbuf[:] = -1
sbuf[:shift + 1] = -1
sbuf[shift] = sseq[0]
fbuf[:shift + 1] = fseq[:shift + 1]
fb_view = fbuf[shift + 1: shift + MMOD]
sb_view = sbuf[shift + 1: shift + MMOD]
sb_view[:] = -1
fb_view[:] = -1
sb_view[aligned[0] != 0] = sseq[1:]
fb_view[aligned[1, :] != 0] = fseq[shift + 1: shift + fslen]
fresult = ''.join(alphabet[charcode] for charcode in fbuf)
sresult = ''.join(alphabet[charcode] for charcode in sbuf)
print fresult
print sresult"""

"""N = fseq.size
fresult = ''.join(alphabet[charcode] for charcode in fseq)
fbuff = np.zeros(N + aligned.shape[0], dtype=np.int8)
sbuff = np.zeros(N, dtype=np.int8)
fbuff[:shift] = fseq[:shift]
sbuff[:] = -1
sbuff[shift] = sseq[0]

fdelta = 0
sdelta = 0
for (i, [f, s]) in enumerate(aligned[1:]):
    if aligned[0, i - 1] == f - 1 and aligned[1, i - 1] == s - 1:
        sbuff[shift + i] = sseq[i - sdelta]
        fbuff[shift + i] = fseq[i - fdelta]
    elif aligned[0, i - 1] == f - 1 and aligned[1, i - 1] == s:
        fbuff[shift + i] = ord('-')
        sbuff[shift + i] = sseq[i - sdelta]
    elif



fresult = ''.join(alphabet[charcode] for charcode in fseq)
print fresult
print sresult"""


#print mtx[-1].argmax()
