from random import randint
import numpy as np

def random_DNA_sequence(length):
    nucleotides = ['A', 'T', 'G', 'C']
    return ''.join([nucleotides[randint(0, 3)] for _ in range(length)])


def print_aligned_seqnences(fseq, sseq, alphabet, aligned, shift):
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
    print sresult