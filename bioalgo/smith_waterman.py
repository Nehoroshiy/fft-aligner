import numpy as np


def mtx_construct(fseq, sseq, score_matrix, gap_penalty=-1):
    way_mtx = np.zeros([sseq.size, fseq.size], dtype=np.int)
    m = sseq.size
    n = fseq.size

    for i in np.arange(1, m, 1):
        for j in np.arange(1, n, 1):
            match = way_mtx[i - 1, j - 1] + score_matrix[fseq[j - 1], sseq[i - 1]]
            delete = way_mtx[i - 1, j] + gap_penalty
            insert = way_mtx[i, j - 1] + gap_penalty
            way_mtx[i, j] = np.max([match, insert, delete, 0])
    return way_mtx


def way_back(fseq, sseq, way_mtx, score_matrix, gap_penalty=-1):
    m = sseq.size
    n = fseq.size
    counter = n - 1
    aligned = np.zeros([2, n], dtype=np.int)
    i = m - 1
    j = way_mtx[-1].argmax()

    while i > 0 and j > 0:
        score, diag, left, up = way_mtx[i, j], way_mtx[i-1, j-1], way_mtx[i, j-1], way_mtx[i-1, j]
        if score == diag + score_matrix[fseq[j - 1], sseq[i-1]]:
            i -= 1
            j -= 1
            aligned[:, counter] = [1, 1]
        elif score == left + gap_penalty:
            j -= 1
            aligned[:, counter] = [0, 1]
        elif score == up + gap_penalty:
            i -= 1
            aligned[:, counter] = [1, 0]
        counter -= 1
    if i > 0:
        raise Exception('second string must be shorter!')
    if j >= 0:
        aligned = aligned[:, counter + 1:]
    return aligned, j