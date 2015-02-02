import numpy as np

def construct_way_mtx(seq1, seq2, sim_mtx, gap_penalty):
    way_mtx = np.zeros(shape=(len(seq1)+1, len(seq2)+1))

    for i in range(len(seq1)+1):
        way_mtx[i][0] = gap_penalty * i

    for j in range(len(seq2)+1):
        way_mtx[0][j] = gap_penalty * j

    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            match = way_mtx[i-1][j-1] + sim_mtx[seq1[i-1]][seq2[j-1]]
            delete = way_mtx[i-1][j] + gap_penalty
            insert = way_mtx[i][j-1] + gap_penalty

            way_mtx[i][j] = max([match, insert, delete])

    return way_mtx

def find_way(seq1, seq2, way_mtx, sim_mtx, gap_penalty):
    def get_elems(a, b): return way_mtx[a][b], way_mtx[a-1][b-1], way_mtx[a][b-1], way_mtx[a-1][b]

    alignment_seq1, alignment_seq2 = '', ''
    i, j = len(seq1), len(seq2)

    way = []

    while i > 0 and j > 0:
        way.append((i, j))

        score, score_diag, score_left, score_up = get_elems(i, j)

        if score == score_diag + sim_mtx[seq1[i-1]][seq2[j-1]]:
            alignment_seq1 += seq1[i-1]
            alignment_seq2 += seq2[j-1]
            i, j = i-1, j-1
        elif score == score_left + gap_penalty:
            alignment_seq1 += seq1[i-1]
            alignment_seq2 += '-'
            j -= 1
        elif score == score_up + gap_penalty:
            alignment_seq1 += '-'
            alignment_seq2 += seq2[j-1]
            i -= 1

    while i > 0:
        alignment_seq1 += seq1[i-1]
        alignment_seq2 += '-'
        i -= 1

    while j > 0:
        alignment_seq1 += '-'
        alignment_seq2 += seq2[j-1]
        j -= 1

    reverse_seq = lambda seq: ''.join(list(reversed(seq)))

    return reverse_seq(alignment_seq1), reverse_seq(alignment_seq2), way

def sequence_alignment(seq1, seq2, sim_mtx, gap_penalty, ret_way=False):
    way_mtx = construct_way_mtx(seq1, seq2, sim_mtx, gap_penalty)
    alignment_seq1, alignment_seq2, way = find_way(seq1, seq2, way_mtx, sim_mtx, gap_penalty)

    if ret_way:
        return alignment_seq1, alignment_seq2, (way_mtx, way)
    else:
        return alignment_seq1, alignment_seq2