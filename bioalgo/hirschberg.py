import numpy as np
from bioalgo import needleman_wunsch

def rev(x):
    return list(reversed(x))

def nw_score(seq1, seq2, sim_mtx, gap_penalty):
    way_mtx = np.zeros(shape=(len(seq1)+1, len(seq2)+1))

    way_mtx[0][0] = 0

    for j in range(1, len(seq2)+1):
        way_mtx[0][j] = way_mtx[0][j-1] + gap_penalty

    for i in range(1, len(seq1)+1):
        way_mtx[i][0] = way_mtx[i-1][0] + gap_penalty
        for j in range(1, len(seq2)+1):
            match = way_mtx[i-1][j-1] + sim_mtx[seq1[i-1]][seq2[j-1]]
            delete = way_mtx[i-1][j] + gap_penalty
            insert = way_mtx[i][j-1] + gap_penalty

            way_mtx[i][j] = max([match, insert, delete])

    return way_mtx[-1]


def find_way(seq1, seq2, sim_mtx, gap_penalty):
    def nil_len(a, b):
        if len(a) == 0: return '-' * len(b), b
        if len(b) == 0: return a, '-' * len(a)

    def partition_y(scorel, scorer):
            return np.argmax(np.array(scorel) + np.array(rev(scorer)))

    if len(seq1) == 0 or len(seq2) == 0:
        alignment_seq1, alignment_seq2 = nil_len(seq1, seq2)
    elif len(seq1) == 1 or len(seq2) == 1:
        alignment_seq1, alignment_seq2 = needleman_wunsch.sequence_alignment(seq1, seq2, sim_mtx, gap_penalty)
    else:
        x_len, x_mid, y_len = len(seq1)-1, len(seq1)/2-1, len(seq2)-1

        score_l = nw_score(seq1[0:x_mid+1], seq2, sim_mtx, gap_penalty)
        score_r = nw_score(rev(seq1[x_mid+1:]), rev(seq2), sim_mtx, gap_penalty)

        y_mid = partition_y(score_l, score_r)

        find_1 = find_way(seq1[0:x_mid+1], seq2[0:y_mid], sim_mtx, gap_penalty)
        find_2 = find_way(seq1[x_mid+1: ], seq2[y_mid:], sim_mtx, gap_penalty)

        alignment_seq1 = find_1[0]+find_2[0]
        alignment_seq2 = find_1[1]+find_2[1]

    return alignment_seq1, alignment_seq2

def sequence_alignment(seq1, seq2, sim_mtx, gap_penalty):
    alignment_seq1, alignment_seq2 = find_way(seq1, seq2, sim_mtx, gap_penalty)

    return alignment_seq1, alignment_seq2