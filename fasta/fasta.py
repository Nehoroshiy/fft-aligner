import numpy as np
from collections import Counter

def fasta_search(alphabet, score_matrix, sequences, k_sequences, k, m, input_sequence):
    input_subsequences = {}
    for begin_idx in 0, len(input_sequence) - k:
        subseq = input_sequence[begin_idx : begin_idx + k]
        input_subsequences[(subseq, begin_idx)] = []
        if subseq in k_sequences:
            input_subsequences[(subseq, begin_idx)] = list(k_sequences[subseq])

    c = Counter()
    for value in input_subsequences.values():
        for (idx, begin_idx) in value:
            c.update([idx])

    indexes = c.most_common(m)



