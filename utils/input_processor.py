__author__ = 'Const'

import re
import numpy as np

def read_score_matrix(path_to_file):
    score_base = [x for x in open(path_to_file)]
    alphabet, score_matrix = score_base[0].split(), map(lambda x: map(lambda x: int(x),  x.split()[1:]), score_base[1:])
    return alphabet, score_matrix


def read_sequence_database(path_to_sequence_database_file):
    sequences = {}
    current_line = ''
    # simple regex for input matching
    input_regex = re.compile(">[^'\n']*")

    for line in open(path_to_sequence_database_file):
        line_trimmed = line[:-1]
        if input_regex.findall(line_trimmed):
            current_line = line_trimmed
            sequences[line_trimmed] = ''
        else:
            sequences[current_line] += line_trimmed
    return sequences


def make_k_partition(sequences, k):
    k_sequences = {}
    for idx, sequence in enumerate(sequences):
        for beg_idx in xrange(len(sequence) - k):
            subseq = sequence[beg_idx:beg_idx + k]
            if subseq in k_sequences:
                k_sequences[subseq].append((idx, beg_idx))
            else:
                k_sequences[subseq] = [(idx, beg_idx)]
    return k_sequences


def read_input_sequence(path_to_input_sequence_file):
    sequences = {}
    current_line = ''
    # simple regex for input matching
    input_regex = re.compile(">[^'\n']*")

    for line in open(path_to_input_sequence_file):
        line_trimmed = line[:-1]
        if input_regex.findall(line_trimmed):
            current_line = line_trimmed
            sequences[line_trimmed] = ''
        else:
            sequences[current_line] += line_trimmed
    return sequences.keys()[0]


def read_nucleotide_sequences(path_to_nucleotide_sequences):
    with open(path_to_nucleotide_sequences, 'r') as seq_file:
        fstrseq = seq_file.readline()[:-1]
        sstrseq = seq_file.readline()[:-1]


def read_nucleotide_sequences(path_to_nucleotide_sequences):
    with open(path_to_nucleotide_sequences, 'r') as seq_file:
        fstrseq = seq_file.readline()[:-1]
        sstrseq = seq_file.readline()[:-1]

def read_nucleotide_raw_sequences(charmap, path_to_nucleotide_sequences):
    fstrseq, sstrseq = read_nucleotide_sequences(path_to_nucleotide_sequences)
    return np.array([charmap[char] for char in fstrseq], dtype=np.int8), np.array([charmap[char] for char in sstrseq], dtype=np.int8)


def read_nucleotide_alphabet_and_matrix(path_to_scoring_matrix):
    with open('score_matrix.txt') as f:
        alphabet = f.readline()[:-1].split(' ')
        alphabet.append('-')
    score_matrix = np.loadtxt('score_matrix.txt', skiprows=1, dtype=int)
    return alphabet, score_matrix