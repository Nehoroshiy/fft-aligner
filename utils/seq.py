from random import randint

def random_DNA_sequence(length):
    nucleotides = ['A', 'T', 'G', 'C']
    return ''.join([nucleotides[randint(0, 3)] for _ in range(length)])