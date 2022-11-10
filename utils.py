# This module defines functions that are used often within the DropSynth repo

def hamming_dist(a, b):
  return sum(i != j for i, j in zip(a, b))

def rev_comp(seq):
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(comp[base] for base in seq[::-1])