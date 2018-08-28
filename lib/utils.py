import os
import re
import numpy as np

REV_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def reverse_sequence(sequence):
    """Return the reversed sequence (5'->3'::3'->5')"""
    return "".join(REV_DICT.get(nuc, nuc) for nuc in sequence[::-1])


# split_len(): Splits the sequence "seq" into "length" amount of seperate sequences.
def split_len(seq, length, position):
    return [seq[i:i+length] for i in range(position, len(seq), length)]


# getSubStrings(): Divides the given sequence "RNA" into 3-mers i.e [0-3,1-4,2-5 etc.]. Position should be 0,1 or 2.
def getSubStrings(RNA, position):
    return [RNA[i:i+3] for i in range(position, len(RNA), 3)]

def list_filter(l1, reflist):
    return [x for y in reflist for x in l1  if y in x]

def _one_hot_encode(one_hot_array, codon):
    one_hot_vector = np.array([])
    for i in range(len(one_hot_array)):
        if one_hot_array[i] == codon:
            one_hot_vector = np.append(one_hot_vector,1)
        else:
            one_hot_vector = np.append(one_hot_vector,0)
    return one_hot_vector

def _create_one_hot_array():
    one_hot_array = []
    one_hot_array.append('AAA')
    one_hot_array.append('AAG')
    one_hot_array.append('AAC')
    one_hot_array.append('AAU')
    one_hot_array.append('AGA')
    one_hot_array.append('AGG')
    one_hot_array.append('AGC')
    one_hot_array.append('AGU')
    one_hot_array.append('ACA')
    one_hot_array.append('ACG')
    one_hot_array.append('ACC')
    one_hot_array.append('ACU')
    one_hot_array.append('AUA')
    one_hot_array.append('AUG')
    one_hot_array.append('AUC')
    one_hot_array.append('AUU')

    one_hot_array.append('GAA')
    one_hot_array.append('GAG')
    one_hot_array.append('GAC')
    one_hot_array.append('GAU')
    one_hot_array.append('GGA')
    one_hot_array.append('GGG')
    one_hot_array.append('GGC')
    one_hot_array.append('GGU')
    one_hot_array.append('GCA')
    one_hot_array.append('GCG')
    one_hot_array.append('GCC')
    one_hot_array.append('GCU')
    one_hot_array.append('GUA')
    one_hot_array.append('GUG')
    one_hot_array.append('GUC')
    one_hot_array.append('GUU')

    one_hot_array.append('CAA')
    one_hot_array.append('CAG')
    one_hot_array.append('CAC')
    one_hot_array.append('CAU')
    one_hot_array.append('CGA')
    one_hot_array.append('CGG')
    one_hot_array.append('CGC')
    one_hot_array.append('CGU')
    one_hot_array.append('CCA')
    one_hot_array.append('CCG')
    one_hot_array.append('CCC')
    one_hot_array.append('CCU')
    one_hot_array.append('CUA')
    one_hot_array.append('CUG')
    one_hot_array.append('CUC')
    one_hot_array.append('CUU')

    one_hot_array.append('UAA')
    one_hot_array.append('UAG')
    one_hot_array.append('UAC')
    one_hot_array.append('UAU')
    one_hot_array.append('UGA')
    one_hot_array.append('UGG')
    one_hot_array.append('UGC')
    one_hot_array.append('UGU')
    one_hot_array.append('UCA')
    one_hot_array.append('UCG')
    one_hot_array.append('UCC')
    one_hot_array.append('UCU')
    one_hot_array.append('UUA')
    one_hot_array.append('UUG')
    one_hot_array.append('UUC')
    one_hot_array.append('UUU')
    return one_hot_array
