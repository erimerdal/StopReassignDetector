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
    one_hot_array.append('AAT')
    one_hot_array.append('AGA')
    one_hot_array.append('AGG')
    one_hot_array.append('AGC')
    one_hot_array.append('AGT')
    one_hot_array.append('ACA')
    one_hot_array.append('ACG')
    one_hot_array.append('ACC')
    one_hot_array.append('ACT')
    one_hot_array.append('ATA')
    one_hot_array.append('ATG')
    one_hot_array.append('ATC')
    one_hot_array.append('ATT')

    one_hot_array.append('GAA')
    one_hot_array.append('GAG')
    one_hot_array.append('GAC')
    one_hot_array.append('GAT')
    one_hot_array.append('GGA')
    one_hot_array.append('GGG')
    one_hot_array.append('GGC')
    one_hot_array.append('GGT')
    one_hot_array.append('GCA')
    one_hot_array.append('GCG')
    one_hot_array.append('GCC')
    one_hot_array.append('GCT')
    one_hot_array.append('GTA')
    one_hot_array.append('GTG')
    one_hot_array.append('GTC')
    one_hot_array.append('GTT')

    one_hot_array.append('CAA')
    one_hot_array.append('CAG')
    one_hot_array.append('CAC')
    one_hot_array.append('CAT')
    one_hot_array.append('CGA')
    one_hot_array.append('CGG')
    one_hot_array.append('CGC')
    one_hot_array.append('CGT')
    one_hot_array.append('CCA')
    one_hot_array.append('CCG')
    one_hot_array.append('CCC')
    one_hot_array.append('CCT')
    one_hot_array.append('CTA')
    one_hot_array.append('CTG')
    one_hot_array.append('CTC')
    one_hot_array.append('CTT')

    one_hot_array.append('TAA')
    one_hot_array.append('TAG')
    one_hot_array.append('TAC')
    one_hot_array.append('TAT')
    one_hot_array.append('TGA')
    one_hot_array.append('TGG')
    one_hot_array.append('TGC')
    one_hot_array.append('TGT')
    one_hot_array.append('TCA')
    one_hot_array.append('TCG')
    one_hot_array.append('TCC')
    one_hot_array.append('TCT')
    one_hot_array.append('TTA')
    one_hot_array.append('TTG')
    one_hot_array.append('TTC')
    one_hot_array.append('TTT')
    return one_hot_array
