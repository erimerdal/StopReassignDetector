import os
import re

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
