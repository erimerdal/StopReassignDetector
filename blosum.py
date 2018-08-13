#!/usr/bin/env python
# Usage: python blosum.py blosum62.txt
#        Then, enter input in "row col" format -- e..g, "s f".
import sys

class InvalidPairException(Exception):
  pass

class Matrix:
  def __init__(self, matrix_filename):
    self._load_matrix(matrix_filename)

  

def main():
  matrix = Matrix("blosum62.txt")

if __name__ == '__main__':
  main()
