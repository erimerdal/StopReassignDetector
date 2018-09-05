import sys
import os
from lib import StopChecker, MasterCollection
import argparse
import shutil


def create_dir(path, force=False):
    if not os.path.exists(os.path.expanduser(path)):
        os.makedirs(path)
    elif force:
        shutil.rmtree(path)
        os.makedirs(path)
    else:
        raise ValueError("Output directory is not empty ! Use --force")
    return path


if __name__ == '__main__':
    # argument parser
    parser = argparse.ArgumentParser(
        description='StopReaDec: a program to detect stop codon reassignments from a list of genomes')

    parser.add_argument('--outdir', '-o', dest="outdir",
                        default=os.path.join(os.getcwd(), "Outdir"), help="Working directory")
    parser.add_argument('--force', '-f', dest="force", action="store_true",
                        help="Erase files in Outdir")
    parser.add_argument('--kmer_length', type=int,
                        default=10,  help="Kmer length for extension")
    parser.add_argument('--threshold', dest="threshold",
                        type=int, default=10,  help="Alignment threshold (?)")
    parser.add_argument('--gap_open_penalty', type=int,
                        default=-2,  help="Gap open penalty for alignment")
    parser.add_argument('--match_score', type=int, default=6,
                        help="Math score between two residues for alignment")
    parser.add_argument('--mismatch_penalty', type=int,
                        default=-4,  help="Mismatch penalty for alignment")
    parser.add_argument('--submat', type=str, choices=('blosum62',),
                        help="Substitutioon matrix to use instead of penalties for the alignment")
    parser.add_argument('seqs', metavar='seqfile1 ... seqfileN',
                        nargs='+', help='List of input files')

    options = parser.parse_args()
    create_dir(options.outdir, options.force)

    infiles = options.seqs

    #files_list = os.listdir(currentWD + "/MasterFiles")
    mcollect = MasterCollection(infiles)
    # for mfile in mcollect.mfilemap:
    #     print('\n--------------------------------------\n')
    #     print(mfile.format_genome())

    schecker = StopChecker(mcollect, options.outdir, kmer_length=options.kmer_length, threshold=options.threshold, gap_open_penalty=options.gap_open_penalty,
                           match_score=options.match_score,    mismatch_penalty=options.mismatch_penalty, submat=options.submat)
