#!/usr/bin/env python

import argparse
from Bio import SeqIO
from locAL import main as locAL
from BandedAlignment import main as BandedAL
#from AffineAlignment import main as AffineAL

def main():
    parser = argparse.ArgumentParser(
        prog="equinox",
        description="Command-line script to perform alignment of reads in FASTQ files"
    )

    # Input
    # accept 1 or more input files
    # input files can be .fa or .fastq 
    parser.add_argument("fa", help="FASTA or FASTQ file(s)", type=str, nargs='+')
    
    # Output
    parser.add_argument("-o", "--out", help="Write output to file. "\
                "Default: stdout", metavar="FILE", type=str, required=False)
    
    # Alignment Penalties
    al = parser.add_argument_group()
    al.title = "required arguments for local alignment"
    al.add_argument("-m", help="reward for sequence match", \
                type=int, metavar="INT",  required=True)
    al.add_argument("-s", help="penalty for a mismatch", \
                type=int, metavar="INT",  required=True)
    # in bwa mem, -s == -d (I think)
    al.add_argument("-d", help="penalty for an indel", \
                type=int, metavar="INT",  required=True)
    
    # Optional Alignment Arguments

    # group alternate alignment args
    alt = parser.add_argument_group()
    alt.title = "arguments for alternate alignment algorithms"
    alt.add_argument("-b", help="allowed bandwidth for banded alignment", \
                type=int, metavar="INT",  required=False)
    alt.add_argument("-g",help="gap open penalty for affine alignment", \
                type=int, metavar="INT",  required=False)
    alt.add_argument("-e", help="gap extend penalty; a gap of size k cost '{-g} + {-e}*k'", \
                type=int, metavar="INT",  required=False)
    
    # Other Optional Arguments
    
    other = parser.add_argument_group()
    other.title = "other optional args"
    other.add_argument("-k", help="minimum seed length", \
                type=int, metavar="INT",  required=False)
    other.add_argument("-c", help="max-occ", \
                type=int, metavar="INT",  required=False)
    other.add_argument("-R",help="RG-line", \
                type=int, metavar="INT",  required=False)
    other.add_argument("-t", help="cut-output", \
                type=int, metavar="INT",  required=False)
    other.add_argument("-C", help="comment-FAST", \
                type=int, metavar="INT",  required=False)
    other.add_argument("-v", help="verbose-level", \
                type=int, metavar="INT",  required=False)
    
    args = parser.parse_args()
    
    # make sure there are not more than two files inputted
    if len(args.fa) > 2:
        parser.error("Input either 1 file for single read or 2 files for paired read")

    # for multiple input files, files should be same format
    if len(args.fa) == 2:
        first = args.fa[0][args.fa[0].index("."):]
        second = args.fa[1][args.fa[1].index("."):]
        if first != second:
            parser.error("Input files must be same format")
    
    # use either bandwidth align OR affine align
    # if -b, no -g -e
    if args.b is not None:
        if (args.g is not None) or (args.e is not None):
            parser.error("Cannot use -b and -g/-e together")
    # if -g, must have -e
    if (args.g is not None and args.e is None) or (args.g is None and args.e is not None):
        parser.error("Must use -g/-e together")

    # Read and parse sequences of fasta and fastq files

    #reference genome, map pair with id (chromosome number) and sequence
    reference_sequences = {}
    for record in SeqIO.parse(args.reference, "fasta"):
        reference_sequences[record.id] = str(record.seq)

    #list of read sequences
    reads = []
    for record in SeqIO.parse(args.reads, "fastq"):
        reads.append(str(record.seq))

    #default for now
    align = locAL("GATA", "GA", args.m, args.s, args.d, True)
    print(align)


if __name__ == "__main__":#
    main()

'''
if -b: run BandedAlignment
if -g and -e: run Affine
else: run locAL

--> other options
'''