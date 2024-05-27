#!/usr/bin/env python

import argparse
from Bio import SeqIO
from locAL import locAL as locAL
from BandedAlignment import BandedAlignment as BandedAL

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
    al.add_argument("-d", help="penalty for an indel. Required for local and banded alignment, but not for affine.", \
                type=int, metavar="INT",  required=False)
    
    # Optional Alignment Arguments
    alt = parser.add_argument_group()
    alt.title = "arguments for alternate alignment algorithms"
    alt.add_argument("-b", help="allowed bandwidth for banded alignment", \
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
    
    # penalties must be negative (-s -d -g -e)
    if args.s >= 0:
        parser.error("Alignment penalties must be negative")

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
    s = "GATA"
    t = "GA"

    # decide which alignment to use
    if args.b is not None:
        if args.d is None:
            parser.error("the following arguments are required: -d")    
        if args.d >= 0:
            parser.error("Alignment penalties must be negative")    
        if args.b <= 0:
            parser.error("Bandwidth must be positive")
        align = BandedAL(s, t, args.m, args.s, args.d, args.b)
    else:
        if args.d is None:
            parser.error("the following arguments are required: -d")
        if args.d >= 0:
            parser.error("Alignment penalties must be negative")
        align = locAL(s, t, args.m, args.s, args.d, True)

    print(align)


if __name__ == "__main__":#
    main()

'''
if -b: run BandedAlignment
if -g and -e: run Affine
else: run locAL

--> other options
'''

#s ref t genome