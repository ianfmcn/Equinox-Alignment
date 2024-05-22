#!/usr/bin/env python

import argparse
from Bio import SeqIO
from locAL import main as locAL
from BandedAlignment import main as BandedAL
#from AffineAlignment import main as AffineAL

def main():
    parser = argparse.ArgumentParser(
        prog="equinox",
        description="Command-line script to perform alignment of FASTQ files"
    )

    # Input
    parser.add_argument("txt", help="Txt file", type=str)
    #for testing parsing, use .fa or .fastq for equinox
    
    # Output
    parser.add_argument("-o", "--out", help="Write output to file. "\
                "Default: stdout", metavar="FILE", type=str, required=False)
    
    # Alignment Penalties
    al = parser.add_argument_group()
    al.title = "Required alignment args"
    al.add_argument("-m", help="match reward", \
                type=int, metavar="INT",  required=True)
    al.add_argument("-s", help="mismatch penalty", \
                type=int, metavar="INT",  required=True)
    al.add_argument("-d", help="indel penalty", \
                type=int, metavar="INT",  required=True)
    
    # Optional Alignment Arguments

    # group alternate alignment args
    alt = parser.add_argument_group()
    alt.title = "Alternate alignment algorithm arguments"
    alt.add_argument("-b", help="allowed bandwidth", \
                type=int, metavar="INT",  required=False)
    alt.add_argument("-g",help="gap open penalty", \
                type=int, metavar="INT",  required=False)
    alt.add_argument("-e", help="gap extend penalty", \
                type=int, metavar="INT",  required=False)
    
    # Other Optional Arguments
    
    other = parser.add_argument_group()
    other.title = "Other optional args"
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
    
    # use either bandwidth align OR affine align
    # -b is mutually exclusive with -g/-e
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