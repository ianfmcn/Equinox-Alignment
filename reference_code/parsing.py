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
    # accept 1 reference file
    parser.add_argument("fa", help="FASTA file of reference genome", type=str)


    # accept 1 or more read files
    parser.add_argument("fq", help="FASTQ file(s) of reads to align", type=str, nargs='+')
    
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
    if len(args.fq) > 2:
        parser.error("Input either 1 file for single read or 2 files for paired read")
    
    # penalties must be negative (-s -d -g -e)
    if args.s >= 0:
        parser.error("Alignment penalties must be negative")

    # Read and parse sequences of fasta and fastq files

    reference_sequences = {}
    for record in SeqIO.parse(args.fa, "fasta"):
        reference_sequences[record.id] = str(record.seq)

    reads = []
    # only single end for now
    for record in SeqIO.parse(args.fq, "fastq"):
        reads.append(str(record.seq))

    # decide which alignment to use
    if args.b is not None:
        if args.d is None:
            parser.error("the following arguments are required: -d")    
        if args.d >= 0:
            parser.error("Alignment penalties must be negative")    
        if args.b <= 0:
            parser.error("Bandwidth must be positive")
        # Perform alignment for each read
        for read in reads:
            best_score = float('-inf')
            best_aligned_ref = ""
            best_aligned_read = ""
            best_ref_id = ""
            best_loc = ""
            for ref_id, ref_seq in reference_sequences.items():
                max_score, rev1, rev2, loc = BandedAL(ref_seq, read, args.m, args.s, args.d, args.b)
                if max_score > best_score:
                    best_score = max_score
                    best_aligned_ref = rev1
                    best_aligned_read = rev2
                    best_ref_id = ref_id
                    best_loc = loc
            print(f"Read: {read}")
            print(f"Score: {best_score}")
            print(f"{best_ref_id})", best_loc)
            print(f"{best_aligned_ref}")
            print(f"{best_aligned_read}")
            print("\n")
    else:
        if args.d is None:
            parser.error("the following arguments are required: -d")
        if args.d >= 0:
            parser.error("Alignment penalties must be negative")
        # Perform alignment for each read
        for read in reads:
            best_score = float('-inf')
            best_aligned_ref = ""
            best_aligned_read = ""
            best_ref_id = ""
            best_loc = ""
            for ref_id, ref_seq in reference_sequences.items():
                max_score, rev1, rev2, loc = locAL(ref_seq, read, args.m, args.s, args.d)
                if max_score > best_score:
                    best_score = max_score
                    best_aligned_ref = rev1
                    best_aligned_read = rev2
                    best_ref_id = ref_id
                    best_loc = loc
            print(f"Read: {read}")
            print(f"Score: {best_score}")
            print(f"{best_ref_id})", best_loc)
            print(f"{best_aligned_ref}")
            print(f"{best_aligned_read}")
            print("\n")

if __name__ == "__main__":#
    main()