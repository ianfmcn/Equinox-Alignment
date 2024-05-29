#!/usr/bin/env python

import argparse
from Bio import SeqIO
import pandas as pd
from locAL import locAL as locAL
from linear_locAL import local_score as linear_locAL
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
                "Default: stdout", type=str, required=True)
    
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

    local = True
    banded = False
    if args.b is not None: #bandwidth align
        if args.d >= 0:
            parser.error("Alignment penalties must be negative")    
        if args.b <= 0:
            parser.error("Bandwidth must be positive")
        banded = True
        local = False
    
    if local: #local align
        if args.d >= 0:
            parser.error("Alignment penalties must be negative")

    reads = pd.DataFrame(columns =  ["id", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "quality"])
    
    # Read and parse sequences of fasta and fastq files
    #f = open("temp.txt", "w")
    for file in args.fq:
        for read_rec in SeqIO.parse(file, "fastq"):
            best_score = float('-inf')
            best_ref_id = ""
            best_loc = ""
            for ref_rec in SeqIO.parse(args.fa, "fasta"):
                if banded:
                    # Perform alignment for each read
                    max_score, rev1, rev2, loc = BandedAL(ref_rec.seq, read_rec.seq, args.m, args.s, args.d, args.b)
                    if max_score > best_score:
                        best_score = max_score
                        best_ref_id = ref_rec.id
                        best_loc = loc
                else:
                    # Perform alignment for each read
                    max_score, max_i, max_j = linear_locAL(ref_rec.seq, read_rec.seq, args.m, args.s, args.d)
                    min_score, min_i, min_j = linear_locAL(ref_rec.seq[:max_i][::-1], read_rec.seq[:max_j][::-1], args.m, args.s, args.d)
                    if max_score > best_score:
                        best_score = max_score
                        best_ref_id = ref_rec.id
                        best_loc = max_i - min_i
            record_score = read_rec.format("fastq").split('\n')[-2]
            reads.loc[-1] = [str(read_rec.id), '*', best_ref_id, int(best_loc), int(best_score), '*', '*', '*', '*', str(read_rec.seq), str(record_score)]
            reads.index = reads.index + 1
            #f.write(str(read_rec.id) + ' * ' + best_ref_id + str(int(best_loc.split('-')[0])) + str(int(best_score)) + ' * '*4 + str(read_rec.seq) + str(record_score) + "\n")
    # f.close()
    #id flag rname(chr) pos mapq cigar rnext pnext tlen seq quality
    #record.id * chr# start_pos max_score * * * * record.seq record.quality
    # rname pos mapq

    reads = reads.sort_values(by=['rname', 'pos'])

    #export DataFrame to text file
    with open(args.out, 'w') as f:
        df_string = reads.to_string(header=True, index=False)
        f.write(df_string)

if __name__ == "__main__":#
    main()

#python3 ./reference_code/parsing.py ./example_files/test_reference.fa ./example_files/test_sequence.fq -m 1 -s -1 -d -1 -o ./example_files/test_local.txt
#python3 ./reference_code/parsing.py ./example_files/test_reference.fa ./example_files/test_sequence.fq -m 1 -s -1 -d -1 -b 5 -o ./example_files/test_banded.txt
