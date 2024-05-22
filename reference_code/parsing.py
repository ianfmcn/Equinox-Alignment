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
    #for testing parsing, use .fa or .fastq
    
    # Output
    parser.add_argument("-o", "--out", help="Write output to file. "\
                "Default: stdout", metavar="FILE", type=str, required=False)
    
    # Alignment Penalties
    parser.add_argument("-m", "--match-reward", \
                help="match reward", \
                type=int, metavar="INT",  required=True)
    parser.add_argument("-s", "--mismatch-penalty", 
                help="mismatch penalty", \
                type=int, metavar="INT",  required=True)
    parser.add_argument("-d", "--indel-penalty", \
                help="indel penalty", \
                type=int, metavar="INT",  required=True)
    
    # Optional Alignment Arguments

    # group gap open and gap extend together
    gap = parser.add_argument_group(title='Affine required args', required=False)
    gap.add_argument("-g", "--gap-open", 
                help="gap open penalty", \
                type=int, metavar="INT",  required=True)
    gap.add_argument("-e", "--gap-ext", \
                help="gap extend penalty", \
                type=int, metavar="INT",  required=True)

    # use either bandwidth align OR affine align
    alt = parser.add_mutually_exclusive_group(title='Alternate alignment args', required=False)
    alt.add_argument("-b", "--bandwidth", \
                help="allowed bandwidth", \
                type=int, metavar="INT",  required=False)
    alt.add_argument(gap)

    '''parser.group("Configuring output format",
        parser:flag "-v --verbose",
        parser:flag "--use-colors",
        parser:option "--encoding"
    )'''
    
    # Other Optional Arguments
    
    other = parser.add_argument_group(title='Other optional args', required=False)
    other.add_argument("-k", "--min-seed-len", \
                help="minimum seed length", \
                type=int, metavar="INT",  required=False)
    other.add_argument("-c", "--max-occ", \
                help="max-occ", \
                type=int, metavar="INT",  required=False)
    other.add_argument("-R", "--RG-line", \
                help="RG-line", \
                type=int, metavar="INT",  required=False)
    other.add_argument("-t", "--cut-output", \
                help="cut-output", \
                type=int, metavar="INT",  required=False)
    other.add_argument("-C", "--comment-FAST", \
                help="comment-FAST", \
                type=int, metavar="INT",  required=False)
    other.add_argument("-v", "--verbose-level", \
                help="verbose-level", \
                type=int, metavar="INT",  required=False)
    
    args = parser.parse_args()

    # Read and parse sequences of fasta and fastq files

    #reference genome, map pair with id (chromosome number) and sequence
    reference_sequences = {}
    for record in SeqIO.parse(args.reference, "fasta"):
        reference_sequences[record.id] = str(record.seq)

    #list of read sequences
    reads = []
    for record in SeqIO.parse(args.reads, "fastq"):
        reads.append(str(record.seq))
    

    
    align = locAL("GATA", "GA", args.match_reward, args.mismatch_penalty, args.indel_penalty, True)
    print(align)


if __name__ == "__main__":
    main()
'''
if -b: run BandedAlignment
if -g and -e: run Affine
else: run locAL

--> other options
'''
    
'''
#from local.py
#establish arguments
    file_index = [idx for idx, s in enumerate(sys.argv) if '.txt' in s][0]
    seq_file = sys.argv[file_index]
    file = open(seq_file)
    content = file.readlines()
    #find two instances of '>'
    #read next line
    index = []
    i = 0
    for line in content:
        l = line.strip()
        if len(l) > 0 and l[0] == '>':
            index.append(i+1)
        i += 1
    s = content[index[0]].strip()
    t = content[index[1]].strip()
    file.close()

    match = int(sys.argv[sys.argv.index("-m")+1])
    mismatch = int(sys.argv[sys.argv.index("-s")+1])
    indel = int(sys.argv[sys.argv.index("-d")+1])

    align_command = [s for s in sys.argv if "-a" in s]
    if len(align_command) == 0:
        align = False
    else:
        align = True
        
    alignment = main(s, t, match, mismatch, indel, align)
    print(alignment)
'''

#vim :set ff=unix :wq 