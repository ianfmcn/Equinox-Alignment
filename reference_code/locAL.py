
#!/usr/bin/env python
from Bio import SeqIO
import argparse
import numpy as np

def locAL(ref, read, match, mismatch, indel):
    #compute local alignment
    #score matrix (ref down, read across)
    n = len(ref)+1
    m = len(read)+1
    score = np.zeros((n, m), dtype=int)

    #initialize i = 0 and j = 0
    for j in range(1, m):
        score[0][j] = 0
    for i in range(1, n):
        score[i][0] = 0

    #fill in other values
    max_score = 0
    max_i = 0
    max_j = 0
    down = 0
    side = 0
    diag = 0
    for i in range(1, n):
        for j in range (1, m):
            down = score[i-1][j] + indel
            side = score[i][j-1] + indel
            if ref[i-1] == read[j-1]:
                diag = score[i-1][j-1] + match
            else:
                diag = score[i-1][j-1] + mismatch

            score[i][j] = max(down, side, diag, 0)
            if max_score < score[i][j]:
                max_score = score[i][j]
                max_i = i
                max_j = j

    #convert to backtrace
    for i in range(n-1, 0, -1):
        for j in range(m-1, 0, -1):
            if score[i][j] == (score[i-1][j] + indel):
                score[i][j] = -1; #down arrow
            elif score[i][j] == (score[i][j-1] + indel):
                score[i][j] = 1; #side arrow
            elif score[i][j] == (score[i-1][j-1] + match):
                score[i][j] = 3; #diagonal arrow
            elif score[i][j] == (score[i-1][j-1] + mismatch):
                score[i][j] = 3; #diagonal arrow
            #if score[i][j] == 0, leave as 0

    #start backtrace at max_i, max_j
    str1 = ''
    str2 = ''
    i = max_i
    j = max_j
    stop = False
    while stop == False:
        if score[i][j] == -1: #down arrow (indel)
            str1 += ref[i-1]
            str2 += "-"
            i -=1
        elif score[i][j] == 1: #side arrow (indel)
            str1 += "-"
            str2 += read[j-1]
            j -=1
        elif score[i][j] == 0: #score of 0 in initial matrix
            stop = True
        else: #diagonal arrow (match/mismatch)
            str1 += ref[i-1]
            str2 += read[j-1]
            i -= 1
            j -= 1
    #print(i, max_i, j, max_j)
    rev1 = str1[::-1]
    rev2 = str2[::-1]

    # We are 1-basing for genome locations
    print(f"length reference: {len(ref)}, length read: {len(read)}, length alignment: {len(rev1)}, end index alignmeent: {max_i}")
    locend = max_i
    locstart =  locend - len(rev1)
    locstr = str(locstart) + "-" + str(locend)
    return max_score, rev1, rev2, locstr

#main function included in this file for testing if the function works properly, delete afterward
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('reference', type=str)
    parser.add_argument('reads', type=str)
    parser.add_argument('-m', '--match', type=float, default=1)
    parser.add_argument('-s', '--mismatch', type=float, default=-1)
    parser.add_argument('-d', '--indel', type=float, default=-1)

    args = parser.parse_args()

    # Read and parse input files
    reference_sequences = {}
    for record in SeqIO.parse(args.reference, "fasta"):
        reference_sequences[record.id] = str(record.seq)

    reads = []
    for record in SeqIO.parse(args.reads, "fastq"):
        reads.append(str(record.seq))

    # Perform alignment for each read
    for read in reads:
        best_score = float('-inf')
        best_aligned_ref = ""
        best_aligned_read = ""
        best_ref_id = ""
        best_loc = ""

        for ref_id, ref_seq in reference_sequences.items():
            max_score, rev1, rev2, loc = locAL(ref_seq, read, args.match, args.mismatch, args.indel)
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

if __name__ == "__main__":
    main()

#py local.py test.txt -m 1 -s -1 -d -1
