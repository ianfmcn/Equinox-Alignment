from Bio import SeqIO
import argparse
import numpy as np

def BandedAlignment(s, t, match_reward, mismatch_penalty, indel_penalty, band_parameter):
  
    n = len(s) + 1
    m = len(t) + 1
    score = np.zeros((n, m), dtype=int)

    for j in range(1, m):
        score[0][j] = 0
    for i in range(1, n):
        score[i][0] = 0

    down = 0
    side = 0
    diag = 0
    diff_down = 0
    diff_side = 0
    max_score = 0
    max_i = 0
    max_j = 0
  
    for i in range(1, n):
        for j in range(1, m):
            diff_down = abs(i - 1 - j)
            if diff_down < band_parameter:
                down = score[i - 1][j] + indel_penalty
            else:
                down = -10000

            diff_side = abs(i - (j - 1))
            if diff_side < band_parameter:
                side = score[i][j - 1] + indel_penalty
            else: 
                side = -10000
            if s[i - 1] == t[j - 1]:
                diag = score[i - 1][j - 1] + match_reward
            else: 
                diag = score[i - 1][j - 1] + mismatch_penalty
            score[i][j] = max(down, side, diag, 0)
            if max_score < score[i][j]:
                max_score = score[i][j]
                max_i = i
                max_j = j
    
    for i in range(n - 1, 0, -1):
        for j in range(m - 1, 0, -1):
            if score[i][j] == (score[i - 1][j] + indel_penalty):
                score[i][j] = -1  # down arrow
            elif score[i][j] == (score[i][j - 1] + indel_penalty):
                score[i][j] = 1  # side arrow 
            elif score[i][j] == (score[i - 1][j - 1] + match_reward):
                score[i][j] = 3  # diagonal arrow
            elif score[i][j] == (score[i - 1][j - 1] + mismatch_penalty):
                score[i][j] = 3  # diagonal arrow

    for j in range(1, m):
        score[0][j] = 1  # side
    for i in range(1, n):
        score[i][0] = -1  # down

  #start backtrace at max_i, max_j
    str1 = ""  
    str2 = "" 
    i = max_i
    j = max_j

    stop = False
    while stop == False:
        if score[i][j] == -1:  # down arrow (indel)
            str1 += s[i - 1]
            str2 += "-"
            i -= 1
        elif score[i][j] == 1:  # side arrow (indel)
            str1 += "-"
            str2 += t[j - 1]
            j -= 1
        elif score[i][j] == 0: #score of 0 in initial matrix
            stop = True
        else:  # diagonal arrow (match/mismatch)
            str1 += s[i - 1]
            str2 += t[j - 1]
            i -= 1
            j -= 1

    rev1 = str1[::-1]
    rev2 = str2[::-1]
    
    return max_score, rev1, rev2

#main function included in this file for testing if the function works properly, delete afterward
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('reference', type=str)
    parser.add_argument('reads', type=str)
    parser.add_argument('bandwidth', type=int)
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

        for ref_id, ref_seq in reference_sequences.items():
            max_score, rev1, rev2 = BandedAlignment(ref_seq, read, args.match, args.mismatch, args.indel, args.bandwidth)
            if max_score > best_score:
                best_score = max_score
                best_aligned_ref = rev1
                best_aligned_read = rev2
                best_ref_id = ref_id

        print(f"Read: {read}")
        print(f"Score: {best_score}")
        print(f"{best_ref_id})")
        print(f"{best_aligned_ref}")
        print(f"{best_aligned_read}")
        print("\n")

if __name__ == "__main__":
    main()

#test usage example
#python3 ./reference_code/BandedAlignment.py ./example_files/test_reference.fa ./example_files/test_sequence.fq 5 -m 1 -s -1 -d -1 > test_banded.txt
