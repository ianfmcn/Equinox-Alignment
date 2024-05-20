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
    
    print(max_score)
    print(rev1)
    print(rev2)


def main():
    BandedAlignment("CATA", "CAT", 1, -1, -1, 3)
    BandedAlignment("ACGTA", "ACT", 1, -1, -1, 3)

if __name__ == "__main__":
    main()

test = """
#possible command implementation using argparse

#will need to update function so it can parse through a file with the sequence instead of having the sequence inserted directly into the command

def main():
    parser = argparse.ArgumentParser(description="Perform local alignment between two sequences with bandwidth")
    parser.add_argument('s', type=str, help="Reference sequence")
    parser.add_argument('t', type=str, help="Other sequence")
    parser.add_argument('bandwidth', type=int, help="Bandwidth parameter")
    parser.add_argument('-m', '--match', type=float, default=1, help="Match reward")  
    parser.add_argument('-s', '--mismatch', type=float, default=-1, help="Mismatch penalty")
    parser.add_argument('-d', '--indel', type=float, default=-1, help="Insert/delete penalty")

    args = parser.parse_args()

    BandedAlignment(args.s, args.t, args.match, args.mismatch, args.indel, args.bandwidth)

if __name__ == "__main__":
    main()

#usage example
#python3 BandedAlignment.py CATA CAT 3 -m 1 -s -1 -d -1 
"""