#!/usr/bin/env python
from Bio import SeqIO
import argparse
import numpy as np

def BandedAlignment(s, t, match, mismatch, indel, bandwidth):
    # compute local alignment with bandwidth constraint
    n = len(s) + 1
    m = len(t) + 1
    
    score = np.zeros((2, m), dtype=int)
    
    max_score = 0
    max_i = 0
    max_j = 0
    
    for i in range(1, n):
        i1 = (i - 1) % 2
        i2 = i % 2
        #print(score[i2, :])
        score[i2, :] = 0  # Reset the current row
        
        for j in range(max(1, i - bandwidth), min(m, i + bandwidth)):
            down = score[i1][j] + indel
            side = score[i2][j - 1] + indel
            if s[i - 1] == t[j - 1]:
                diag = score[i1][j - 1] + match
            else:
                diag = score[i1][j - 1] + mismatch

            score[i2, j] = max(down, side, diag, 0)
            if max_score < score[i2, j]:
                max_score = score[i2, j]
                max_i = i
                max_j = j
        #print(i2)
    return max_score, max_i, max_j

