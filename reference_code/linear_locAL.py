#!/usr/bin/env python
from Bio import SeqIO
import argparse
import numpy as np

def local_score(s, t, match, mismatch, indel):
    #compute local alignment
    #score matrix (s down, t across)
    n = len(s)+1
    m = len(t)+1
    
    score = np.zeros((2, m), dtype=int)
    
    max_score = 0
    max_i = 0
    max_j = 0
    down = 0
    side = 0
    diag = 0
    for i in range(1, n):
        for j in range (1, m):
            i1 = (i-1)%2
            i2 = i%2
            
            down = score[i1][j] + indel
            side = score[i2][j-1] + indel
            if s[i-1] == t[j-1]:
                diag = score[i1][j-1] + match
            else:
                diag = score[i1][j-1] + mismatch

            score[i2, j] = max(down, side, diag, 0)
            if max_score < score[i2, j]:
                max_score = score[i2, j]
                max_i = i
                max_j = j
    return max_score, max_i, max_j

def main(s, t, match, mismatch, indel, align):
    max_score, max_i, max_j = local_score(s, t, match, mismatch, indel)

    s_rev = s[:max_i][::-1]
    t_rev = t[:max_j][::-1]
    min_score, min_i, min_j = local_score(s_rev, t_rev, match, mismatch, indel)
    
    print(min_i, max_i, min_j, max_j)