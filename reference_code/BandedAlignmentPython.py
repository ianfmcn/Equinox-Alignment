import sys
import numpy as np

def BandedAlignment(match_reward, mismatch_penalty, indel_penalty, band_parameter, s, t):
    maxScore = 0
  
    n = len(s) + 1
    m = len(t) + 1
    score = np.zeros((n, m), dtype=int)

    for j in range(1, m):
        score[0][j] = score[0][j - 1] - indel_penalty

    for i in range(1, n):
        score[i][0] = score[i - 1][0] - indel_penalty

    str1 = ""
    str2 = ""
    down = 0
    side = 0
    diag = 0
    diff_down = 0
    diff_side = 0
  
    for i in range(1, n):
        for j in range(1, m):
            diff_down = abs(i - 1 - j)
            if diff_down < band_parameter:
                down = score[i - 1][j] - indel_penalty
            else:
                down = -10000

            diff_side = abs(i - (j - 1))
            if diff_side < band_parameter:
                side = score[i][j - 1] - indel_penalty
            else: 
                side = -10000
            if s[i - 1] == t[j - 1]:
                diag = score[i - 1][j - 1] + match_reward
            else: 
                diag = score[i - 1][j - 1] + mismatch_penalty
            score[i][j] = max(down, side, diag)
    
  maxScore = score[n - 1][m - 1]
    for i in range(n - 1, 0, -1):
        for j in range(m - 1, 0, -1):
            if score[i][j] == (score[i - 1][j] - indel_penalty):
                score[i][j] = -1  # down arrow
            elif score[i][j] == (score[i][j - 1] - indel_penalty):
                score[i][j] = 1  # side arrow (-->)
            elif score[i][j] == (score[i - 1][j - 1] + match_reward):
                score[i][j] = 0  # diagonal arrow
            elif score[i][j] == (score[i - 1][j - 1] - mismatch_penalty):
                score[i][j] = 0  # diagonal arrow

    for j in range(1, m):
        score[0][j] = 1  # side
    for i in range(1, n):
        score[i][0] = -1  # down

    str1 = ""  # s down
    str2 = ""  # t across
    i = n - 1
    j = m - 1
    while i != 0 or j != 0:
        if score[i][j] == -1:  # down arrow (indel)
            str1 += s[i - 1]
            str2 += "-"
            i -= 1
        elif score[i][j] == 1:  # side arrow (indel)
            str1 += "-"
            str2 += t[j - 1]
            j -= 1
        else:  # diagonal arrow (match/mismatch)
            str1 += s[i - 1]
            str2 += t[j - 1]
            i -= 1
            j -= 1

    rev1 = str1[::-1]
    rev2 = str2[::-1]

    return maxScore, rev1, rev2
  
