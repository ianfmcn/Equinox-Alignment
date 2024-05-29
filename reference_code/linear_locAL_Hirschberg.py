import sys
import numpy as np
import math
sys.setrecursionlimit(100000000)

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
    return [max_score, max_i, max_j]

def global_score(s, t, match, mismatch, indel):
    #compute global alignment
    #score matrix (s down, t across)
    n = len(s)+1
    m = len(t)+1
    
    score = np.zeros((2, m), dtype=int)
    
    down = 0
    side = 0
    diag = 0
    i1 = 0
    i2 = 0
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
    return score[i2]

def globAL(s, t, match, mismatch, indel):
    #compute global alignment
    #score matrix (s down, t across)
    n = len(s)+1
    m = len(t)+1
    score = np.zeros((n, m), dtype=int)

    #initialize i = 0 and j = 0
    for j in range(1, m):
        score[0][j] = 0
    for i in range(1, n):
        score[i][0] = 0

    #fill in other values
    down = 0
    side = 0
    diag = 0
    for i in range(1, n):
        for j in range (1, m):
            down = score[i-1][j] + indel
            side = score[i][j-1] + indel
            if s[i-1] == t[j-1]:
                diag = score[i-1][j-1] + match
            else:
                diag = score[i-1][j-1] + mismatch

            score[i][j] = max(down, side, diag)
    max_score = score[n- 1][m-1]

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
    
    #update backtrace for i = 0 and j = 0
    for j in range(1, m):
        score[0][j] = 1 # side
    for i in range(1, n):
        score[i][0] = -1 #down

    str1 = ''
    str2 = ''
    i = n-1
    j = m-1
    while i != 0 or j != 0:
        if score[i][j] == -1: #down arrow (indel)
            str1 += s[i-1]
            str2 += "-"
            i -=1
        elif score[i][j] == 1: #side arrow (indel)
            str1 += "-"
            str2 += t[j-1]
            j -=1
        else: #diagonal arrow (match/mismatch)
            str1 += s[i-1]
            str2 += t[j-1]
            i -= 1
            j -= 1

    rev1 = str1[::-1]
    rev2 = str2[::-1]

    return [str(max_score), str1, str2]

def Hirschberg(s, t, match, mismatch, indel):
    #bc this is called on the substr containing the alignment, we can use global hirschberg to find the alignment sequences
    str1 = ''
    str2 = ''
    score = 0
    if len(s) == 0:
        for i in range (0, len(t)):
            str1 += '-'
            str2 += t[i]
    elif len(t) == 0:
        for i in range (0, len(s)):
            str1 += s[i]
            str2 += '-'
    elif len(s) == 1 and len(t) == 1:
        alignment = globAL(s, t, match, mismatch, indel)
        str1 = alignment[1]
        str2 = alignment[2]
    else:
        s_mid = math.floor(len(s) / 2)
        s_first = s[0:s_mid]
        score_l = global_score(s_first, t, match, mismatch, indel)

        s_second = s[s_mid:]
        s_second_rev = s_second[::-1]
        t_rev = t[::-1]
        score_r = global_score(s_second_rev, t_rev, match, mismatch, indel)

        stop = max(len(score_l), len(score_r))
        score_sum = [0]*stop
        if stop == len(score_l):
            #extend score_r
            to_add = stop - len(score_r)
            for i in range (0, to_add):
                score_r.append(0)
        elif stop == len(score_r):
            #extend score_l
            to_add = stop - len(score_l)
            for i in range (0, to_add):
                score_l.append(0)

        score_r = score_r[::-1]

        y_mid = 0
        y_score = -math.inf
        for i in range (0, stop):
            score_sum[i] = score_l[i] + score_r[i]
            if y_score < score_sum[i]:
                y_mid = i
                y_score = score_sum[i]

        t_first = t[0:y_mid]
        strs1 = Hirschberg(s_first, t_first, match, mismatch, indel)
        t_second = t[y_mid:]
        strs2 = Hirschberg(s_second, t_second, match, mismatch, indel)

        str1 = str1 + strs1[0] + strs2[0]
        str2 = str2 + strs1[1] + strs2[1]

    return [str1, str2]

def main(s, t, match, mismatch, indel, align):
    maxes = local_score(s, t, match, mismatch, indel)
    max_score = maxes[0]
    max_i = maxes[1]
    max_j = maxes[2]

    s_rev = s[:max_i][::-1]
    t_rev = t[:max_j][::-1]
    mins = local_score(s_rev, t_rev, match, mismatch, indel)
    min_i = len(s_rev) - mins[1]
    min_j = len(t_rev) - mins[2]
    #print(max_score)
    print(min_i, max_i, min_j, max_j)

    #local(s, t) == global(s[min_i:max_i], t[min_j:max_j])
    #to help, call Hirschberg on substr of len 5, then piece together
    #for i in range
    ##
    aligned = Hirschberg(s[min_i:max_i], t[min_j:max_j], match, mismatch, indel)
    
    return_values = [str(max_score), str(len(aligned[0]))]
    if align == True:
        return_values.append(aligned[0])
        return_values.append(aligned[1])
        
    return return_values

if __name__ == "__main__":
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
