'''
python locAl.py <seq-files> -m <match> -s <mismatch> -d <indel> -a

seq-files: FASTA-formatted txt file containing two seq
match, mismatch, indel: int scores

output:
score of best local alignment
length of best local alignment
alignment if '-a'

python locAl.py p1seqs.txt -m 1 -s -10 -d -1
'''

'''
file = open("datafile.txt")
sequence = ''
for line in file.readlines():
    l = line.strip()

sys.argv[0] = locAL.py

'''
#!/usr/bin/env python
import sys
import numpy as np

def main(ref, read, match, mismatch, indel):
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
            if s[i-1] == t[j-1]:
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

    return_values = [str(max_score), str(len(rev1))]
    
    return str(max_score), str(len(rev1)), rev1, rev2

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

#py local.py test.txt -m 1 -s -1 -d -1 -a
'''
https://www.geeksforgeeks.org/command-line-arguments-in-python/
https://stackoverflow.com/questions/2170900/get-first-list-index-containing-sub-string
'''
