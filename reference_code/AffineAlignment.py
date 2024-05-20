import io
import sys

def AffineAlignment(match_reward, mismatch_penalty, gap_opening_penalty, gap_extension_penalty, s, t):
    maxScore = None
    maxMiddle = None
    maxLower = None
    maxUpper = None

    # s down, t across
    n = len(s)+1
    m = len(t)+1

    # score matrix
    row = [0] * m
    upper = [row] * n # gaps in t
    middle = [row] * n # match-mismatch
    lower = [row] * n # gaps in s
    
    upper_backtrack = [row] * n
    middle_backtrack = [row] * n
    lower_backtrack = [row] * n

    # initialize i = 0 and j = 0
    for i in range(1,n):
        upper[i][0] = -(gap_opening_penalty) - (gap_extension_penalty * i)
        middle[i][0] = -(gap_opening_penalty) - (gap_extension_penalty * (i-1)) 
        lower[i][0] = -(-10000) # not used
        middle_backtrack[i][0] = -1 # down
    for j in range(1,m):
        upper[0][j] = -10000 # not used
        middle[0][j] =  - gap_opening_penalty - (gap_extension_penalty * (j-1))
        lower[0][j] = - gap_opening_penalty - (gap_extension_penalty * j)
        middle_backtrack[0][j] = 1 # side

    # fill in other table values
    str1 = None
    str2 = None
    down = None # lower
    side = None # upper
    diag = None # midde
    for i in range(1,n):
        for j in range(1,m):
            # lower
            down = lower[i-1][j] - gap_extension_penalty
            diag = middle[i-1][j] - gap_opening_penalty
            lower[i][j] = max(down, diag)
            if (lower[i][j] == down):
                lower_backtrack[i][j] = -1 # lower
            elif (lower[i][j] == diag):
                lower_backtrack[i][j] = 0 # middle

            # upper
            side = upper[i][j-1] - gap_extension_penalty
            diag = middle[i][j-1] - gap_opening_penalty
            upper[i][j] = max(side, diag)

            if (upper[i][j] == side):
                upper_backtrack[i][j] = 1 # upper
            elif (upper[i][j] == diag):
                upper_backtrack[i][j] = 0 # middle
            
            # middle
            down = lower[i][j]
            side = upper[i][j]
          
            # May need to change slice bounds
            if (s[i-1:1] == t[j-1:1]):
                diag = middle[i-1][j-1] + match_reward
            else:
                diag = middle[i-1][j-1] - mismatch_penalty
            
            middle[i][j] = max(down, side, diag)

            # backtrack
            if (middle[i][j] == down):
                middle_backtrack[i][j] = -1 # lower
            elif (middle[i][j] == diag):
                middle_backtrack[i][j] = 0 # middle
            elif (middle[i][j] == side):
                middle_backtrack[i][j] = 1 #upper

    # change backtrack values
    # backtrack always appends to str then changes matrices
    # doesn't work for the initial row/col
    # change backtrack values to middle such that initial row/col never called
    for i in range(1,n):
        if (upper_backtrack[i][1] == 1):
            upper_backtrack[i][1] = 0
    for j in range(1,m):
        if (lower_backtrack[1][j] == -1):
            lower_backtrack[1][j] = 0
    # for i in range(1,n):
    #     for j in range(1,m):
    #         cout << upper[i][j] << "\t"
    #     cout << "\n"
    # cout << "\n"

    # for i in range(1,n):
    #     for j in range(1,m):
    #         cout << upper_backtrack[i][j] << "\t"
    #     cout << "\n"
    # cout << "\n"

    # maxUpper = upper[n-1][m-1]
    # maxMiddle = middle[n-1][m-1]
    # maxLower = lower[n-1][m-1]
    # maxScore = max(maxUpper, maxMiddle, maxLower)
    
    # start = None;
    # matrix = None;
    # if (maxScore == maxUpper):
    #     matrix = 1 # upper
    # else if (maxScore == maxMiddle):
    #     matrix = 0 # middle
    # else if (maxScore == maxLower):
    #     matrix = -1 # lower

    # backtrace to get string
    # str1 = "" # s down
    # str2 = "" # t across
    # i = n-1
    # j = m-1
    # while (i > 0 || j > 0):
    #     if(matrix == 1):
    #       start = upper_backtrack[n-1][m-1]
    #     else if(matrix == 0):
    #       start = middle_backtrack[n-1][m-1]
    #     else if(matrix == -1):
    #       start = lower_backtrack[n-1][m-1]

    #     if (start == -1): # down (lower)
    #         if(matrix != 0):
    #             str1.append(s.substr(i-1, 1))
    #             str1 = str1 + s[i-1:1]
    #             str2 = str2 + "-"
    #             i -= 1
    #         else:
    #             matrix = -1
    #     else if(start == 1): # side (upper)
    #         if(matrix != 0):
    #             str1 = str1 + "-"
    #             str2 = str2 + t[i-1:1]
    #             j -=1
    #         else:
    #             matrix = 1     
    #     else: # diag (middle)
    #         str1 = str1 + s[i-1:1]
    #         str2 = str2 + t[i-1:1]
    #         i -= 1
    #         j -= 1
    #     #set start

    maxScore = middle[n-1][m-1]

    # backtrace to get string
    str1 = "" # s down
    str2 = "" # t across
    i = n-1
    j = m-1
    matrix = 0
    while(i > 0 or j > 0):
        # cout << i << " " << j << " m: " << matrix

        if (matrix == 0):
            # std::cout << " b: " << middle_backtrack[i][j] << " ";
            if (middle_backtrack[i][j] == -1): # down (lower)
                matrix = -1
            elif (middle_backtrack[i][j] == 1): #side (upper)
                matrix = 1
            elif (middle_backtrack[i][j] == 0): #diagonal arrow (match/mismatch)
                # May need to change slice bounds
                str1 = str1 + s[i-1:1]
                str2 = str2 + t[j-1:1]
                i -= 1
                j -= 1
        elif (matrix == -1): #lower
            # cout << " b: " << lower_backtrack[i][j] << " ";
            if (lower_backtrack[i][j] == 1): #side (upper)
                matrix = 1
            elif(lower_backtrack[i][j] == 0): #diagonal arrow (match/mismatch)
                matrix = 0
            str1 = str1 + s[i-1:1] # may need to change slice bounds
            str2 = str2 + "-"
            i -= 1
        elif(matrix == 1): #upper
            # std::cout << " b: " << upper_backtrack[i][j] << " "
            if (upper_backtrack[i][j] == -1): #down (lower)
                matrix = -1
            elif(upper_backtrack[i][j] == 0): #diagonal arrow (match/mismatch)
                matrix = 0
            str1 = str1 + "-"
            str2 = str2 + t[j-1:1] # may need to change slice bounds
            j -=1
        # std::cout << str1 << " " << str2 << " m: " << matrix << "\n"

    rev1 = str1[::-1]
    rev2 = str2[::-1]

    g = (maxScore, rev1, rev2)
    return g

    # middle: diag of weight score(si, tj) 
    # lower: down of weight -extension
    # upper: side of weight -extension

    # connect middle[i][j] to lower[i-1][j] and upper [i][j-1] w weight -opening
    # lower[i][j] and upper [i][j] to middle[i][j] w weight 0

    # lower[i][j] = max{(lower[i-1][j] - extension), (middle[i-1][j] - opening)}
    # upper[i][j] = max{(upper[i-1][j] - extension), (middle[i-1][j] - opening)}
    # middle[i][j] = max{(lower[i][j]), (middle[i-1][j-1] + score), (upper[i][j])} //score = match or mismatch

def main():
    # g = AffineAlignment(1, 5, 2, 1, "CCAT", "GAT")
    # sys.stdout.write(str(g[0]) + " " + g[1] + " " + g[2])
    g = AffineAlignment(5, 2, 15, 5, "ACGTA", "ACT")
    sys.stdout.write(str(g[0]) + " " + g[1] + " " + g[2])
    # ACGTA ACT--

if __name__ == "__main__":
    main()
