#include <iostream>
#include <string>
#include <tuple>
#include <sstream>
#include <ranges>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cassert>

int max(int a, int b, int c){ //position of max of three ints
    if (a >= b && a >= c) return a;
    if (b >= a && b >= c) return b;
    return c;
}

std::tuple<int, std::string, std::string>
AffineAlignment(int match_reward, int mismatch_penalty, int gap_opening_penalty, int gap_extension_penalty,
		std::string s, std::string t){
    int maxScore;
    int maxMiddle;
    int maxLower;
    int maxUpper;

    //s down, t across
    int n = s.length()+1;
    int m = t.length()+1;

    //score matrix
    std::vector<int> row(m, 0);
    std::vector<std::vector<int>> upper(n, row); //gaps in t
    std::vector<std::vector<int>> middle(n, row); //match-mismatch
    std::vector<std::vector<int>> lower(n, row); //gaps in s

    
    std::vector<std::vector<int>> upper_backtrack(n, row);
    std::vector<std::vector<int>> middle_backtrack(n, row);
    std::vector<std::vector<int>> lower_backtrack(n, row);

    //initialize i = 0 and j = 0
    for(int i = 1; i < n; i++){
        upper[i][0] = - gap_opening_penalty - (gap_extension_penalty * i);
        middle[i][0] = - gap_opening_penalty - (gap_extension_penalty * (i-1)); 
        lower[i][0] = - -10000; //not used
        middle_backtrack[i][0] = -1; //down
    }
    for(int j = 1; j < m; j++){
        upper[0][j] = -10000; //not used
        middle[0][j] =  - gap_opening_penalty - (gap_extension_penalty * (j-1)); 
        lower[0][j] = - gap_opening_penalty - (gap_extension_penalty * j); 
        middle_backtrack[0][j] = 1; //side
    }

    //fill in other table values
    std::string str1;
    std::string str2;
    int down; //lower
    int side; //upper
    int diag; //midde
    for(int i = 1; i < n; i++){
        for(int j = 1; j < m; j++){

            //lower
            down = lower[i-1][j] - gap_extension_penalty;
            diag = middle[i-1][j] - gap_opening_penalty;
            lower[i][j] = std::max(down, diag);

            if(lower[i][j] == down){
                lower_backtrack[i][j] = -1; //lower
            }
            else if(lower[i][j] == diag){
                lower_backtrack[i][j] = 0; //middle
            }

            //upper
            side = upper[i][j-1] - gap_extension_penalty;
            diag = middle[i][j-1] - gap_opening_penalty;
            upper[i][j] = std::max(side, diag);

            if(upper[i][j] == side){
                upper_backtrack[i][j] = 1; //upper
            }
            else if(upper[i][j] == diag){
                upper_backtrack[i][j] = 0; //middle
            }

            //middle
            down = lower[i][j];
            side = upper[i][j];
            if(s.substr(i-1, 1).compare(t.substr(j-1, 1)) == 0){
                diag = middle[i-1][j-1] + match_reward;
            }
            else{
                diag = middle[i-1][j-1] - mismatch_penalty;
            }
            
            middle[i][j] = max(down, side, diag);

            //backtrack
            if(middle[i][j] == down){
                middle_backtrack[i][j] = -1; //lower
            }
            else if(middle[i][j] == diag){
                middle_backtrack[i][j] = 0; //middle
            }
            else if(middle[i][j] == side){
                middle_backtrack[i][j] = 1; //upper
            }
        }
    }

    /** change backtrack values
     * backtrack always appends to str then changes matrices
     * doesn't work for the initial row/col
     * change backtrack values to middle such that initial row/col never called
    */
    for(int i = 1; i < n; i++){
        if(upper_backtrack[i][1] == 1){
            upper_backtrack[i][1] = 0;
        }
    }
    for(int j = 1; j < m; j++){
        if(lower_backtrack[1][j] == -1){
            lower_backtrack[1][j] = 0;
        }
    }
    /*for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            std::cout << upper[i][j] << "\t";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            std::cout << upper_backtrack[i][j] << "\t";
        }
        std::cout << "\n";
    }
    std::cout << "\n";*/

    /*maxUpper = upper[n-1][m-1];
    maxMiddle = middle[n-1][m-1];
    maxLower = lower[n-1][m-1];
    maxScore = max(maxUpper, maxMiddle, maxLower);
    
    int start;
    int matrix;
    if(maxScore == maxUpper){
        matrix = 1; //upper
    }
    else if(maxScore == maxMiddle){
        matrix = 0; //middle
    }
    else if(maxScore == maxLower){
        matrix = -1; //lower
    }

    //backtrace to get string
    str1 = ""; //s down
    str2 = ""; //t across
    int i = n-1;
    int j = m-1;
    while(i > 0 || j > 0){
        if(matrix == 1){start = upper_backtrack[n-1][m-1];}
        else if(matrix == 0){start = middle_backtrack[n-1][m-1];}
        else if(matrix == -1){start = lower_backtrack[n-1][m-1];}

        if(start == -1){ //down (lower)
            if(matrix != 0){
                str1.append(s.substr(i-1, 1));
                str2.append("-");
                i -=1;
            }
            else{
                matrix = -1;
            }
        }
        else if(start == 1){ //side (upper)
            if(matrix != 0){
                str1.append("-");
                str2.append(t.substr(j-1, 1));           
                j -=1;
            } 
            else{
                matrix = 1;
            }         
            
        }
        else{ //diag (middle)
            str1.append(s.substr(i-1, 1));
            str2.append(t.substr(j-1, 1));
            i -= 1;
            j -= 1;
        }
        //set start
    }*/

    maxScore = middle[n-1][m-1];

    //backtrace to get string
    str1 = ""; //s down
    str2 = ""; //t across
    int i = n-1;
    int j = m-1;
    int matrix = 0;
    while(i > 0 || j > 0){
        //std::cout << i << " " << j << " m: " << matrix;

        if(matrix == 0){
            //std::cout << " b: " << middle_backtrack[i][j] << " ";
            if(middle_backtrack[i][j] == -1){ //down (lower)
                matrix = -1;
            }
            else if(middle_backtrack[i][j] == 1){ //side (upper)
                matrix = 1;
            }
            else if(middle_backtrack[i][j] == 0){ //diagonal arrow (match/mismatch)
                str1.append(s.substr(i-1, 1));
                str2.append(t.substr(j-1, 1));
                i -= 1;
                j -= 1;
            }
        }
        else if(matrix == -1){ //lower
            //std::cout << " b: " << lower_backtrack[i][j] << " ";
            if(lower_backtrack[i][j] == 1){ //side (upper)
                matrix = 1;
            }
            else if(lower_backtrack[i][j] == 0){ //diagonal arrow (match/mismatch)
                matrix = 0;
            }

            str1.append(s.substr(i-1, 1));
            str2.append("-");
            i -=1;
        }
        else if(matrix == 1){ //upper
            //std::cout << " b: " << upper_backtrack[i][j] << " ";
            if(upper_backtrack[i][j] == -1){ //down (lower)
                matrix = -1;
            }
            else if(upper_backtrack[i][j] == 0){ //diagonal arrow (match/mismatch)
                matrix = 0;
            }

            str1.append("-");
            str2.append(t.substr(j-1, 1));
            j -=1;
        }
        //std::cout << str1 << " " << str2 << " m: " << matrix << "\n";
    }

    std::string rev1 = std::string(str1.rbegin(), str1.rend());
    std::string rev2 = std::string(str2.rbegin(), str2.rend());

    std::tuple<int, std::string, std::string> g(maxScore, rev1, rev2);
    return g;

    //middle: diag of weight score(si, tj) 
    //lower: down of weight -extension
    //upper: side of weight -extension

    //connect middle[i][j] to lower[i-1][j] and upper [i][j-1] w weight -opening
    //lower[i][j] and upper [i][j] to middle[i][j] w weight 0

    /*
    lower[i][j] = max{(lower[i-1][j] - extension), (middle[i-1][j] - opening)}
    upper[i][j] = max{(upper[i-1][j] - extension), (middle[i-1][j] - opening)}
    middle[i][j] = max{(lower[i][j]), (middle[i-1][j-1] + score), (upper[i][j])} //score = match or mismatch
    */
}

int main(){
    //std::tuple<int, std::string, std::string> g = AffineAlignment(1, 5, 2, 1, "CCAT", "GAT");
    std::tuple<int, std::string, std::string> g = AffineAlignment(5, 2, 15, 5, "ACGTA", "ACT");
    std::cout << std::get<0>(g) <<  " " << std::get<1>(g) << " " << std::get<2>(g);

    //ACGTA ACT--
}