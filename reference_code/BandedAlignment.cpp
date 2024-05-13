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
BandedAlignment(int match_reward, int mismatch_penalty, int indel_penalty,
        int band_parameter, std::string s, std::string t){

    //CREATE SCORE MATRIX

    int maxScore;

    //s down, t across
    int n = s.length()+1;
    int m = t.length()+1;

    //score matrix
    std::vector<int> row(m, 0);
    std::vector<std::vector<int>> score(n, row);

    //initialize i = 0 and j = 0
    for(int j = 1; j < m; j++){
        score[0][j] = score[0][j-1] - indel_penalty;
    }
    for(int i = 1; i < n; i++){
        score[i][0] = score[i-1][0] - indel_penalty;
    }

    //fill in other table values
    std::string str1;
    std::string str2;
    int down;
    int side;
    int diag;
    int diff_down;
    int diff_side;
    for(int i = 1; i < n; i++){
        for(int j = 1; j < m; j++){
            //only allow edges from valid nodes

            //down
            diff_down = i-1 - j;
            if(diff_down < 0){ diff_down *= -1;} //abs value 
            if(diff_down < band_parameter){
                down = score[i-1][j] - indel_penalty;
            }
            else{
                down = -10000;
            }

            //side
            diff_side = i - j-1;
            if(diff_side < 0){ diff_side *= -1;} //abs value 
            if(diff_side < band_parameter){
                side = score[i][j-1] - indel_penalty;
            }
            else{
                side = -10000;
            }


            if(s.substr(i-1, 1).compare(t.substr(j-1, 1)) == 0){
                diag = score[i-1][j-1] + match_reward;
            }
            else{
                diag = score[i-1][j-1] - mismatch_penalty;
            }
            
            score[i][j] = max(down, side, diag);
        }
    }

    maxScore = score[n-1][m-1];

    //FOLLOW VALID BACKTRACK
    //update score to backtrace values
    for (int i = n-1; i > 0; i--) {
        for (int j = m-1; j > 0; j--) {
            if(score[i][j] == (score[i-1][j] - indel_penalty)){
                score[i][j] = -1; //down arrow
            }
            else if(score[i][j] == (score[i][j-1] - indel_penalty)){
                score[i][j] = 1; //side arrow (-->)
            }
            else if(score[i][j] == (score[i-1][j-1] + match_reward)){
                score[i][j] = 0; //diagonal arrow
            }
            else if(score[i][j] == (score[i-1][j-1] - mismatch_penalty)){
                score[i][j] = 0; //diagonal arrow
            }
        }
    }

    
    //update backtrace for i = 0, j = 0
    for(int j = 1; j < m; j++){
        score[0][j] = 1; //side
    }
    for(int i = 1; i < n; i++){
        score[i][0] = -1; //down
    }

    /*for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            std::cout << score[i][j] << "\t";
        }
        std::cout << "\n";
    }
    std::cout << "\n";*/

    //backtrace to get string
    str1 = ""; //s down
    str2 = ""; //t across
    int i = n-1;
    int j = m-1;
    while(i != 0 || j != 0){
        if(score[i][j] == -1){ //down arrow (indel)
            str1.append(s.substr(i-1, 1));
            str2.append("-");
            i -=1;
        }
        else if(score[i][j] == 1){ //side arrow (indel)
            str1.append("-");
            str2.append(t.substr(j-1, 1));
            j -=1;
        }
        else{ //diagonal arrow (match/mismatch)
            str1.append(s.substr(i-1, 1));
            str2.append(t.substr(j-1, 1));
            i -= 1;
            j -= 1;
        }
    }
    std::string rev1 = std::string(str1.rbegin(), str1.rend());
    std::string rev2 = std::string(str2.rbegin(), str2.rend());

    /*score[i][j] = max(
        score[i-1][j] - indel penalty (down)
        score[i][j-1] - indel penalty (side)
        score[i-1][j-1] + match reward if vi = wj (diag)
        score[i-1][j-1] - mismatch penalty if vi != wj (diag)
    )
    */

    std::tuple<int, std::string, std::string> g(maxScore, rev1, rev2);
    return g;
    //return tuple (max alignment score, t, backtrace string) //strs that achieve the max score
}