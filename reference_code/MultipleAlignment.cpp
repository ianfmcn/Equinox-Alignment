#include <iostream>
#include <string>
#include <tuple>
#include <sstream>
#include <ranges>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cassert>

//DID NOT PASS TEST CASE 7

int max(int a, int b, int c){ //position of max of three ints
    if (a >= b && a >= c) return a;
    if (b >= a && b >= c) return b;
    return c;
}

int score(std::string s, std::string t){
    if(s.compare(t) == 0){
        if(s.compare("-") == 0){
            return 0;
        }
        else{
            return 1;
        }
    }
    else{
        return 0;
    }
}

int spScore(std::string v, std::string w, std::string u){
    int spScore = score(v, w) + score(w, u) + score(v, u);
    return spScore;
}

std::tuple<int, std::string, std::string, std::string> MultipleAlignment(std::string s1, std::string s2, std::string s3){

    int maxScore = 0;
    
    int n = s1.length()+1;
    int m = s2.length()+1;
    int o = s3.length()+1;

    //score matrix
    std::vector<int> row(o, 0);
    std::vector<std::vector<int>> row2(m, row);
    std::vector<std::vector<std::vector<int>>> matrix(n, row2);
    std::vector<std::vector<std::vector<int>>> backtrack(n, row2);

    int m1;
    int m2;
    
    std::string vi;
    std::string wj;
    std::string uk;
    for(int i = 1; i < n; i++){
        for(int j = 1; j < m; j++){
            for(int k = 1; k < o; k++){
                vi = s1.substr(i-1, 1);
                wj = s2.substr(j-1, 1);
                uk = s3.substr(k-1, 1);

                m1 = max(matrix[i-1][j][k] + spScore(vi, "-", "-"), matrix[i][j-1][k] + spScore("-", wj, "-"), matrix[i][j][k-1] + spScore("-", "-", uk));
                m2 = max(matrix[i-1][j-1][k] + spScore(vi, wj, "-"), matrix[i-1][j][k-1] + spScore(vi, "-", uk), matrix[i][j-1][k-1] + spScore("-", wj, uk));
                matrix[i][j][k] = max(m1, m2, matrix[i-1][j-1][k-1] + spScore(vi, wj, uk));
                /**matrix[i][j][k] = max
                 * 
                 * matrix[i-1][j][k] + spScore(vi, -, -)
                 * matrix[i][j-1][k] + spScore(-, wj, -)
                 * matrix[i][j][k-1] + spScore(-, -, uk)
                 * 
                 * matrix[i-1][j-1][k] + spScore(vi, wj, -)
                 * matrix[i-1][j][k-1] + spScore(vi, -, uk)
                 * matrix[i][j-1][k-1] + spScore(-, wj, uk)
                 * 
                 * matrix[i-1][j-1][k-1] + spScore(vi, wj, uk)
                */
               
                //set backtrack values
                if(matrix[i][j][k] == m1){
                    if(matrix[i][j][k] == matrix[i-1][j][k] + spScore(vi, "-", "-")){
                        backtrack[i][j][k] = 1;
                    }
                    else if(matrix[i][j][k] == matrix[i][j-1][k] + spScore("-", wj, "-")){
                        backtrack[i][j][k] = 2;
                    }
                    else if(matrix[i][j][k] == matrix[i][j][k-1] + spScore("-", "-", uk)){ //matrix[i][j][k] == matrix[i][j][k-1] + spScore("-", "-", uk)
                        backtrack[i][j][k] = 3;
                    }
                }
                else if(matrix[i][j][k] == m2){
                    if(matrix[i][j][k] == matrix[i-1][j-1][k] + spScore(vi, wj, "-")){
                        backtrack[i][j][k] = 4;
                    }
                    else if(matrix[i][j][k] == matrix[i-1][j][k-1] + spScore(vi, "-", uk)){
                        backtrack[i][j][k] = 5;
                    }
                    else if(matrix[i][j][k] == matrix[i][j-1][k-1] + spScore("-", wj, uk)){ //matrix[i][j][k] == matrix[i][j-1][k-1] + spScore("-", wj, uk)
                        backtrack[i][j][k] = 6;
                    }
                }
                else if(matrix[i][j][k] = matrix[i-1][j-1][k-1] + spScore(vi, wj, uk)){ //matrix[i][j][k] = matrix[i-1][j-1][k-1] + spScore(vi, wj, uk)
                    backtrack[i][j][k] = 7;
                }
            }
        }
    }

    
    for(int j = 0; j < m; j++){
        for(int k = 0; k < o; k++){
            backtrack[0][j][k] = 2; //[i][j-1][k]
        }
    }
    for(int i = 0; i < n; i++){
        for(int k = 0; k < o; k++){
            backtrack[i][0][k] = 1; // [i-1][j][k]
        }
    }
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            backtrack[i][j][0] = 4; //[i-1][j-1][k]
        }
    }
    for(int k = 0; k < o; k++){
        backtrack[0][0][k] = 3; //[i][j][k-1]
    }

    /*for(int i = 1; i < n; i++){
        for(int j = 1; j < m; j++){
            for(int k = 1; k < o; k++){
                std::cout << matrix[i][j][k] << " ";
            }
            std::cout << "\t";
        }
        std::cout << "\n";
    }*/

    //backtrack
    std::string t1 = "";
    std::string t2 = "";
    std::string t3 = "";
    int i = n-1;
    int j = m-1;
    int k = o-1;
    while(i > 0 || j > 0 || k > 0){
        //std::cout << i << " " << j << " " << k << " " << backtrack[i][j][k] << "\n";

        if(backtrack[i][j][k] == 1){
            i--;
            t1.append(s1.substr(i, 1));
            t2.append("-");
            t3.append("-");
        }
        else if(backtrack[i][j][k] == 2){
            j--;
            t1.append("-");
            t2.append(s2.substr(j, 1));
            t3.append("-");
        }
        else if(backtrack[i][j][k] == 3){
            k--;
            t1.append("-");
            t2.append("-");
            t3.append(s3.substr(k, 1));
        }
        else if(backtrack[i][j][k] == 4){
            i--;
            j--;
            t1.append(s1.substr(i, 1));
            t2.append(s2.substr(j, 1));
            t3.append("-");
        }
        else if(backtrack[i][j][k] == 5){
            i--;
            k--;
            t1.append(s1.substr(i, 1));
            t2.append("-");
            t3.append(s3.substr(k, 1));
        }
        else if(backtrack[i][j][k] == 6){
            j--;
            k--;
            t1.append("-");
            t2.append(s2.substr(j, 1));
            t3.append(s3.substr(k, 1));
        }
        else if(backtrack[i][j][k] == 7 || backtrack[i][j][k] == 0){
            i--;
            j--;
            k--; 

            t1.append(s1.substr(i, 1));
            t2.append(s2.substr(j, 1));
            t3.append(s3.substr(k, 1));
            maxScore++;
        }
        else{
            break;
        }
        
        //std::cout << t1 << "\n" << t2 << "\n" << t3 << "\n" << "\n";
    }
    
    std::cout << backtrack[3][4][6] << " " << backtrack[3][4][7] << "\n";
    std::cout << matrix[3][4][7] << "\n";

    std::cout << matrix[3-1][4-1][7-1] << " + " << spScore(s1.substr(2,1), s2.substr(3,1), s3.substr(6,1)) <<"\n";

    std::cout << matrix[3][4-1][7-1] << " + " << spScore("-", s2.substr(3,1), s3.substr(6,1)) <<"\n";

    std::cout << matrix[3][4][7-1] << " + " << spScore("-", "-", s3.substr(6,1)) <<"\n";

    //std::cout << s1.substr(2,1) << " " <<  s2.substr(3,1) << " " << s3.substr(6,1) <<"\n";
    //<< matrix[3][4-1][7-1] << " " << matrix[3][4][7-1] << "\n";

    t1 = std::string(t1.rbegin(), t1.rend());
    t2 = std::string(t2.rbegin(), t2.rend());
    t3 = std::string(t3.rbegin(), t3.rend());

    std::tuple<int, std::string, std::string, std::string> align(maxScore, t1, t2, t3);
    return align;

    //score = alignment of str1, str2, str3
    //for each pos, score ++ iff ch1 == ch2 == ch3
}

int main(){
    std::tuple<int, std::string, std::string, std::string> align = MultipleAlignment("TGTTTAAAAATGTCCGCAACCATTTC", "GATATAAAACAGGGATAACTGCAATGG", "CCTGCTACTTTATGCCGTCTCCATATGCG");
    /*
    11
    T---G--T---T----------T-AAAAA--TGTCCGC-------A-ACCATTTC-----
    ----GA-TA--TAAAACAGGGAT-A----ACTG-C----------A-A---T---GG---
    -CCTG-CT-ACT----------TTA------TG-C---CGTCTCCATA---T-----GCG
    */
    std::cout << std::get<0>(align) << "\n";
    std::cout << std::get<1>(align) << "\n";
    std::cout << std::get<2>(align) << "\n";
    std::cout << std::get<3>(align) << "\n";
}