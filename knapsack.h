#pragma once 
#include <iostream>
#include <vector>
#include "rim.h"
using namespace std;
 
 
double max(double a, double b) { return (a > b) ? a : b; }
 
/**
Implementation of an Knapsack Dynaminc Programm.
@param im: IM variable 
@param budget: knapsack constraint
@param policy: the policy to examine
@param ret_sol: boolean value that determines if we want the solution or only the score \n 
  | if ret_sol is True then the algorithm has space complexity O(N * budget) \n 
  | else if ret_sol is False then the algorithm has space complexity O(2 * budget)
@return pair<vector<int, double> that includes the solution and its knapsack score 
*/  
pair<vector<int>, double> knapSack(IM *im, int budget, int policy, bool ret_sol, bool rev_res)
{
    pair<vector<int>, double> solution;
    int n = im->mdp->N;
    if(ret_sol){
        float **K;
        K = new float *[n+1];
        for (int i = 0; i < n+1; i++)
        {
            K[i] = new float[budget+1];
        }
        for(int i = 0; i<(n+1); i++){
            for(int w = 0; w<(budget+1); w++){
                if(i == 0 || w == 0){
                    K[i][w] = 0;
                }
                else if(im->mdp->findCost(i-1) <= w){
                    K[i][w] = max(im->score()[policy][i-1] + K[i-1][w-im->mdp->findCost(i-1)], K[i-1][w]);
                }
                else{
                    K[i][w] = K[i-1][w];
                }
            }
        }
        if(rev_res) solution.second = im->all_scores(policy) - K[n][budget];
        else solution.second = K[n][budget];
        int w = budget;
        for(int i=n; i>=1; i--){
            if(K[i][w] != K[i-1][w]){
                if(rev_res) solution.first.push_back(0);
                else solution.first.push_back(1);
                w -= im->mdp->findCost(i-1);
            }
            else{
                if(rev_res) solution.first.push_back(1);
                else solution.first.push_back(0);
            }
        }
        for (int i=0; i<n+1; i++){
            delete[] K[i];
        }
        delete K;
    }
    else{
        int n = im->mdp->N;
        double K[2][budget+1];
        for(int i = 0; i<(n+1); i++){
            for(int w = 0; w<(budget+1); w++){
                if(i == 0 || w == 0){
                    K[i%2][w] = 0;
                }
                else if(im->mdp->findCost(i-1) <= w){
                    K[i%2][w] = max(im->score()[policy][i-1] + K[(i+1)%2][w-im->mdp->findCost(i-1)], 
                                  K[(i+1)%2][w]);
                }
                else{
                    K[i%2][w] = K[(i+1)%2][w];
                }
            }
        }
        if(rev_res) solution.second = im->all_scores(policy) - K[n%2][budget];
        else solution.second = K[n%2][budget];
        
    }
    reverse(solution.first.begin(),solution.first.end());
    return solution;
}
 
