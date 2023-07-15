#pragma once 
#include <iostream>
#include <vector>
#include "rim.h"
using namespace std;
 
/**
* @brief Max between doubles a and b
* @param a 
* @param b 
* @return max(a,b) 
*/
double max(double a, double b) { return (a > b) ? a : b; }
 

/**
* @brief Implementation of an Knapsack Dynaminc Programm.
* @param im: IM variable 
* @param budget: knapsack constraint
* @param policy: the policy to examine
* @param ret_sol: boolean value that determines if we want the solution or only the score \n 
  | if ret_sol is True then the algorithm has space complexity O(N * budget) \n 
  | else if ret_sol is False then the algorithm has space complexity O(2 * budget)
* @return pair<vector<int, double> that includes the solution and its knapsack score 
*/  
void knapSack(IM *im, int budget, int policy, bool ret_sol, bool rev_res, pair<vector<int>, double>* solution)
{
    //pair<vector<int>, double> solution;
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
        if(rev_res) solution->second = im->all_scores(policy) - K[n][budget];
        else solution->second = K[n][budget];
        int w = budget;
        for(int i=n; i>=1; i--){
            if(K[i][w] != K[i-1][w]){
                if(rev_res) solution->first.push_back(0);
                else solution->first.push_back(1);
                w -= im->mdp->findCost(i-1);
            }
            else{
                if(rev_res) solution->first.push_back(1);
                else solution->first.push_back(0);
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
        if(rev_res) solution->second = im->all_scores(policy) - K[n%2][budget];
        else solution->second = K[n%2][budget];
        
    }
    reverse(solution->first.begin(),solution->first.end());
    //return solution;
}
 

/**
* @brief Sum of vectors.
* @param v1 pointer to float vector
* @param v2 pointer to float vector
* @return v1+v2
*/
vector<float> SUM_vectors(vector<float>* v1, vector<float>* v2){
    vector<float> sum_vec;
    for(int i=0; i<v1->size(); i++){
        sum_vec.push_back(v1->at(i)+v2->at(i));
    }
    return sum_vec;
}


/**
* @brief Min value of a vector.
* @param v1 pointer to float vector
* @return Min of vector v1.
*/
float min(vector<float>* v1){
    float min_val=v1->at(0);
    for(auto el: *v1){
        if(el<min_val) min_val=el;
    }
    return min_val;
}



/**
* @brief Implementation of an Min-Max Knapsack Heuristic Dynaminc Programm.
* @param im: IM variable 
* @param budget: knapsack constraint
* @return pair<vector<int, double> that includes the solution and its knapsack score 
*/  
pair<vector<int>, double> MinMaxknapSack(IM *im, int budget)
{
    pair<vector<int>, double> solution;
    int n = im->mdp->N;
    if(true){
        vector<vector<float>> **K;
        K = new vector<vector<float>> *[n+1];
        for (int i = 0; i < n+1; i++)
        {
            K[i] = new vector<vector<float>>[budget+1];
        }
        float **K2;
        K2 = new float *[n+1];
        for (int i = 0; i < n+1; i++)
        {
            K2[i] = new float[budget+1];
        }
        vector<float> current_v;
        for(int policy = 0; policy<im->mdp->P; policy++){
            current_v.push_back(0);
        }
        //printVector(current_v);
        for(int i = 0; i<(n+1); i++){
            //cout<<"i: "<<i;
            if(i>0){
                for(int policy = 0; policy<im->mdp->P; policy++){
                    current_v.at(policy) = im->score()[policy][i-1];
                }
                //printVector(current_v);
                //cout<<"policy<im->mdp->P "<<im->mdp->P<<endl;
            }
            for(int w = 0; w<(budget+1); w++){
                //cout<<" and w: "<<w<<endl;
                //cout<<"-------------------"<<endl;
                if(i == 0 || w == 0){
                    K[i][w].push_back({});
                    for(int policy = 0; policy<im->mdp->P; policy++){
                         K[i][w][0].push_back(0);
                    }
                    K2[i][w] = 0;
                }
                else if(im->mdp->findCost(i-1) <= w){
                    K[i][w] = K[i-1][w];
                    K2[i][w] = K2[i-1][w];
                    
                    //for(int z=0; z<K[i-1][w-im->mdp->findCost(i-1)].size(); z++){
                    int zupperbound = 1;
                    int zmax = zupperbound;
                    if(K[i-1][w-im->mdp->findCost(i-1)].size()<zmax) zmax = K[i-1][w-im->mdp->findCost(i-1)].size();
                    for(int z=0; z<zmax; z++){
                        vector<float> current_to_check = SUM_vectors(&current_v, &K[i-1][w-im->mdp->findCost(i-1)][z]);
                        float new_min =  min(&current_to_check);
                        if(zmax<zupperbound){
                            if(new_min > K2[i][w]){
                                K2[i][w] = new_min;
                            }
                            K[i][w].push_back(current_to_check);
                        }
                        else{
                            int worst_j = 0;
                            float minall = min(& K[i][w][0]);
                            for(int jj=0; jj<K[i][w].size(); jj++){
                                float minj =  min(& K[i][w][jj]);
                                if(minj<minall){
                                    minall = minj;
                                    worst_j = jj;
                                }
                            }
                            if(new_min>minall){
                                K[i][w][worst_j] = current_to_check;
                                if(new_min> K2[i][w]){
                                    K2[i][w] = new_min;
                                }
                            }
                        }
                    }
                }
                else{
                    K[i][w] = K[i-1][w];
                    K2[i][w] = K2[i-1][w];
                }
            }
        }
        solution.second = K2[n][budget];
        int w = budget;
        for(int i=n; i>=1; i--){
            if(K2[i][w] != K2[i-1][w]){
                solution.first.push_back(1);
                w -= im->mdp->findCost(i-1);
            }
            else{
                solution.first.push_back(0);
            }
        }
        for (int i=0; i<n+1; i++){
            delete[] K2[i];
            delete[] K[i];
        }
        delete K2;
        delete K;
    }
    reverse(solution.first.begin(),solution.first.end());
    return solution;
}
 
