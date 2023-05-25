#pragma once 
using namespace std;
#include <cmath>
#include "mdp.h"
#include "structs.h"
#include <set>
#include <queue>
#include<ios> //used to get stream size
#include<limits>
#include <unordered_set>
#include <unordered_map>
#include <numeric> // for accumulate




void* find_thread_scores(void *par){
    //cout<<"find scores"<<endl;
    unordered_set<int> reward_states;
    //float **S;
    //int **visited;
    int N = ((param *)par)->mdp->N;
    //S: for each state we save [current_value, updated_value]
    ((param *)par)->S = new float *[N];
    for (int i = 0; i < N; i++) ((param *)par)->S[i] = new float[2];
    //visited: for each state we save [current_step, previous_step]
    ((param *)par)->visited = new int *[N];
    for (int i = 0; i < N; i++) ((param *)par)->visited[i] = new int[2];
    for(int policy = 0; policy<((param *)par)->mdp->P; policy++){
        int new_counter = 0;
        for(auto s:((param *)par)->to_compute){
            reward_states.insert(s);
            double influence = 0.0;
            int k = ((param *)par)->mdp->MCs->at(policy).Kmax;
            for(int i=0; i<((param *)par)->mdp->N; i++){
                //visited | current_step = k = remaining steps
                //visited | previous_step = k+1 = the total steps
                ((param *)par)->visited[i][0] = k;
                ((param *)par)->visited[i][1] = k+1;
                if(reward_states.find(i) != reward_states.end()){
                    influence += (*(((param *)par)->mdp->MCs))[policy].I[i];
                    //current reward is 1 if it is in the reward states to be checked 
                    ((param *)par)->S[i][0] = 1;
                    ((param *)par)->S[i][1] = 0;
                }
                else{
                    //current reward is 0 if it is not in the reward states to be checked 
                    ((param *)par)->S[i][0] = 0;
                    ((param *)par)->S[i][1] = 0;
                }
            }     
            //Algorithms starts from the reward states and calculates influence backwards by following reverse edges.
            //Init a queue, saving the state to be checked and the remaining steps
            queue<node_step> q;
            for(auto state: reward_states){
                q.push(node_step{state, k});
            }
            //Variable to check the remaining steps
            int current_k = k;
            //While queue is not empty
            
            while(!q.empty()){
                //next state to be examined and the remaining steps from this state
                node_step nps = q.front(); q.pop();
                //in neighs of node
                int in_neighs = (*(((param *)par)->mdp->MCs))[policy].RevT[nps.node].size();
                //if there are no remaining steps else updated current remaining steps if needed
                if(nps.step==0) break;
                else if(nps.step < current_k){
                    current_k = nps.step;
                }
                //for each of the in-neghbors of the state
                for(int i=0; i<in_neighs; i++){
                    Edge edge = (*(((param *)par)->mdp->MCs))[policy].RevT[nps.node][i];
                    int from_state = edge.first; 
                    //if the in-neighbor has already been visited for the current step then we should 
                    //not consider into account that state, thus we continue with the next state in the queue.
                    if(((param *)par)->visited[from_state][0] <= current_k-1){
                        continue;
                    }
                    //the current score for inclunding from_state
                    double current_score = 0.0;
                    int out_neighs = (*(((param *)par)->mdp->MCs))[policy].T[from_state].size();
                    for(int j=0; j<out_neighs; j++){
                        edge = (*(((param *)par)->mdp->MCs))[policy].T[from_state][j];
                        int to_state = edge.first;
                        double prob = edge.second;
                        if(((param *)par)->visited[to_state][0] == nps.step){
                            current_score += prob*(((param *)par)->S[to_state][0]);
                        }
                        else if(((param *)par)->visited[to_state][1] == nps.step){
                            current_score += prob*(((param *)par)->S[to_state][1]);
                        }
                    }
                    //update S and visited
                    ((param *)par)->visited[from_state][1] = ((param *)par)->visited[from_state][0];
                    ((param *)par)->S[from_state][1] = ((param *)par)->S[from_state][0];
                    ((param *)par)->visited[from_state][0] = nps.step-1;
                    ((param *)par)->S[from_state][0] = current_score;
                    //insert from state to queue
                    q.push(node_step{from_state, nps.step-1});
                    //increase influence
                    influence += ((param *)par)->mdp->MCs->at(policy).I[from_state] * ((param *)par)->mdp->MCs->at(policy).K[from_state][nps.step-1]*current_score;
                }    
            }
            reward_states.clear();
            (((param *)par)->scores)[policy][new_counter] = influence;
            new_counter++;     
        }
    } 
    for (int i=0; i<N; i++){
        delete[] ((param *)par)->S[i];
        delete ((param *)par)->visited[i];
    }
    //cout<<"***end***"<<endl;
    delete ((param *)par)->S;
    delete ((param *)par)->visited;
    pthread_exit(NULL);
    return NULL;                
};
