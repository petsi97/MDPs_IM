#pragma once 
using namespace std;
#include <cmath>
#include "mdp.h"
#include "thread_functionality.h"
#include "structs.h"
#include <set>
#include <queue>
#include<ios> //used to get stream size
#include<limits>
#include <unordered_set>
#include <unordered_map>
#include <numeric> // for accumulate



/*
Implementation base class called IM. 
This class is used to compute (and/or save) the influence of a state reward or a collection of state rewards.
@param mdp: an MDP 
*/
class IM{
    public:
        MDP* mdp; 
    
    IM(){
        this->mdp = new MDP();
    }

    IM(MDP* model){
        this->mdp = model;
    }

    void my_destructor(){
        mdp->my_destructor();
    }

    virtual vector<double>* score(){
        return new vector<double>();
    }

    virtual double all_scores(int policy){
        return 0.0;
    }

    virtual double F(int policy){return 0;}

    virtual void setReward(int reward_state, int reward){this->mdp->setRewardState(reward_state, reward);}

    virtual void setRewards(vector<int>& rewards){*(this->mdp->Reward) = rewards;}

    virtual void emptyRewardStates(){this->mdp->emptyRewardStates();}

    virtual double f(int policy, unordered_set<int> &reward_states){return 0;}

    virtual double G(double current_reward, double state_reward, int t, int k){return 0;}

};

/*
Implementation of class that inherits from IM. 
This class is used to compute (and/or save) the influence and change (addremove) 
the reward states or a collection of state rewards according to the IMPRESSION COUNTS.
@param mdp: an MDP 
@param scores: individual score for each state and for each policy
@param current_score: current score of the selected reward states for each policy 
*/
class IM_IC : public IM {
    using IM::IM;
    public:
        vector<double>* scores;
        vector<double> current_score;

        //Empty Constructor
        IM_IC(){
            this->mdp = new MDP();
            scores = new vector<double>();
        }

        //Constructor that uses the given MDP
        //@param model: an MDP 
        IM_IC(MDP* model){
            this->mdp = model;
            scores = new vector<double>[this->mdp->P];
            cout<<"Calculating scores..."<<endl;
            //this->find_scores();
            this->find_scores_multithreading();
        }

        //Returns the saved scores 
        //@return model: vector<double>* dynamic array (for each policy) of 
        //vectors (for each state) where scores are saved
        vector<double>* score(){
            return this->scores;
        }

        //Calculating the scores
        void find_scores(){
            //cout<<"find scores"<<endl;
            unordered_set<int> reward_states;
            //S: for each state we save [current_value, updated_value]
            float **S;
            S = new float *[this->mdp->N];
            for (int i = 0; i < this->mdp->N; i++) S[i] = new float[2];
            //visited: for each state we save [current_step, previous_step]
            int **visited;
            visited = new int *[this->mdp->N];
            for (int i = 0; i < this->mdp->N; i++) visited[i] = new int[2];
            for(int s = 0; s<mdp->N; s++){
                for(int policy = 0; policy<mdp->P; policy++){
                    //if((s%100 == 0)) cout<<s<<" - "<<policy<<endl;
                    if(s==0){
                        current_score.push_back(0);
                    }
                    //cout<<"Calculating score for policy "<<policy<<" and state "<<s<<endl;
                    reward_states.insert(s);
                    double score =  f(policy, reward_states, S, visited);
                    reward_states.clear();
                    //if((s%100 == 0)) cout<<"score "<<score<<endl;
                    //cout<<"score is "<<score<<endl;
                    scores[policy].push_back(score);
                }
            } 
            for (int i=0; i<this->mdp->N; i++){
                delete[] S[i];
                delete visited[i];
            }
            delete S;
            delete visited;
        }
        
        //Calculating the scores
        void find_scores_multithreading(){
            //cout<<"For multithreding"<<endl;
            int num_threads = 6;
            unsigned int threadIter = (1.0*this->mdp->N)/num_threads;
            unsigned int remainIter = this->mdp->N%num_threads;
            pthread_t thread_id[num_threads];
            param pars[num_threads];
            int counter = 0;
            for(int i = 0; i < num_threads; i++){
                pars[i].mdp = this->mdp;
                if(i<remainIter){
                    for(int j=0; j<threadIter+1; j++){
                        pars[i].to_compute.push_back(counter);
                        counter++;
                    }
                    pars[i].scores = new float *[this->mdp->P];
                    for (int j = 0; j < this->mdp->P; j++) pars[i].scores[j] = new float[threadIter+1];
                }
                else{
                    for(int j=0; j<threadIter; j++){
                        pars[i].to_compute.push_back(counter);
                        counter++;

                    }
                    pars[i].scores = new float *[this->mdp->P];
                    for (int j = 0; j < this->mdp->P; j++) pars[i].scores[j] = new float[threadIter+1];
                }
                pthread_create(&thread_id[i], NULL, find_thread_scores, (void *)&(pars[i]));
            }
            //cout<<"We are here"<<endl;
            for(int j=0; j < num_threads; j++) 
                pthread_join( thread_id[j], NULL); 
            //cout<<"We end here"<<endl;
            counter = 0;
            for(unsigned int i = 0; i < num_threads; i++)
            {
                for (int policy = 0; policy < this->mdp->P; policy++){
                    if(i==0) current_score.push_back(0);
                    for(int j=0; j<pars[i].to_compute.size();j++){
                        this->scores[policy].push_back(pars[i].scores[policy][j]);
                    }
                    delete[] pars[i].scores[policy];
                }
                delete pars[i].scores;
            }
        }

        //The score for a given policy 
        //@param policy: the policy for which we want to find the score
        //@return The influence score
        double F(int policy){
            return current_score[policy];
        }

        //Function that calculates the score for a given policy and reward states
        //@param policy: the policy for which we want to find the score
        //@param reward_states: the states for which we want to compute the score
        //@return The total score
        double f(int policy, unordered_set<int> &reward_states, float **S, int **visited){
            double influence = 0.0;
            int k = this->mdp->MCs->at(policy).Kmax;
            for(int i=0; i<this->mdp->N; i++){
                //visited | current_step = k = remaining steps
                //visited | previous_step = k+1 = the total steps
                visited[i][0] = k;
                visited[i][1] = k+1;
                if(reward_states.find(i) != reward_states.end()){
                    influence += (*(this->mdp->MCs))[policy].I[i];
                    //current reward is 1 if it is in the reward states to be checked 
                    S[i][0] = 1;
                    S[i][1] = 0;
                }
                else{
                    //current reward is 0 if it is not in the reward states to be checked 
                    S[i][0] = 0;
                    S[i][1] = 0;
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
                int in_neighs = (*(this->mdp->MCs))[policy].RevT[nps.node].size();
                //if there are no remaining steps else updated current remaining steps if needed
                if(nps.step==0) break;
                else if(nps.step < current_k){
                    current_k = nps.step;
                }
                //for each of the in-neghbors of the state
                for(int i=0; i<in_neighs; i++){
                    Edge edge = (*(this->mdp->MCs))[policy].RevT[nps.node][i];
                    int from_state = edge.first; 
                    //if the in-neighbor has already been visited for the current step then we should 
                    //not consider into account that state, thus we continue with the next state in the queue.
                    if(visited[from_state][0] <= current_k-1){
                        continue;
                    }
                    //the current score for inclunding from_state
                    double current_score = 0.0;
                    int out_neighs = (*(this->mdp->MCs))[policy].T[from_state].size();
                    for(int j=0; j<out_neighs; j++){
                        edge = (*(this->mdp->MCs))[policy].T[from_state][j];
                        int to_state = edge.first;
                        double prob = edge.second;
                        if(visited[to_state][0] == nps.step){
                            current_score += prob*(S[to_state][0]);
                        }
                        else if(visited[to_state][1] == nps.step){
                            current_score += prob*(S[to_state][1]);
                        }
                    }
                    //update S and visited
                    visited[from_state][1] = visited[from_state][0];
                    S[from_state][1] = S[from_state][0];
                    visited[from_state][0] = nps.step-1;
                    S[from_state][0] = current_score;
                    //insert from state to queue
                    q.push(node_step{from_state, nps.step-1});
                    //increase influence
                    influence += this->mdp->MCs->at(policy).I[from_state] * this->mdp->MCs->at(policy).K[from_state][nps.step-1]*current_score;
                }    
            }
            return influence;
        }

        //This function can be modified to express different (additive) function.
        //In its current form every time we visit a state we increase our cumulative reward by the state reward.
        //@param current_reward: currently total reward
        //@param state_reward: the reward of the new state
        //@param t: timestamp
        //@param k: remaining steps
        //@return The new total reward
        double G(double current_reward, double state_reward, int t, int k){
            return current_reward + state_reward; //pow(this->mdp->gamma, (k - t)) * state_reward;
        }

        //Function that changes the reward of the reward_state. We have made the assumption that the 
        //reward can be either 1 or 0, corresponding to be or not a reward state respectively.
        //@param reward_state: state to be updated
        //@param reward: new reward
        virtual void setReward(int reward_state, int reward){
            this->mdp->setRewardState(reward_state, reward);
            for(int policy = 0; policy<mdp->P; policy++){
                if(reward == 0){
                    current_score[policy] -= this->scores[policy].at(reward_state);
                }
                else{
                    current_score[policy] += this->scores[policy].at(reward_state);
                }
            }
        }

        //Function that sets new reward_states.
        //@param rewards: states to be reward states
        virtual void setRewards(vector<int>& rewards){
            //this->mdp->emptyRewardStates();
            for(int policy = 0; policy<mdp->P; policy++) current_score[policy] = 0;
            for(int i=0; i<rewards.size(); i++){
                this->mdp->setRewardState(i, rewards[i]);
                if(rewards[i] == 0) continue;
                for(int policy = 0; policy<mdp->P; policy++){
                    current_score[policy] += this->scores[policy].at(i);
                }
                
            }
        }

        //Empty reward states.
        virtual void emptyRewardStates(){
            this->mdp->emptyRewardStates();
            for(int policy = 0; policy<mdp->P; policy++){
                current_score[policy] = 0;
            }
        }
};

/*
Implementation of class that inherits from IM. 
This class is used to compute (and/or save) the influence and change (addremove) 
the reward states or a collection of state rewards according to the HITTING PROBABILITIES.
@param mdp: an MDP 
@param reward_states: the reward states
*/
class IM_HT : public IM {
    using IM::IM;
    public:
        unordered_set<int> reward_states;
        float **S;
        int **visited;


        IM_HT(){
            this->mdp = new MDP();
            //S: for each state we save [current_value, updated_value]
            this->S = new float *[this->mdp->N];
            for (int i = 0; i < this->mdp->N; i++) S[i] = new float[2];
            //visited: for each state we save [current_step, previous_step]
            this->visited = new int *[this->mdp->N];
            for (int i = 0; i < this->mdp->N; i++) visited[i] = new int[2];
        }

        IM_HT(MDP* model){
            this->mdp = model;
            //S: for each state we save [current_value, updated_value]
            this->S = new float *[this->mdp->N];
            for (int i = 0; i < this->mdp->N; i++) S[i] = new float[2];
            //visited: for each state we save [current_step, previous_step]
            this->visited = new int *[this->mdp->N];
            for (int i = 0; i < this->mdp->N; i++) visited[i] = new int[2];
        }

        void my_destructor(){
            for (int i=0; i<this->mdp->N; i++){
                delete[] this->S[i];
                delete this->visited[i];
            }
            delete this->S;
            delete this->visited;
            delete this->mdp;
        }

        //The score for a given policy 
        //@param policy: the policy for which we want to find the score
        //@return The influence score
        double F(int policy){
            return f(policy, reward_states);
        }

        //Function that changes the reward of the reward_state. We have made the assumption that the 
        //reward can be either 1 or 0, corresponding to be or not a reward state respectively.
        //@param reward_state: state to be updated
        //@param reward: new reward
        virtual void setReward(int reward_state, int reward){
            this->mdp->setRewardState(reward_state, reward);
            if(reward == 0) this->reward_states.erase(reward_state);   
            if(reward > 0) this->reward_states.insert(reward_state);   
        }

        //Function that sets new reward_states.
        //@param rewards: states to be reward states
        virtual void setRewards(vector<int>& rewards){
            this->emptyRewardStates();
            for(int i=0; i<rewards.size(); i++){
                this->mdp->setRewardState(i, rewards[i]);
                if(rewards[i] == 0) continue;
                this->reward_states.insert(i);   
            }
        }
        
        //Empty reward states.
        virtual void emptyRewardStates(){
            this->mdp->emptyRewardStates();
            this->reward_states.clear();
        }

        //Function that calculates the score for a given policy and reward states
        //@param policy: the policy for which we want to find the score
        //@param reward_states: the states for which we want to compute the score
        //@return The influence score
        double f(int policy, unordered_set<int> &reward_states){
            double influence = 0.0;
            int k = this->mdp->MCs->at(policy).Kmax;
            for(int i=0; i<this->mdp->N; i++){
                //visited | current_step = k = remaining steps
                //visited | previous_step = k+1 = the total steps
                this->visited[i][0] = k;
                this->visited[i][1] = k+1;
                if(reward_states.find(i) != reward_states.end()){
                    influence += (*(this->mdp->MCs))[policy].I[i];
                    //current reward is 1 if it is in the reward states to be checked 
                    this->S[i][0] = 1;
                    this->S[i][1] = 0;
                }
                else{
                    //current reward is 0 if it is not in the reward states to be checked 
                    this->S[i][0] = 0;
                    this->S[i][1] = 0;
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
                int in_neighs = (*(this->mdp->MCs))[policy].RevT[nps.node].size();
                //if there are no remaining steps else updated current remaining steps if needed
                if(nps.step==0) break;
                if(nps.step < current_k){
                    current_k = nps.step;
                }
                //for each of the in-neghbors of the state
                for(int i=0; i<in_neighs; i++){
                    Edge edge = (*(this->mdp->MCs))[policy].RevT[nps.node][i];
                    int from_state = edge.first;
                    //if the in-neighbor has already been visited for the current step or it is a reward (absorption) state 
                    //then we should not consider into account that state, thus we continue with the next state in the queue.
                    if((this->visited[from_state][0] <= current_k-1)||(reward_states.find(from_state) != reward_states.end())){
                        continue;
                    }
                    //the current score for inclunding from_state
                    double current_score = 0.0;
                    int out_neighs = (*(this->mdp->MCs))[policy].T[from_state].size();
                    for(int j=0; j<out_neighs; j++){
                        edge = (*(this->mdp->MCs))[policy].T[from_state][j];
                        int to_state = edge.first;
                        double prob = edge.second;
                        if(this->visited[to_state][0] == nps.step){
                            current_score += prob*(S[to_state][0]);
                        }
                        else if(this->visited[to_state][1] == nps.step){
                            current_score += prob*(S[to_state][1]);
                        }           
                    }
                    //update S and visited
                    this->visited[from_state][1] = this->visited[from_state][0];
                    this->S[from_state][1] = this->S[from_state][0];
                    this->visited[from_state][0] = nps.step-1;
                    this->S[from_state][0] = current_score;
                    //insert from state to queue
                    q.push(node_step{from_state, nps.step-1});
                    //increase influence
                    influence += this->mdp->MCs->at(policy).I[from_state] * this->mdp->MCs->at(policy).K[from_state][nps.step-1]*current_score;
                }
            }
            return influence;
        }

};