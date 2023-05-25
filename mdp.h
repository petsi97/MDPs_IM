#pragma once 
using namespace std;
#include "mc.h"
#include <unordered_set>


/*
Implementation of an MDP:
@param N: number of nodes
@param P: number of policies
@param MCs: different markov chains of the MDP
@param Cost: cost for each different state
@param Reward: reward for each state
@param total_budget: sum of all costs
*/
class MDP{
    public:
        int N, P;
        double gamma;
        vector<MC>* MCs;     //pointer to a vector of MCs
        vector<int>* Cost;   //pointer to a vector of Costs
        vector<int>* Reward; //pointer to a vector of Rewards
        int total_budget;


        //Empty Constructor
        MDP(){
            printf("Initialize an MDP...\n");
            this->N = 0;
            this->P = 0;
            this->gamma = 1.0;
            this->total_budget = 0;
            this->MCs = new vector<MC>();
            this->Cost = new vector<int>();
            this->Reward = new vector<int>();
        }

        //Constructor that initialize the MDP with the given MCs
        //@params[in] MCs: the MCs
        MDP(vector<MC>* MCs){
            printf("Initialize an MDP with multiple Markov Chains...\n");
            this->gamma = 1.0;
            this->MCs = MCs;
            this->Cost = new vector<int>();
            this->Reward = new vector<int>();
            this->P = this->MCs->size();
            this->total_budget = 0;
            this->N = 0;
            if(this->P > 0){
                this->N = this->MCs->begin()->N;
            }
            for(int i = 0; i < this->N; i++){
                this->Reward->emplace_back(0);
            }
            this->setUniformCosts();
        }

        //My Destructor
        void my_destructor(){
            for(auto mc:*MCs) mc.my_destructor();
            delete this->Cost;
            delete this->Reward;
        }
        

        //Set uniform costs (equal to one)
        void setUniformCosts() {
            for(int i = 0; i < this->N; i++){
                this->Cost->emplace_back(1);
            }
        }

        //Set a given reward to the given state
        //@params[in] state: state to update its reward
        //@params[in] reward: the new reward
        void setRewardState(int state, int reward){
            (*(this->Reward))[state] = reward;
        }

        //Set all states as non-reward states (empty reward states)
        void emptyRewardStates(){
            for(int i=0; i<this->N; i++){
                (*(this->Reward))[i] = 0;
            }
        }

        //Set all states as reward states
        void fillRewardStates(){
            for(int i=0; i<this->N; i++){
                (*(this->Reward))[i] = 1;
            }
        }
        
        //Read costs
        //@param[in] file: filename to read the costs for each state
        //@throw Throws an std::invalid_argument error if the number of states differ from the ones in the file
        void read_costs(const string& file) {
            if(this->N != read_N(file)) throw std::invalid_argument("Cost size doen't match the number of States!!!" );
            this->Cost->clear();
            int counter = 0;
            ifstream input(file);
            for (string line; getline(input, line);) {
                counter++;
                if(counter == 1) continue;
                int cost = stoi(line);
                this->total_budget += cost;
                this->Cost->emplace_back(cost);
                counter++;
            }
        }

        //Find cost of a reward state
        //@param[in] s: reward states
        //@return Cost of state s
        int findCost(int s){
            return (*Cost)[s];
        }

        //Find the total cost of the given reward states
        //@param[in] Reward: pointer to reward states
        //@return Total cost
        int findCost(vector<int>* Reward){
            int cost = 0;
            for(int i=0; i<Reward->size(); i++){
                cost += (*Cost)[i]*((*Reward)[i]);
            }
            return cost;
        }

        //Find the total cost of the given reward states
        //@param[in] Reward: pointer to reward states
        //@return Total cost
        int findCost(unordered_set<int>* Reward){
            int cost = 0;
            for(auto reward: *Reward){
                cost += reward;
            }
            return cost;
        }

        //Print all data
        void printAll(){
            for(int i=0; i<this->P; i++){
                this->MCs->at(i).printEdges();
                printf("----------------------------------------\n");
            }
        }

        //Print costs per node
        void printCosts(){
            printf("Printing Costs: \n");
            for(int i=0; i<this->N; i++){
                printf("For %d: %d\n", i, (*(this->Cost))[i]);
            }
        }

        //Print rewards per node
        void printRewards(){
            printf("Printing Rewards: \n");
            for(int i=0; i<this->N; i++){
                printf("For %d: %d\n", i, (*(this->Reward))[i]);
            }
        }      

};
