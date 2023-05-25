#pragma once 
using namespace std;
#include<stdio.h>
#include <iostream>
#include <vector>
#include "functionality.h"


/*
Implementation of a Markov Chain:
@param N: number of nodes
@param M: number of edges
@param Kmax: max number of steps in the MC
@param I: initial state probabilities
@param T: out-going transition probabilities
@param RevT: in-going transition probabilities
@param K: step probabilities
*/
class MC{
    public:
        int N, M, Kmax;
        vector<double> I;
        vector<Edge>* T;
        vector<Edge>* RevT;
        vector<double>* K;
        
        //Empty Constructor
        MC(){
            printf("Initialize Markov Chain...");
            this->N = 0;
            this->M = 0;
            this->T = new vector<Edge>[1];
            this->RevT = new vector<Edge>[1];
            this->K = new vector<double>();
        }

        //Constructor reading from file
        MC(const string& file){
            printf("Initialize Markov Chain reading %s...\n", file.c_str());
            read_graph(file);
            this->K = new vector<double>();
            this->Kmax = 2;
            setUniformSteps();
        }

        //MC(MC const &) = delete;
        //void operator=(MC const &) = delete;

        //Destructor
        void my_destructor(){
            delete[] this->T;
            delete[] this->RevT;
        }
        

        //Function that sets uniform steps i.e. if Kmax = 3 then probabilities are [1/3, 1/3, 1/3]
        //Note that we save probabilites [1/3, 2/3, 3/3] instead of [1/3, 1/3, 1/3].
        void setUniformSteps(){
            if(this->Kmax <= 0) return;
            this->K = new vector<double>[this->N];
            for(int i = 0; i<this->N; i++){
                for(int j = 0; j<this->Kmax; j++){
                    this->K[i].emplace_back((j+1)*1.0/Kmax);
                }
            } 
        }

        //Function that sets steps according to normal distribution i.e. if Kmax = 5 then probabilities are [1/6, 2/6, 3/6, 2/6, 1/6]
        void setNormalDistribution(){
            if(this->Kmax <= 0) return;
            this->K = new vector<double>[this->N];
            vector<double> steps;
            for(int j = 0; j<this->Kmax; j++){
                if(j<=((Kmax)/2))
                    steps.emplace_back((j+1));
                else
                    steps.emplace_back(Kmax+1 - (j+1));
            } 
            this->setSteps(steps);
        }

        //Set Kmax = k and set uniform steps
        //@param[in] k: steps
        void setUniformSteps(int k){
            this->Kmax = k;
            setUniformSteps();
        }

        //Set Kmax = k and set steps according to normal distribution
        //@param[in] k: steps
        void setNormalDistribution(int k){
            this->Kmax = k;
            setNormalDistribution();
        }

        //Normalize vector to probabilities
        //@param[in] k: steps
        void setSteps(vector<double>& steps){
            this->Kmax = steps.size();
            this->K = new vector<double>[this->N];
            double sumall = 0.0;
            for(auto step: steps) sumall+=step;
            for(int i = 0; i<this->N; i++){
                this->K[i].emplace_back(steps[0]/sumall);
                for(int j = 1; j<this->Kmax; j++){
                    this->K[i].emplace_back(this->K[i][j-1]+steps[j]/sumall);
                }
            } 
        }

        //Function that reads the MC from a file
        //@param[in] file: filename to read the MC
        void read_graph(const string& file) {
            //read num of nodes (first row of text)
            this->N = read_N(file); 
            this->M = 0;
            this->T = new vector<Edge>[this->N];
            this->RevT = new vector<Edge>[this->N];
            for(int i = 0; i<this->N; i++) this->I.emplace_back(0);
            string line;
            ifstream input(file);
            getline(input, line);
            for (string line; getline(input, line);) {
                //Three numbers per line
                auto [a, b, c] = read_line(line);    
                //Here we assume no self-loops are allowed (self-weights = initial prob.)  
                if (a == b) {
                    this->I[a] = c;
                } else {
                    this->M++;
                    this->T[a].emplace_back(b, c);
                }
            }
            input.close();
        }

        //Normalize edge weigths
        void normalizeEdgeWeight(){
            this->M = 0;
            for(int i = 0; i < this->N; i++){
                if(this->T[i].size() == 0){
                    this->T[i].emplace_back(i, 1);
                    this->RevT[i].emplace_back(i, 1);
                    this->M++;
                }
                else{
                    double sum = 0;
                    for(auto edge: this->T[i]) sum += edge.second;
                    //if an edge does not have out-going edges then we set a self-loop
                    if(sum==0){
                        this->T[i].clear();
                        this->T[i].emplace_back(i, 1);
                        this->RevT[i].emplace_back(i, 1);
                        this->M++;
                        continue;
                    }
                    for(int j = 0; j < this->T[i].size(); j++){ 
                        this->T[i][j].second /= sum;
                        this->RevT[this->T[i][j].first].emplace_back(i, this->T[i][j].second);
                        this->M++;
                    }
                }
            }
            this->normalizeInitialProbs();
        }

        //Normalize initial probabilities
        void normalizeInitialProbs(){
            int sum = 0;
            for(int i = 0; i < this->N; i++) sum += this->I[i];
            for(int i = 0; i < this->N; i++) this->I[i] /= sum;
        }

        //Set uniform initial probabilities
        void set_uniform_distribution(){
            if(N==0) return;
            for(int i=0; i<this->N; i++){
                this->I[i] = 1.0/N;
            }
        }

        //Print nodes
        void printNodes(){
            printf("There are in total %d nodes: ",N);
            for(int i = 0; i < this->N; i++){
                printf("(%d) ",i);
            }
            printf("\n");
        }

        //Print edges
        void printEdges(){
            this->printNodes();
            for(int i = 0; i < this->N; i++){
                printf("For node %d: ", i);
                printVector(this->T[i]);
            }
            printf("\n");
        }

        //Print in-coming probabilities
        void printReverseEdges(){
            this->printNodes();
            for(int i = 0; i < this->N; i++){
                printf("For node %d: ", i);
                printVector(this->RevT[i]);
            }
            printf("\n");
        }

        //Print steps
        void printSteps(){
            this->printNodes();
            for(int i = 0; i < this->N; i++){
                printf("For node %d steps are:", i);
                printVector(this->K[i]);
            }
            printf("\n");
        }

};