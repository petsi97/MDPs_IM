#pragma once 
using namespace std;
#include "mdp.h"

//struct for node, step
struct node_step {
  int node;
  int step;
};

//struct for score, step
struct score_step {
  double score;
  int step;
};

//struct for node, prob, step
struct node_prob_step {
  int node;
  double prob;
  int step;
};

//struct used for threading
struct param {
    vector<int> to_compute;
    MDP* mdp;
    float** scores;
    float **S;
    int ** visited;
};