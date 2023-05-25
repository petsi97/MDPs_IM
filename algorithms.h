#pragma once 
using namespace std;
#include "rim.h"
#include <queue>
#include <unordered_map>
#include <limits>
#include <random>
#include <algorithm>
#include "knapsack.h"
#include<stdio.h>


//Random int 
//@param[in] i: max int - 1
//@return Random int [0, max-1)
int myrandom (int i) { return std::rand()%i;}



//Init a priority queue for the lazy greedy method
//@param[in] im: IM 
//@param[in] budget: budget constraint
//@param[in] policy: policy
//@param[in,out] V_minus_S: states that not included in the solution (or exceed the constraint)
//@param[in,out] current_cost: current cost
//@param[in,out] current_score: current score
//@param[in,out] pq: the priority queue
void init_priority_queue(IM *im, int budget, int policy, unordered_set<int> *V_minus_S, int *current_cost, 
                         double *current_score, priority_queue<pair<double, int> > *pq){
    //-----------------------------------------------------------------------------------
    for (auto s = V_minus_S->begin(); s != V_minus_S->end(); ){
        int cost_s = im->mdp->findCost(*s);
        if ((cost_s + *current_cost) > budget){
            V_minus_S->erase(s++);
            continue;
        }
        else{
            im->setReward(*s, 1);
            pq->push(make_pair(im->F(policy)/cost_s, *s)); //Both are ok marginal/cost_s and score/cost_s
            im->setReward(*s, 0);
            ++s;
        }
    }
    //-----------------------------------------------------------------------------------
    if(pq->empty()) return;
    im->setReward(pq->top().second, 1);

    int cost_top_s = im->mdp->findCost(pq->top().second); 
    *current_cost += cost_top_s;

    double ret_score = pq->top().first * cost_top_s; //This is due to the division by cost
    *current_score = ret_score;                      //Current score

    V_minus_S->erase(pq->top().second); pq->pop();
    //-----------------------------------------------------------------------------------   
}



//One step of the lazy greedy method
//@param[in] im: IM 
//@param[in] budget: budget constraint
//@param[in] policy: policy
//@param[in,out] V_minus_S: states that not included in the solution (or exceed the constraint)
//@param[in,out] current_cost: current cost
//@param[in,out] current_score: current score
//@param[in,out] pq: the priority queue
void lazy_step(IM *im, int budget, int policy, unordered_set<int> *V_minus_S, int *current_cost, 
               double *current_score, priority_queue<pair<double, int> > *pq){
    //-----------------------------------------------------------------------------------
    unordered_map<int, double> checked;
    while(pq->size()>0){
        int top_s = pq->top().second;
        pq->pop();
        if(V_minus_S->find(top_s) == V_minus_S->end()) continue;
        float cost_s = im->mdp->findCost(top_s);
        if ((cost_s + *current_cost) > budget){
            V_minus_S->erase(top_s);
            continue;
        }
        if(checked.find(top_s) != checked.end()){
            im->setReward(top_s, 1);
            *current_cost += cost_s;
            *current_score = checked[top_s]; //Current score
            V_minus_S->erase(top_s);
            return;
        }
        im->setReward(top_s, 1);
        double score = im->F(policy);
        pq->push(make_pair(score/cost_s, top_s));
        im->setReward(top_s, 0);
        checked[top_s] = score; //Score after including node
    }
    //-----------------------------------------------------------------------------------
}



//Lazy Greedy implementation
//@param[in] im: IM 
//@param[in] budget: budget constraint
//@param[in] policy: policy
//@return pair<vector<int>, double> pair of solution with the influence score
pair<vector<int>, double> lazyGreedy(IM *im, int budget, int policy){
    //-----------------------------------------------------------------------------------
    im->emptyRewardStates();
    unordered_set<int> *V_minus_S = new unordered_set<int>();
    for(int i=0; i<im->mdp->N; i++) V_minus_S->insert(i);
    int *current_cost = new int();
    double *current_score = new double();
    *current_cost = 0.0;
    *current_score = 0.0;
    priority_queue<pair<double, int> > *pq;
    pq = new priority_queue<pair<double, int> >();
    //-----------------------------------------------------------------------------------
    init_priority_queue(im, budget, policy, V_minus_S, current_cost, current_score, pq);
    while(V_minus_S->size()>0){
        lazy_step(im, budget, policy, V_minus_S, current_cost, current_score, pq);
    }
    double score = *current_score;
    //-----------------------------------------------------------------------------------
    V_minus_S->clear(); delete V_minus_S;
    delete current_cost; delete current_score;
    *pq = priority_queue<pair<double, int> >(); delete pq;

    return make_pair(*(im->mdp->Reward), score);
    //-----------------------------------------------------------------------------------
}



//Find optimal (best) solution (or approximate it)
//@param[in] im: IM 
//@param[in] budget: budget constraint
//@param[in] policy: policy
//@param[in] ret_sol: boolean variable that determines if we return the solution or not
//@return pair<vector<int>, double> pair of solution with the influence score
pair<vector<int>, double> findOptimal(IM *im, int budget, int policy, bool ret_sol){
    //-----------------------------------------------------------------------------------
    if(im->score()->size()>0){
        return knapSack(im, budget, policy, ret_sol, false);
    }
    return lazyGreedy(im, budget, policy);
    //-----------------------------------------------------------------------------------
}



//Find non-optimal (worst) solution (or approximate it)
//@param[in] im: IM 
//@param[in] budget: budget constraint
//@param[in] policy: policy
//@param[in] ret_sol: boolean variable that determines if we return the solution or not
//@return pair<vector<int>, double> pair of solution with the influence score
pair<vector<int>, double> findNonOptimal(IM *im, int budget, int policy, bool ret_sol){
    //-----------------------------------------------------------------------------------
    pair<vector<int>, double> ret;
    ret.second = 0;
    return ret; //comment to use it
    if(im->score()->size()>0){
        return knapSack(im, im->mdp->total_budget - budget, policy, ret_sol, true);
    }
    pair<vector<int>, double> to_ret = lazyGreedy(im, im->mdp->total_budget - budget, policy);
    to_ret.second = lazyGreedy(im, im->mdp->total_budget, policy).second - to_ret.second;
    if(!ret_sol) return to_ret;
    for(int i=0; i<to_ret.first.size(); i++){
        if(to_ret.first[i]==0) to_ret.first[i]=1;
        else to_ret.first[i]=0;
    }
    return to_ret;
    //-----------------------------------------------------------------------------------
}




//Evaluation of the minmax objective
//@param[in] im: IM 
//@param[in] budget: the budget constraint
//@param The minmax score objective
double findMinMaxScore(IM *im, int budget){
    //-----------------------------------------------------------------------------------
    vector<double> scores, best_scores;
    vector<int> rewards = *(im->mdp->Reward);
    for(int i=0; i<im->mdp->P; i++){
        pair<vector<int>, double> non_state_score = findNonOptimal(im, budget, i, false);
        scores.push_back(im->F(i)-non_state_score.second);
        pair<vector<int>, double> state_score = findOptimal(im, budget, i, false);
        best_scores.push_back(state_score.second - non_state_score.second);
    }
    double minMaxScore = std::numeric_limits<double>::max();
    for(int i=0; i<im->mdp->P; i++){
        if((scores[i]/best_scores[i]) < minMaxScore) minMaxScore = (scores[i]/best_scores[i]);
    }
    *(im->mdp->Reward) = rewards;
    return minMaxScore;
    //-----------------------------------------------------------------------------------
}





//Random Algorithm implementation
//@param[in] im: IM 
//@param[in] budget: budget constraint
//@param[in] alpha: parameter alpha for the bicriteria approximation
//@return pair<vector<int>, double> pair of solution with the influence score
pair<vector<int>, double> randomSelection(IM *im, int budget, double alpha){
    im->emptyRewardStates();
    srand(time(NULL));
    vector<int> rewards;
    for(int i=0; i<im->mdp->N; i++) rewards.push_back(i);
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(rewards.begin(), rewards.end(), g);
    int current_cost = 0;
    for(int i=0; i<im->mdp->N; i++){
        int cost = im->mdp->findCost(rewards[i]);
        if((cost + current_cost) > alpha*budget) continue;
        im->setReward(rewards[i], 1);
        current_cost += cost;
    }
    return make_pair(*(im->mdp->Reward), findMinMaxScore(im, budget));
}



//All Greedy Algorithm implementation
//@param[in] im: IM 
//@param[in] budget: budget constraint
//@param[in] alpha: parameter alpha for the bicriteria approximation
//@return pair<vector<int>, double> pair of solution with the influence score
pair<vector<int>, double> allGreedy(IM *im, int budget, double alpha){
    vector<pair<vector<int>, double>> optimal_rewards_with_scores;
    vector<double> non_optimal_scores;
    for(int i=0; i<im->mdp->P; i++){  
        optimal_rewards_with_scores.push_back(findOptimal(im, budget, i, true));
        //cout<<"COST -----****----> "<<im->mdp->findCost(&(optimal_rewards_with_scores.at(optimal_rewards_with_scores.size()-1).first))<<" - "<<budget<<endl;
        non_optimal_scores.push_back(findNonOptimal(im, budget, i, false).second);
    }
    vector<pair<vector<int>, double>> alpha_optimal_rewards_with_scores;
    if(alpha==1){
        alpha_optimal_rewards_with_scores = optimal_rewards_with_scores;
    }
    else{
        for(int i=0; i<im->mdp->P; i++){  
            alpha_optimal_rewards_with_scores.push_back(findOptimal(im, alpha*budget, i, true));
        }
    }
    im->emptyRewardStates();
    pair<int, double> best_result;
    best_result.second = -1;
    for(int i=0; i<im->mdp->P; i++){
        im->setRewards(alpha_optimal_rewards_with_scores[i].first);
        double min_score = 1000000000000;
        for(int j=0; j<im->mdp->P; j++){
            double score_div_best = (im->F(j)-non_optimal_scores[i])/(optimal_rewards_with_scores[j].second-non_optimal_scores[i]);
            if(score_div_best <= min_score){
                min_score = score_div_best;
            }
        }
        if(min_score >= best_result.second){
            best_result.first = i;
            best_result.second = min_score;
        }
    }
    im->setRewards(alpha_optimal_rewards_with_scores[best_result.first].first);
    double minMaxScore = findMinMaxScore(im, budget);
    cout<<"All Greedy: "<<minMaxScore<<" ===== "<<best_result.second<<endl;
    return make_pair(*(im->mdp->Reward), best_result.second);
}



//Single Greedy Algorithm implementation
//@param[in] im: IM 
//@param[in] budget: budget constraint
//@param[in] alpha: parameter alpha for the bicriteria approximation
//@return pair<vector<int>, double> pair of solution with the influence score
pair<vector<int>, double> singleGreedy(IM *im, int budget, double alpha){
    im->emptyRewardStates();
    vector<double> optimal_scores;
    vector<double> non_optimal_scores;
    vector<double> prev_f;
    for(int i=0; i<im->mdp->P; i++){  
        optimal_scores.push_back(findOptimal(im, budget, i, false).second);
        non_optimal_scores.push_back(findNonOptimal(im, budget, i, false).second);
        prev_f.push_back(0);
        //cout<<"OPTIMAL: "<<optimal_rewards_with_scores[i].second<<endl;
    }
    unordered_set<int> V_minus_S;
    for(int i=0; i<im->mdp->N; i++) V_minus_S.insert(i);
    int current_cost = 0;
    double best_score = 0.0;
    
    while(V_minus_S.size()>0){
        double current_score = 0.0;
        int top_s = -1;
        for (auto s = V_minus_S.begin(); s != V_minus_S.end();){
            float cost_s = im->mdp->findCost(*s);
            if((current_cost+cost_s)>alpha*budget){
                V_minus_S.erase(s++);
                continue;
            }
            im->setReward(*s, 1);
            double min_score = 100000000000;
            for(int i=0; i<im->mdp->P; i++){ 
                double score_for_i = (im->F(i)-non_optimal_scores[i])/(optimal_scores[i]-non_optimal_scores[i]);
                if(score_for_i < min_score){
                    //cout<<"min_score for i "<<i<<"| "<<score_for_i<<endl;
                    min_score = score_for_i;
                }
                if(im->F(i)<(prev_f[i]-0.0001)){
                    //cout<<i<<" "<<im->F(i)<<" "<<prev_f[i]<<"-------"<<endl;
                }
            }
            //cout<<"-----> "<<min_score/cost_s<<" "<<*s<<endl;


            im->setReward(*s, 0);
            if(min_score/cost_s >= current_score){ //cost_s
                current_score = min_score/cost_s;
                //cout<<min_score/cost_s<<" for "<<*s<<endl;
                best_score = min_score;
                top_s = *s;
            }
            ++s;
        }
        if(top_s<0) break;
        //cout<<"Inserting tops "<<top_s<<" ***** ";
        //cout<<best_score<<" | "<<im->mdp->findCost(top_s)<<" | "<<best_score/im->mdp->findCost(top_s)<<" | "<<top_s<<endl;
        im->setReward(top_s, 1);
        current_cost += im->mdp->findCost(top_s);
        V_minus_S.erase(top_s);
        for(int i=0; i<im->mdp->P; i++){ 
            prev_f[i] = im->F(i);
        }
    }
    double minMaxScore = findMinMaxScore(im, budget);
    cout<<"Single Greedy: "<<minMaxScore<<" <=====> "<<best_score<<endl;
    return make_pair(*(im->mdp->Reward), best_score);
}



//Init a priority queue for the Myopic method
//@param[in] im: IM 
//@param[in] budget: budget constraint
//@param[in] policy: policy
//@param[in,out] V_minus_S: states that not included in the solution (or exceed the constraint)
//@param[in,out] current_cost: current cost
//@param[in,out] pq: the priority queue
//@param[in] alpha: parameter alpha for the bicriteria approximation
void single_init_priority_queue(IM *im, 
                                int budget, 
                                int policy, 
                                unordered_set<int> *V_minus_S, 
                                int *current_cost, 
                                priority_queue<pair<double, int> > *pq,
                                double alpha){
    //-----------------------------------------------------------------------------------
    for (auto s = V_minus_S->begin(); s != V_minus_S->end(); ){
        int cost_s = im->mdp->findCost(*s);
        if ((cost_s + *current_cost) > alpha*budget){
            V_minus_S->erase(s++);
            continue;
        }
        else{
            im->setReward(*s, 1);
            pq->push(make_pair(im->F(policy)/cost_s, *s)); //Both are ok marginal/cost_s and score/cost_s
            im->setReward(*s, 0);
            ++s;
        }
    }  
}



//Step of Myopic method
//@param[in] im: IM 
//@param[in] budget: budget constraint
//@param[in] policy: policy
//@param[in,out] V_minus_S: states that not included in the solution (or exceed the constraint)
//@param[in,out] current_cost: current cost
//@param[in,out] pq: the priority queue
//@param[in] alpha: parameter alpha for the bicriteria approximation
void single_lazy_step(IM *im, 
                      int budget, 
                      int policy, 
                      unordered_set<int> *V_minus_S, 
                      int *current_cost, 
                      priority_queue<pair<double, int> > *pq,
                      double alpha){
    //-----------------------------------------------------------------------------------
    unordered_map<int, double> checked;
    //cout<<"--> "<<pq->size()<<endl;
    while(pq->size()>0){
        int top_s = pq->top().second;
        pq->pop();
        if(V_minus_S->find(top_s) == V_minus_S->end()){
            continue;
        }
        float cost_s = im->mdp->findCost(top_s);
        if ((cost_s + *current_cost) > alpha*budget){
            V_minus_S->erase(top_s);
            continue;
        }
        if(checked.find(top_s) != checked.end()){
            pq->push(make_pair(10000, top_s));
            return;
        }
        im->setReward(top_s, 1);
        double score = im->F(policy);
        pq->push(make_pair(score/cost_s, top_s));
        im->setReward(top_s, 0);
        checked[top_s] = score; //Score after including node
    }
    //-----------------------------------------------------------------------------------
}



//Myopic Algorithm implementation
//@param[in] im: IM 
//@param[in] budget: budget constraint
//@param[in] alpha: parameter alpha for the bicriteria approximation
//@return pair<vector<int>, double> pair of solution with the influence score
pair<vector<int>, double> single_lazyGreedy(IM *im, int budget, double alpha){
    //-----------------------------------------------------------------------------------
    im->emptyRewardStates();
    vector<double> optimal_scores;
    vector<double> non_optimal_scores;
    vector<priority_queue<pair<double, int> >*> pqs;
    unordered_set<int> *V_minus_S = new unordered_set<int>();
    for(int i=0; i<im->mdp->N; i++) V_minus_S->insert(i);
    int current_cost = 0;
    double best_score = 0.0;
    //-----------------------------------------------------------------------------------
    for(int i=0; i<im->mdp->P; i++){  
        optimal_scores.push_back(findOptimal(im, budget, i, false).second);
        non_optimal_scores.push_back(findNonOptimal(im, budget, i, false).second);
        pqs.push_back(new priority_queue<pair<double, int> >());
    }
    int i = 0;
    for(auto pq: pqs){  
        single_init_priority_queue(im, budget, i, V_minus_S, &current_cost, pq, alpha);
        i++;
    }

    while(V_minus_S->size()>0){
        //cout<<"len "<<V_minus_S->size()<<endl;
        double current_score = 0.0;
        int top_s = -1;
        for(auto pq: pqs){  
            int s = pq->top().second;
            int cost_s = im->mdp->findCost(s);
            im->setReward(s, 1);
            double min_score = 100000000000;
            for(int i=0; i<im->mdp->P; i++){ 
                double score_for_i = (im->F(i)-non_optimal_scores[i])/(optimal_scores[i]-non_optimal_scores[i]);
                if(score_for_i <= min_score){
                    min_score = score_for_i;
                }
            }
            im->setReward(s, 0);
            if(min_score/cost_s >= current_score){ //cost_s
                current_score = min_score/cost_s;
                best_score = min_score;
                top_s = s;
            }
        }

        if(top_s<0) break;
        im->setReward(top_s, 1);
        current_cost += im->mdp->findCost(top_s);
        //cout<<"Inserting "<<top_s<<" new cost "<<current_cost<<endl;
        V_minus_S->erase(top_s);

        int policy = 0;
        for(auto pq: pqs){
            single_lazy_step(im, budget, policy, V_minus_S, &current_cost, pq, alpha);
            policy++;
        }





    }


    //-----------------------------------------------------------------------------------

    V_minus_S->clear(); delete V_minus_S;
    for(auto pq: pqs){  
        *pq = priority_queue<pair<double, int> >(); delete pq;
    }
    double minMaxScore = findMinMaxScore(im, budget);
    cout<<"Myopic: "<<minMaxScore<<" <=====> "<<best_score<<endl;
    return make_pair(*(im->mdp->Reward), best_score);
    //-----------------------------------------------------------------------------------
}




//Init a priority queue for the lazy greedy method used in the Saturate Greedy
//@param[in] im: IM 
//@param[in,out] current_cost: current cost
//@param[in,out] current_score: current score
//@param[in,out] pq: the priority queue
//@param[in] c: the value c for setting f' = min(f,c)
//@param[in] optimal_scores: optimal (best) scores for each policy 
//@param[in] non_optimal_scores: non-optimal (worst) scores for each policy [the default values for these scores are 0's]
void init_priority_queue_GPC(IM *im, int *current_cost, double *current_score, 
                             priority_queue<pair<double, int> > *pq, 
                             double c, vector<double>* optimal_scores, vector<double>* non_optimal_scores){
    //-----------------------------------------------------------------------------------
    for(int s=0; s<im->mdp->N; s++){
        int cost_s = im->mdp->findCost(s);
        double score = 0.0;
        im->setReward(s, 1);
        for(int i=0; i<im->mdp->P; i++){
            double score_for_i = (im->F(i)-(*non_optimal_scores)[i])/((*optimal_scores)[i] - (*non_optimal_scores)[i]);
            if(c <= score_for_i) score += c;
            else score += score_for_i;
        }
        im->setReward(s, 0);
        pq->push(make_pair(score/cost_s, s));
    }
    //-----------------------------------------------------------------------------------
    if(pq->empty()) return;
    im->setReward(pq->top().second, 1);

    int cost_top_s = im->mdp->findCost(pq->top().second); 
    *current_cost += cost_top_s;

    double ret_score = pq->top().first * cost_top_s; //This is due to the division by cost
    *current_score = ret_score;                      //Update current score

    pq->pop();
    //-----------------------------------------------------------------------------------
}



//Lazy greedy step of Saturate Greedy
//@param[in] im: IM 
//@param[in,out] current_cost: current cost
//@param[in,out] current_score: current score
//@param[in,out] pq: the priority queue
//@param[in] c: the value c for setting f' = min(f,c)
//@param[in] optimal_scores: optimal (best) scores for each policy 
//@param[in] non_optimal_scores: non-optimal (worst) scores for each policy [the default values for these scores are 0's]
void lazy_step_GPC(IM *im, int *current_cost, double *current_score, 
                   priority_queue<pair<double, int> > *pq, 
                   double c, vector<double>* optimal_scores, vector<double>* non_optimal_scores){
    //-----------------------------------------------------------------------------------
    unordered_map<int, double> checked;
    while(!pq->empty()){
        int top_s = pq->top().second;
        pq->pop();
        //if(V_minus_S->find(top_s) == V_minus_S->end()) continue;
        int cost_s = im->mdp->findCost(top_s);
        if(checked.find(top_s) != checked.end()){
            im->setReward(top_s, 1);
            *current_cost += cost_s;
            *current_score = checked[top_s];
            return;
        }
        im->setReward(top_s, 1);
        double score = 0.0;
        for(int i=0; i<im->mdp->P; i++){
            double score_for_i = (im->F(i)-(*non_optimal_scores)[i])/((*optimal_scores)[i]-(*non_optimal_scores)[i]);
            if(c <= score_for_i) score += c;
            else score += score_for_i;
        }
        pq->push(make_pair(score/cost_s, top_s));
        im->setReward(top_s, 0);
        checked[top_s] = score;
    }
    //-----------------------------------------------------------------------------------
}



//Lazy greedy of Saturate Greedy
//@param[in] im: IM 
//@param[in] c: the value c for setting f' = min(f,c)
//@param[in] optimal_scores: optimal (best) scores for each policy 
//@param[in] non_optimal_scores: non-optimal (worst) scores for each policy [the default values for these scores are 0's]
int lazyGPC(IM *im, double c, vector<double>* optimal_scores, vector<double>* non_optimal_scores){
    //-----------------------------------------------------------------------------------
    im->emptyRewardStates();
    int *current_cost = new int();
    double *current_score = new double();
    *current_cost = 0.0;
    *current_score = 0.0;
    priority_queue<pair<double, int> > *pq;
    pq = new priority_queue<pair<double, int> >();
    init_priority_queue_GPC(im, current_cost, current_score, pq, c, optimal_scores, non_optimal_scores);
    while(((*current_score)/(im->mdp->P) < (c - 1/100000))){
        lazy_step_GPC(im, current_cost, current_score, pq, c, optimal_scores, non_optimal_scores);
        if(pq->empty()) break;
        //cout<<(*current_score)/(im->mdp->P)<<" _ "<<*current_cost<<" "<<pq->empty()<<endl;
    }
    double cost = *current_cost;
    delete current_cost; delete current_score;
    *pq = priority_queue<pair<double, int> >(); delete pq;
    return cost;
    //-----------------------------------------------------------------------------------
}



//Saturate Greedy implementation
//@param[in] im: IM 
//@param[in] budget: budget constraint
//@param[in] alpha: parameter alpha for the bicriteria approximation
//@return pair<vector<int>, double> pair of solution with the influence score
pair<vector<int>, double> Saturate(IM *im, int budget, double alpha){
    im->emptyRewardStates();
    vector<double>* optimal_scores = new vector<double>();
    vector<double>* non_optimal_scores = new vector<double>();
    for(int i=0; i<im->mdp->P; i++){
        optimal_scores->push_back(findOptimal(im, budget, i, false).second);
        non_optimal_scores->push_back(findNonOptimal(im, budget, i, false).second);
    }
    double cmin = 0;
    int current_cost = 0;
    double cmax = 2;
    int no_of_policies = im->mdp->P;
    int big_number = 1000;
    while((cmax-cmin)>=1.0/(no_of_policies*big_number)){
        double c = (cmax+cmin)/2.0;
        //cout<<"c -----> "<<c<<endl;
        current_cost = lazyGPC(im, c, optimal_scores, non_optimal_scores);
        //cout<<"The real is: "<<findMinMaxScore(im, budget)<<endl;
        //cout<<current_cost<<" "<<alpha * budget<<" "<<alpha<<endl;
        if(current_cost > alpha * budget) cmax = c;
        else cmin = c;
    }
        if(current_cost > alpha * budget){
            //cout<<"INSIIIIIIIDE"<<endl;
            //cout<<findMinMaxScore(im, budget)<< " vs "<<current_cost<<endl;
            current_cost = lazyGPC(im, cmin, optimal_scores, non_optimal_scores);
        }
        //cout<<findMinMaxScore(im, budget)<< " vs "<<current_cost<<endl;
        //-----------------------------------------------------------------------------------
        //double current_score = -1;
        //current_score = -1;
        //-----------------------------------------------------------------------------------
        vector<priority_queue<pair<double, int> >*> pqs;
        unordered_set<int> *V_minus_S = new unordered_set<int>();
        for(int i=0; i<im->mdp->N; i++){
            if((*(im->mdp->Reward))[i] == 0)
                V_minus_S->insert(i);
        }
        //-----------------------------------------------------------------------------------
        for(int i=0; i<im->mdp->P; i++){  
            optimal_scores->push_back(findOptimal(im, budget, i, false).second);
            non_optimal_scores->push_back(findNonOptimal(im, budget, i, false).second);
            pqs.push_back(new priority_queue<pair<double, int> >());
        }
        int i = 0;
        for(auto pq: pqs){  
            single_init_priority_queue(im, budget, i, V_minus_S, &current_cost, pq, alpha);
            i++;
        }
        
        while(V_minus_S->size()>0){
            //cout<<"len "<<V_minus_S->size()<<endl;
            int top_s = -1;
            double current_score = -1;
            for(auto pq: pqs){  
                int s = pq->top().second;
                int cost_s = im->mdp->findCost(s);
                im->setReward(s, 1);
                double min_score = 100000000000;
                for(int i=0; i<im->mdp->P; i++){ 
                    double score_for_i = (im->F(i)-(*non_optimal_scores)[i])/((*optimal_scores)[i]-(*non_optimal_scores)[i]);
                    if(score_for_i <= min_score){
                        min_score = score_for_i;
                    }
                }
                im->setReward(s, 0);
                //cout<<min_score/cost_s <<" __ "<<current_score<<endl;
                if(min_score/cost_s >= current_score){ //cost_s
                    current_score = min_score/cost_s;
                    top_s = s;
                }
            }

            if(top_s<0) break;
            im->setReward(top_s, 1);
            current_cost += im->mdp->findCost(top_s);
            //cout<<"Inserting "<<top_s<<" new cost "<<current_cost<<endl;
            V_minus_S->erase(top_s);

            int policy = 0;
            for(auto pq: pqs){
                single_lazy_step(im, budget, policy, V_minus_S, &current_cost, pq, alpha);
                policy++;
            }
        }
        //-----------------------------------------------------------------------------------
    
    double minMaxScore = findMinMaxScore(im, budget);
    optimal_scores->clear(); delete optimal_scores;
    non_optimal_scores->clear(); delete non_optimal_scores;
    //cout<<"*****************************-------------------------> "<<current_cost<<endl;
    return make_pair(*(im->mdp->Reward), minMaxScore);
}
