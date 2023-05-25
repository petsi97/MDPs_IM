#include <iostream>
#include "mc.h"
#include "mdp.h"
#include "rim.h"
#include "algorithms.h"
#include <chrono>
#include <fstream>
#include <filesystem>
namespace fs = std::filesystem;


void init_results(string foldername){
    fs::create_directories(foldername);
    vector<string> filenames = {foldername+"/saturate.txt", 
                                //foldername+"/singleGreedy.txt", 
                                foldername+"/allGreedy.txt", 
                                foldername+"/myopic.txt"}; 
    for(auto filename: filenames){
        ofstream foutput; 
        foutput.open (filename,ios::app); 
        foutput<<"N P S budget alpha Score Cost Time Size PreTime"<<endl;
        foutput.close(); 
    }
}


void write_new_line(string filename, vector<double> new_line){
    ofstream foutput; 
    foutput.open (filename,ios::app); 
    string line;
    for(int i=0; i<new_line.size(); i++){
        if(i==0){
            line = to_string(int(new_line[i]));
        }
        else if((i==3)||(i==4)||(i==5)||(i==7)||(i==9)){
            line = line + " " + to_string(new_line[i]);
        }
        else{
            line = line + " " + to_string(int(new_line[i]));
        }
    }
    foutput<<line<<endl;
    foutput.close(); 
}


vector<MC> make_im(string read_folder, string name, int P, int step, bool uniform_steps, string name_end){
    vector<MC> MCs;
    for(int i = 0; i<P; i++){
        string toname =read_folder+"/"+name+name_end+"/"+name+name_end+"-"+to_string(i)+".txt";
        MC *mc = new MC(toname);
        mc->normalizeEdgeWeight();
        if(uniform_steps) mc->setUniformSteps(step);
        else mc->setNormalDistribution(step);
        MCs.push_back(*mc);
    }
    return MCs;
}


void experiment1(string foldername, vector<int> Ns, vector<int> Ps, vector<int> Steps, vector<double> budgets,
                 vector<double> alphas, bool use_IC, bool uniform_steps, string read_folder, string name_end, 
                 bool inres){ 
    std::chrono::steady_clock::time_point begin, end;
    IM *im;
    vector<MC> mcs;
    pair<vector<int>, double> sol;
    if(inres) init_results(foldername);
    vector<double> res;
    for(int i=0; i<10; i++){
        res.push_back(0);
    }
    //"0 1 2   3      4    5     6     7    8      9"
    //"N P S budget alpha Score Cost Time Size PreTime";
    int diftime, pretime;
    for(auto N: Ns){
        res[0] = N;
        for(auto P: Ps){
            res[1] = P;
            for(auto step: Steps){
                res[2] = step;
                for(auto alpha: alphas){
                    res[4] = alpha;
                    mcs = make_im(read_folder, to_string(N), P, step, uniform_steps, name_end);
                    MDP mdp(&mcs);
                    cout<<read_folder+"/"+to_string(N)+name_end+"/costs_"+to_string(N)+name_end+".txt"<<endl;
                    mdp.read_costs(read_folder+"/"+to_string(N)+name_end+"/costs_"+to_string(N)+name_end+".txt");
                    begin = std::chrono::steady_clock::now();
                    if(use_IC) im = new IM_IC(&mdp);
                    else im = new IM_HT(&mdp);
                    end = std::chrono::steady_clock::now();
                    pretime =  std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count(); 
                    res[9] = pretime;
                    for(auto b: budgets){
                        int budget = (int)(mdp.total_budget * b);
                        budget = (int)im->mdp->N*b;
                        //int budget = (int)im->mdp->N*b;
                        res[3] = b;
                        cout<<"Budget is "<<budget<<endl;
                        cout<<"-----------------"<<endl;
                        begin = std::chrono::steady_clock::now();
                        sol = Saturate(im, budget, alpha);
                        end = std::chrono::steady_clock::now();
                        diftime = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
                        //cout<<"Saturate | solution: ";printSolution(sol.first);
                        cout<<"Saturate | score: "<<sol.second<<endl; res[5] = sol.second;
                        cout<<"Saturate | cost: "<<im->mdp->findCost(&sol.first)<<endl; res[6] = im->mdp->findCost(&sol.first);
                        cout<<"Saturate | time: "<<diftime<<" [msecs]"<<endl; res[7] = diftime; 
                        cout<<"Saturate | size: "<<sizeOfSolution(sol.first)<<endl; res[8] = sizeOfSolution(sol.first);
                        write_new_line(foldername+"/saturate.txt", res);
                        cout<<"-----------------"<<endl;

                        cout<<"-----------------"<<endl;
                        begin = std::chrono::steady_clock::now();
                        sol = allGreedy(im, budget, alpha);
                        end = std::chrono::steady_clock::now();
                        diftime = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
                        //cout<<"AllGreedy | solution: ";printSolution(sol.first);
                        cout<<"AllGreedy | score: "<<sol.second<<endl; res[5] = sol.second;
                        cout<<"AllGreedy | cost: "<<im->mdp->findCost(&sol.first)<<endl; res[6] = budget;
                        cout<<"AllGreedy | time: "<<diftime<<" [msecs]"<<endl; res[7] = diftime; 
                        cout<<"AllGreedy | size: "<<sizeOfSolution(sol.first)<<endl; res[8] = sizeOfSolution(sol.first);
                        write_new_line(foldername+"/allGreedy.txt", res);
                        cout<<"-----------------"<<endl;
                        
                        cout<<"-----------------"<<endl;
                        begin = std::chrono::steady_clock::now();
                        sol = single_lazyGreedy(im, budget, alpha);
                        end = std::chrono::steady_clock::now();
                        diftime = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
                        //cout<<"Myopic | solution: ";printSolution(sol.first);
                        cout<<"Myopic | score: "<<sol.second<<endl; res[5] = sol.second;
                        cout<<"Myopic | cost: "<<im->mdp->findCost(&sol.first)<<endl; res[6] = budget;
                        cout<<"Myopic | time: "<<diftime<<" [msecs]"<<endl; res[7] = diftime; 
                        cout<<"Myopic | size: "<<sizeOfSolution(sol.first)<<endl; res[8] = sizeOfSolution(sol.first);
                        write_new_line(foldername+"/myopic.txt", res);
                        cout<<"-----------------"<<endl;   
                    }
                    delete im;
                }
            }
        }
    }

    

    

}


void N_experiments(string foldername, string read_folder, string name_end, bool use_IC, bool uniform_steps){
    vector<int> Ns = {2500, 5000, 7500, 10000, 12500};
    vector<int> Ps = {10};
    vector<int> Steps = {6};
    vector<double> budgets = {0.25};
    vector<double> alphas = {1.0};
    experiment1(foldername, Ns, Ps, Steps, budgets, alphas, use_IC, uniform_steps, read_folder, name_end, true);
}


void pi_experiments(string foldername, string read_folder, string name_end, bool use_IC, bool uniform_steps){
    vector<int> Ns = {10000};
    vector<int> Ps = {1, 3, 5, 10, 15, 20};
    vector<int> Steps = {6};
    vector<double> budgets = {0.25};
    vector<double> alphas = {1.0};
    experiment1(foldername, Ns, Ps, Steps, budgets, alphas, use_IC, uniform_steps, read_folder, name_end, true);
}


void budget_experiments(string foldername, string read_folder, string name_end, bool use_IC, bool uniform_steps){
    vector<int> Ns = {10000};
    vector<int> Ps  {10};
    vector<int> Steps = {6};
    vector<double> budgets = {0.1, 0.25, 0.5, 0.75};
    vector<double> alphas = {1.0};
    experiment1(foldername, Ns, Ps, Steps, budgets, alphas, use_IC, uniform_steps, read_folder, name_end, true);
}


void steps_experiments(string foldername, string read_folder, string name_end, bool use_IC, bool uniform_steps){
    vector<int> Ns = {10000};
    vector<int> Ps =  {10};
    vector<int> Steps = {2, 4, 6, 8, 10};
    vector<double> budgets = {0.25};
    vector<double> alphas = {1.0};
    experiment1(foldername, Ns, Ps, Steps, budgets, alphas, use_IC, uniform_steps, read_folder, name_end, true);
}


void alphas_experiments(string foldername, string read_folder, string name_end, bool use_IC, bool uniform_steps){
    vector<int> Ns = {10000};
    vector<int> Ps =  {10};
    vector<int> Steps = {6};
    vector<double> budgets = {0.25};
    vector<double> alphas = {1.0, 1.25, 1.5, 1.75, 2};
    experiment1(foldername, Ns, Ps, Steps, budgets, alphas, use_IC, uniform_steps, read_folder, name_end, true);
}


void avg_deg_experiments(string foldername, string read_folder, string name_end){
    vector<int> Ns = {10000};
    vector<int> Ps =  {10};
    vector<int> Steps = {6};
    vector<double> budgets = {0.25};
    vector<double> alphas = {1.0};
    cout<<foldername<<endl<<read_folder<<endl<<name_end<<endl;
    cout<<read_folder+"/erdos_avg_"+to_string(3)<<endl;
    cout<<"_"+name_end+"_"+to_string(3)<<endl;
    //return;
    for(int iter = 0; iter<10; iter++){
        bool inres = false;
        if(iter == 0) inres = true;
        for(auto avg:{3,6,9,12})
            experiment1(foldername+"/"+to_string(avg), Ns, Ps, Steps, budgets, alphas, true, true, 
            read_folder+"/erdos_avg_"+to_string(avg)+"_id_"+to_string(iter), 
            "_"+name_end+"_avg_"+to_string(avg)+"_id_"+to_string(iter), inres);
    }
}


void distr_experiments(string foldername, string read_folder, string name_end){
    vector<int> Ns = {10000};
    vector<int> Ps =  {10};
    vector<int> Steps = {6};
    vector<double> budgets = {0.25};
    vector<double> alphas = {1.0};
    cout<<foldername<<endl<<read_folder<<endl<<name_end<<endl;
    cout<<read_folder+"/scale_beta_"+to_string(8)<<endl;
    cout<<"_"+name_end+"_"+to_string(8)<<endl;
    //return;
    for(int iter = 0; iter<10; iter++){
        bool inres = false;
        if(iter == 0) inres = true;
        for(auto distr:{6,7,8,9})
            experiment1(foldername+"/"+to_string(distr), Ns, Ps, Steps, budgets, alphas, true, true, 
            read_folder+"/scale_beta_"+to_string(distr)+"_id_"+to_string(iter), 
            "_"+name_end+"_beta_"+to_string(distr)+"_id_"+to_string(iter), inres);
    }
}



int main()
{
    string save_foler ="./24_05_2023/";
    string read_folder = "./Datasets";
    vector<string> parameters = {"N", "Pi", "Steps", "Budget", "Alpha"};
    vector<string> graphnames = {"erdos", "scale"};
    vector<bool> use_ICs = {true};
    vector<bool> uniform_steps = {true};

    for(int id=0; id<10; id++){
        for(auto parameter: parameters){
            for(auto graphname: graphnames){
                if(graphname == "erdos") graphname += "_avg_6_id_"+to_string(id);
                if(graphname == "scale") graphname += "_beta_8_id_"+to_string(id);
                for(auto use_IC: use_ICs){
                    for(auto uniform_step: uniform_steps){
                        string foldername =  save_foler+graphname;
                        foldername+="/"+parameter;
                        string name_end = "_"+graphname;
                        //N
                        if(parameter == "N") N_experiments(foldername, read_folder+"/"+graphname, name_end, use_IC, uniform_step);
                        //Pi
                        if(parameter == "Pi") pi_experiments(foldername, read_folder+"/"+graphname, name_end, use_IC, uniform_step);
                        //Steps
                        if(parameter == "Steps") steps_experiments(foldername, read_folder+"/"+graphname, name_end, use_IC, uniform_step);
                        //Budget
                        if(parameter == "Budget") budget_experiments(foldername, read_folder+"/"+graphname, name_end, use_IC, uniform_step);
                        //Alpha
                        if(parameter == "Alpha") alphas_experiments(foldername, read_folder+"/"+graphname, name_end, use_IC, uniform_step);
                    }
                }
            }
        }
    }

    //Avg - Scale
    string graphname = "scale";
    distr_experiments(save_foler+graphname+"/distr", read_folder, graphname);

    //Beta - Erdos
    graphname = "erdos";
    avg_deg_experiments(save_foler+graphname+"/avg", read_folder, graphname);

    return 0;
}