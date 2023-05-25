#pragma once 
using namespace std;
#include<stdio.h>
#include <iostream>
#include <vector>
#include <iostream>
#include <vector>
#include <fstream>
#include <tuple>
#include <string>
#include <sstream>
#include <unordered_set>

typedef pair<int, double> Edge;
typedef tuple<int, int, int> triple;

void printVector(const vector<Edge> &edges){
    for(auto edge: edges){
        printf("(%d|%f) ", edge.first, edge.second);
    }
    printf("\n");
}

void printVector(const vector<double> &steps){
    for(int k=0; k<steps.size(); k++){
        printf("(%d|%f) ", k, steps[k]);
    }
    printf("\n");
}

void printVector(const vector<int> &rewards){
    for(int i=0; i<rewards.size(); i++){
        printf("(%d|%d) ", i, rewards[i]);
    }
    printf("\n");
}

void printSolution(const vector<int> &rewards){
    for(int i=0; i<rewards.size(); i++){
        if(rewards[i] == 1){
            printf("%d | ", i);
        }
    }
    printf("\n");
}

int sizeOfSolution(const vector<int> &rewards){
    int counter = 0;
    for(int i=0; i<rewards.size(); i++){
        if(rewards[i] == 1){
            counter++;
        }
    }
    return counter;
}

void printSet(unordered_set<int> &S){
    for(auto s: S){
        printf("(%d) ", s);
    }
    printf("\n");
}

int read_N(const string& file){
    string first_line;
    ifstream infile(file);
    getline(infile, first_line);
    infile.close();
    return stoi(first_line);
}

triple read_line(const string &line) {
    char delimiter = ' ';
    istringstream ss(line);
    string token;
    vector<int> numbers;
    while (getline(ss, token, delimiter)) {
        numbers.push_back(stoi(token));
    }

    return {numbers[0], numbers[1], numbers[2]};
}