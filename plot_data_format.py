import pandas as pd
import matplotlib.pyplot as plt 
import functools as ft
import sys
import os
from os import path

def make_csv():
    folder_name = ['24_05_2023']  
    dataset_type = ['erdos']
    iters = [str(i) for i in range(10)]
    parameters = ['Alpha', 'Budget', 'N', 'Pi', 'Steps']
    algorithms = ['allGreedy', 'myopic', 'saturate']
    
    for name in folder_name:
        for datatype in dataset_type:
            for param in parameters:
                for alg in algorithms:
                    dfs = []
                    for iter in iters:
                        if datatype == 'erdos': filename = f'./{name}/{datatype}_avg_6_id_{iter}/{param}/{alg}.txt'
                        else: filename = f'./{name}/{datatype}_beta_8_id_{iter}/{param}/{alg}.txt'
                        print(filename)
                        dfs.append(pd.read_csv(filename, sep=" "))
                    dfall = dfs[0]/len(iters)
                    for i in range(1,len(dfs)):
                        dfall = dfall.add(dfs[i]/len(iters), fill_value=0)
                    print('_________Parameter ', param,'__________')
                    print(dfall)
                    print(f'./{name}/{datatype}/{param}')
                    print(f'./{name}/{datatype}/{param}/{alg}.txt')
                    isExist = os.path.exists(f'./{name}/{datatype}/{param}')
                    if not isExist:
                        os.makedirs(f'./{name}/{datatype}/{param}')
                    with open(f'./{name}/{datatype}/{param}/{alg}.txt', 'w') as f:
                        f.write(dfall.to_string(header=True, index=False))

def make_csv_distribution_plots():
    folder_name = ['24_05_2023']  
    dataset_type = ['erdos']
    avgs = ['3', '6', '9', '12']
    parameters = ['avg']
    algorithms = ['allGreedy', 'myopic', 'saturate']
    
    for name in folder_name:
        for datatype in dataset_type:
            for param in parameters:
                for alg in algorithms:
                    dfs = []
                    for avg in avgs:
                        filename = f'./{name}/{datatype}/{param}/{avg}/{alg}.txt'
                        print(filename)
                        dfs.append(pd.read_csv(filename, sep=" "))
                        dfs[-1]['avg'] = int(avg)
                        dfs[-1] = dfs[-1].groupby(['N']).mean()
                    dfall = dfs[0]
                    for i in range(1,len(dfs)):
                        dfall = dfall.append(dfs[i], ignore_index=True)
                    print('_________Parameter ', param,'__________')
                    print(dfall)
                    isExist = os.path.exists(f'./{name}/{datatype}/{param}')
                    if not isExist:
                        os.makedirs(f'./{name}/{datatype}/{param}')
                    with open(f'./{name}/{datatype}/{param}/{alg}.txt', 'w') as f:
                        f.write(dfall.to_string(header=True, index=False))
    dataset_type = ['scale']
    avgs = ['6', '7', '8', '9']
    parameters = ['distr']                
    for name in folder_name:
        for datatype in dataset_type:
            for param in parameters:
                for alg in algorithms:
                    dfs = []
                    for avg in avgs:
                        filename = f'./{name}/{datatype}/{param}/{avg}/{alg}.txt'
                        print(filename)
                        dfs.append(pd.read_csv(filename, sep=" "))
                        dfs[-1]['avg'] = int(avg)
                        dfs[-1]['avg'] = int(avg)
                        dfs[-1] = dfs[-1].groupby(['N']).mean()
                        print(dfs[-1])
                    dfall = dfs[0]
                    for i in range(1,len(dfs)):
                        dfall = dfall.append(dfs[i], ignore_index=True)
                    print('_________Parameter ', param,'__________')
                    print(dfall)
                    isExist = os.path.exists(f'./{name}/{datatype}/{param}')
                    if not isExist:
                        os.makedirs(f'./{name}/{datatype}/{param}')
                    with open(f'./{name}/{datatype}/{param}/{alg}.txt', 'w') as f:
                        f.write(dfall.to_string(header=True, index=False))

if __name__ == "__main__":   
    make_csv()
    make_csv_distribution_plots()