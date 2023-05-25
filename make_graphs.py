import sys
from random import randint, shuffle
import os
from os import path
import numpy as np


#return random integer in range [0, max_val]
def random_int(max_val):
    return randint(0, max_val)


#sample edge-weigths from normal distribution
def sample_edge_weight_normal_distrib(G, edges, sigma = 0):
    edges = list(edges)
    for i in range(len(edges)):
        out_deg = G.out_degree(edges[i][0])
        edge_weight = np.random.normal(out_deg, 2*(out_deg)*sigma, 1)
        #print(edges[i][0], edges[i][1])
        edges[i] = (edges[i][0], edges[i][1], edge_weight[0])
    return edges
 
    
#make erdos-renyi 
def make_erdos_graphs(no_of_nodes, num_of_graphs, folder = 'data_new2', name = 'erdos', avg = 10):
    import networkx as nx
    prob = (avg*no_of_nodes)/(no_of_nodes*(no_of_nodes-1))
    G = nx.erdos_renyi_graph(no_of_nodes, prob, seed=None, directed=True)
    G = nx.DiGraph(G)
    print('Graph Characteristics: ', G.number_of_nodes(), G.number_of_edges(), G.number_of_edges()/G.number_of_nodes())
    name_graph = f'{name}'
    folder += "/"+name
    for i in range(num_of_graphs):        
        edges = sample_edge_weight_normal_distrib(G, nx.to_edgelist(G), (i+1)/num_of_graphs)
        full = None
        direc = fr"./{folder}/{no_of_nodes}_{name_graph}"
        isExist = os.path.exists(direc)
        if not isExist:
            os.makedirs(direc)
        
        # Allow up to 100 graphs
        for i in range(100):
            filename = f"{no_of_nodes}_{name_graph}-{i}"
            full = fr"./{folder}/{no_of_nodes}_{name_graph}/{filename}.txt"
            if not path.isfile(full):
                break

        print(f'Create file at {full}')
        with open(full, 'w') as file:
            file.write(f'{no_of_nodes}\n')
            
            init_probs = [0 for _ in range(no_of_nodes)]
            
            for edge in edges:
                if(edges[0] == edges[1]): continue
                if edge[2] < 0:
                    file.write(f'{edge[0]} {edge[1]} {0}\n')
                else: 
                    file.write(f'{edge[0]} {edge[1]} {1+edge[2]}\n')
                    init_probs[edge[0]] += edge[2]
                
            for index in range(no_of_nodes):
                file.write(f'{index} {index} {init_probs[index]}\n')
    
    costs = []
    for index in range(no_of_nodes):
        costs.append(1+G.in_degree(index))

    filename = f"costs_{no_of_nodes}_{name_graph}"
    full = fr"./{folder}/{no_of_nodes}_{name_graph}/{filename}.txt"
    if path.isfile(full):
        sys.exit()

    print(f'Create file at {full}')
    with open(full, 'w') as file:
        file.write(f'{no_of_nodes}\n')
        for i in range(no_of_nodes):
            file.write(f'{costs[i]}\n')


#make scale-free
def make_scale_graphs(no_of_nodes, num_of_graphs, folder = 'data_new2', name = 'scale', beta = 6):
    import networkx as nx
    G = nx.scale_free_graph(no_of_nodes, alpha=2*(1-beta/10)/3, beta=beta/10, gamma=(1-beta/10)/3, delta_in=0.7, delta_out=0.3)
    G = nx.DiGraph(G)
    print('Graph Characteristics: ', G.number_of_nodes(), G.number_of_edges(), G.number_of_edges()/G.number_of_nodes())
    name_graph = f'{name}'
    folder += "/"+name
    for i in range(num_of_graphs):        
        edges = sample_edge_weight_normal_distrib(G, nx.to_edgelist(G), (i+1)/num_of_graphs)
        full = None
        direc = fr"./{folder}/{no_of_nodes}_{name_graph}"
        isExist = os.path.exists(direc)
        if not isExist:
            os.makedirs(direc)
        
        # Allow up to 100 graphs
        for i in range(100):
            filename = f"{no_of_nodes}_{name_graph}-{i}"
            full = fr"./{folder}/{no_of_nodes}_{name_graph}/{filename}.txt"
            if not path.isfile(full):
                break

        print(f'Create file at {full}')
        with open(full, 'w') as file:
            file.write(f'{no_of_nodes}\n')
            
            init_probs = [0 for _ in range(no_of_nodes)]
            
            for edge in edges:
                if(edges[0] == edges[1]): continue
                if edge[2] < 0:
                    file.write(f'{edge[0]} {edge[1]} {0}\n')
                else: 
                    file.write(f'{edge[0]} {edge[1]} {edge[2]}\n')
                    init_probs[edge[0]] += edge[2]
                
            for index in range(no_of_nodes):
                file.write(f'{index} {index} {1+init_probs[index]}\n')
    
    costs = []
    for index in range(no_of_nodes):
        costs.append(1+G.in_degree(index))

    filename = f"costs_{no_of_nodes}_{name_graph}"
    full = fr"./{folder}/{no_of_nodes}_{name_graph}/{filename}.txt"
    if path.isfile(full):
        sys.exit()

    print(f'Create file at {full}')
    with open(full, 'w') as file:
        file.write(f'{no_of_nodes}\n')
        for i in range(no_of_nodes):
            file.write(f'{costs[i]}\n')


if __name__ == "__main__":    
    folder = 'Datasets'
    
    for t in range(10):
        for num_of_nodes in [2500, 5000, 7500, 10000, 12500]:
            for beta in [6, 7, 8, 9]:
                if beta == 8:
                    make_scale_graphs(num_of_nodes, 20, folder, 'scale_beta_'+str((int)(beta))+"_id_"+str(t), beta) 
                elif num_of_nodes == 10000:
                    make_scale_graphs(num_of_nodes, 20, folder, 'scale_beta_'+str((int)(beta))+"_id_"+str(t), beta)
                    
    for t in range(10):
        for num_of_nodes in [2500, 5000, 7500, 10000, 12500]:
            for avg in [3,6,9,12]:
                if avg == 6:
                    make_erdos_graphs(num_of_nodes, 20, folder, 'erdos_avg_'+str((int)(avg))+"_id_"+str(t), avg) 
                elif num_of_nodes == 10000:
                    make_erdos_graphs(num_of_nodes, 20, folder, 'erdos_avg_'+str((int)(avg))+"_id_"+str(t), avg) 

    
    
    
    
    