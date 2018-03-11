#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  10 10:26:12 2018

@author: xuecho
"""
import pandas as pd
from igraph import *
import networkx as nx

import pickle
import time
import numpy as np
import matplotlib.pyplot as plt
import collections
import pickle
from random import randint
from random import choice
from random import shuffle


########################################
#Functions
########################################\
def save_obj(obj, name):
    f = open(name + ".dat", "wb")
    pickle.dump(obj, f)
    f.close()

def load_obj(name):
    return pickle.load(open(name + ".dat", "rb"))


def plot_infected_per_time(nodesD):
    avg = {}
    var = {}
    for t in nodesD[0]:   # for all times
        #print t
        timeList = [d[t] for d in nodesD]
        avg[t] = np.mean(timeList)
        var[t] = np.std(timeList)

    #save_obj(avg, "average9")
    #save_obj(var, "var9")
    plot_Dict(avg, var ,"lightblue")
    #plot_Dict(var, "green")
    plt.show()
"""
def plot_Dict(d, err, mycolor):
    for t in d:
        #plt.scatter(t,d[t], color = mycolor, s=50)
        plt.errorbar(t, d[t], err[t],  fmt='o', color=mycolor)
"""      
def plot_Dict(d, err, mycolor):
    for t in d:
        #plt.scatter(t,d[t], color = mycolor, s=50)
        plt.plot(t, d[t], 'ro', markersize=0.4)
        plt.errorbar(t, d[t], yerr=err[t],  fmt='', ecolor=mycolor)

def plot_Dict_2(D,mycolor):
    plt.bar(range(len(D)), list(D.values()), align='center',color = mycolor)
    plt.xticks(range(len(D)), list(D.keys()))

  

def create_aggr_network(filename):
    G=nx.Graph()    
    data = pd.read_table(filename,header=0,delim_whitespace=True)
    extCol = np.array(data.iloc[:,[0,1]]).tolist()
    for row in extCol:
        vertex1 = int(row[0])
        vertex2 = int(row[1])
        G.add_edge(vertex1,vertex2)
    return G  


def create_edgelist_per_time(filename, timeD): #create edgelist 
    data = pd.read_table(filename,header=0,delim_whitespace=True)
    extCol = np.array(data.iloc[:,[0,1,2]]).tolist()
    for row in extCol:
        vertex1 = int(row[0])
        vertex2 = int(row[1])
        t = int(row[2])
        
        # Create edgelist per time - Adding one edge at a time   
        if t in timeD:
            timeD[t].append((vertex1, vertex2))
        else:
            timeD[t]=[]
            timeD[t].append((vertex1, vertex2))

def create_graph(edgeList):
    G=nx.Graph()    
    for edge in edgeList:
        G.add_edge(edge[0], edge[1])
    return G

def rewiring(G,filename): #if deg==0 positive else negative
    data = pd.read_table(filename,header=0,delim_whitespace=True)
    timet = np.array(data.iloc[:,2]).tolist()
    shuffle(timet);
    edgeList=[]
    edges = list(G.edges())
    for i in range(1,len(edges)+1):
        e1 = choice(edges)
        edgeList.append(e1);
    return create_graph(edgeList)

def degree_rewiring(G,deg): #if deg==0 positive else negative
    tup=[]
    edges = G.edges()
    e1 = choice(edges)
    e2 = choice(edges)
    while e1 == e2:
        e2 = choice(edges)
    tup.append((e1[0],G.degree(e1[0])))
    tup.append((e1[1],G.degree(e1[1])))
    tup.append((e2[0],G.degree(e2[0])))
    tup.append((e2[1],G.degree(e2[1])))

    tup.sort(key=lambda x:x[1])
    if deg==1:
        newEdge1 = (tup[0][0],tup[3][0])
        newEdge2 = (tup[1][0],tup[2][0])
    else:
        newEdge1 = (tup[0][0],tup[1][0])
        newEdge2 = (tup[2][0],tup[3][0])
    if newEdge1 not in G.edges() and newEdge2 not in G.edges():
        G.add_edge(*newEdge1)
        G.add_edge(*newEdge2)
        G.remove_edge(*e1)
        G.remove_edge(*e2)
    return G


def create_edgelist_reshuffle_time(filename): #create edgelist 
    timeD={} 
    data = pd.read_table(filename,header=0,delim_whitespace=True)
    extCol = np.array(data.iloc[:,[0,1,2]]).tolist()
    timet = np.array(data.iloc[:,2]).tolist()
    shuffle(timet);
    for i in range(len(extCol)):
        vertex1 = int(extCol[i][0])
        vertex2 = int(extCol[i][1])
        t = timet[i]
        # Create edgelist per time - Adding one edge at a time     
        if t in timeD:
            timeD[t].append((vertex1, vertex2))
        else:
            timeD[t]=[]
            timeD[t].append((vertex1, vertex2))
    return timeD


def create_edgelist_per_random_time(filename): #create edgelist 
    timeD={} 
    data = pd.read_table(filename,header=0,delim_whitespace=True)
    for row in data:
        vertex1 = int(row[0])
        vertex2 = int(row[1])
        t = randint(1,7375)
        # Create edgelist per time - Adding one edge at a time
        
        if t in timeD:
            timeD[t].append((vertex1, vertex2))
        else:
            timeD[t]=[]
            timeD[t].append((vertex1, vertex2))
    return timeD


def create_edgelist_from_graph_random_time(filename,G,timeD):
    data = pd.read_table(filename,header=0,delim_whitespace=True)
    timet = np.array(data.iloc[:,2]).tolist()
    shuffle(timet);
    edgeList = list(G.edges())
    
    for i in range(len(edgeList)):
        vertex1 = int(edgeList[i][0])
        vertex2 = int(edgeList[i][1])
        t = timet[i]
       # Create edgelist per time - Adding one edge at a time
        if t in timeD:
           timeD[t].append((vertex1, vertex2))
        else:
           timeD[t]=[]
           timeD[t].append((vertex1, vertex2))
    return timeD

def stimulation(G,timeD):
    nodesD = []
    total_number_of_nodes = G.number_of_nodes()
    for node in range(1, total_number_of_nodes+1):
        #for node in range
        print("node", node)
        allInf = [node] # List with all the currently infected nodes
        inf = {}  # Dictionary with the infected nodes per time
      
        
    #    for t in timeD: #Creating the network topology per time
        o = max(timeD)+1
        for t in range(1,o): #Creating the network topology per time
            #print("time",t)
            # if all nodes infected stop
            if t not in timeD:
                inf[t] = len(allInf[:])
            elif len(allInf)==total_number_of_nodes:
                inf[t] = total_number_of_nodes
            else:
                G = create_graph(timeD[t]) # Create the graph that corresponds to this time
                
                graphs = list(nx.connected_component_subgraphs(G)) # find the connected subgraphs
                
                for subgraph in graphs:
                    for sub_node in subgraph.nodes():             # for each node in each connected subgraph
                        if sub_node in allInf:                    # check if node is in the infected list and if it is, add all the nodes
                            
                            #inf[t].extend(subgraph.nodes())   # of the subgraph to the infected list.
                            #inf[t] = [x for x in inf[t] if x not in allInf] # put in the infected dictionary only the new infected nodes
                            
                            allInf.extend(subgraph.nodes())
                            allInf = list(set(allInf))
                            break
                
                inf[t] = len(allInf[:]) # copy variable, not bind
                #seedNodeInf.append(allInf)
        nodesD.append(inf)
    return nodesD
###############################################################################3


Graph = create_aggr_network("/Users/xuecho/Desktop/Data_Highschool.txt")
total_number_of_nodes = Graph.number_of_nodes()
degAss =  nx.degree_assortativity_coefficient(Graph)

G3 = rewiring(Graph,"/Users/xuecho/Desktop/Data_Highschool.txt")
print (nx.degree_assortativity_coefficient(G3))

timeD3 = {}
nodeD3 = {}


#Creation of G2,G3

# G2 is exactly the same as G_data except that the time stamps describing 
#when each temporal link (contact) appears in G_data are randomlized in G2.
time_shuffle = {}
time_shuffle = create_edgelist_reshuffle_time("/Users/xuecho/Desktop/Data_Highschool.txt")
nodes_time_shuffle = stimulation(time_shuffle)
save_obj(nodes_time_shuffle, "nodes_time_shuffle")
plot_infected_per_time(nodes_time_shuffle)


#G3 assign the time stamps in G_data to the linked node pairs (links) randomly.
#A link may receive more than one time stamps, meaning that the two nodes contact more than once. 
#A link G∗3 receives no time stamp means that there is no contact between the corresponding two nodes. 
#G3 is composed of all these contacts.

timeD3 = create_edgelist_from_graph_random_time("/Users/xuecho/Desktop/Data_Highschool.txt",G3, timeD3)

nodeD3 = stimulation(G3,timeD3)

save_obj(nodeD3, "nodeD3")

plot_infected_per_time(nodeD3)
