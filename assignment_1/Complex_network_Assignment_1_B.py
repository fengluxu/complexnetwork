#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 20:36:17 2018

@author: xuecho
"""
import matplotlib.pylab as pylab
from igraph import *
import networkx as nx
import pandas as pd
#import csv
import pickle
import time
import numpy as np
import matplotlib.pyplot as plt
from random import randint


start = time.time()

"""
FUNCTIONS
"""
def save_obj(obj, name):
    f = open(name + ".dat", "wb")
    pickle.dump(obj, f)
    f.close()

def load_obj(name):
    return pickle.load(open(name + ".dat", "rb"))


"""
def save_obj(obj, name):
    with open('obj/'+ name + '.pkl', 'w+') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)
"""
def create_graph(edgeList):
    G=nx.Graph()    
    for edge in edgeList:
        G.add_edge(edge[0], edge[1])
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

"""
def create_edgelist_per_random_time(filename): #create edgelist 
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

    #create graph from csv
    data = pd.read_table(filename,header=0,delim_whitespace=True)
    extCol = np.array(data.iloc[:,[0,1]]).tolist()
    for row in extCol:
        vertex1 = int(row[0])
        vertex2 = int(row[1])
        G.add_edge(vertex1,vertex2)
    return G   

def intersect(a, b):
    return list(set(a) & set(b))

def temporal_degree(time_dict):
    for t in range(1,max(timeD)+1):
        Gt = create_graph(timeD[t])
        tempDegreeList = list(Gt.degree(Gt.nodes()).values())
        print (tempDegreeList)

def compute_closeness_metric(clossD):    

    for i in clossD:
        clossTupleList.append((i, clossAgg[i]))    

    Cl = [int(tup[0]) for tup in sorted(clossTupleList, key=lambda x: x[1], reverse=True)]    

    for f in np.arange(0.05,0.55,0.05):
        #print f
        lastElement = int(f*len(Cl))
        rRXf[f] = float(len(intersect(R[:lastElement], Cl[:lastElement])))/float(len(R[:lastElement]))
    return rRXf
#######################################################################################################################################################################
################################################################################################################################################################

#Global Variavbles
timeD = {}  # Dictionary with the edgelist per time  {"0":[{1,2},(5,6)...], "1":[(4,7)...]...}
nodesD = [] # Dictionary which contains an inf{} per node

create_edgelist_per_time("/Users/xuecho/Desktop/Data_Highschool.txt", timeD)      # this could have been done only once. However, it is faster o build the edgelist each time 
#save_obj(timeD, "edgelist_dictionary")   #than loading it.
#save_obj(inf, "infected_dictionary")


########################################
#Creation of aggregate graph 
########################################
aggr_graph = create_aggr_network("/Users/xuecho/Desktop/Data_Highschool.txt")
total_number_of_nodes = aggr_graph.number_of_nodes()
#nx.draw(aggr_graph)
#########################################
## Stimulation of the temporal network
#########################################
for node in range(1, total_number_of_nodes+1):
    #for node in range
    print("node", node)
    allInf = [node] # List with all the currently infected nodes
    inf = {}  # Dictionary with the infected nodes per time
  
    
#    for t in timeD: #Creating the network topology per time
    for t in range(1,max(timeD)+1): #Creating the network topology per time
        #print("timesme",t)
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
    #print nodesD    
end = time.time()
print(end - start) #8107.738505363464

###############
# Computations#
###############

###################################################################
# (9) Computation of the Average Number of Infected nodes E|I(t)| #
###################################################################
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
#plt.legend(['Standard Deviation'],
#           loc='upper left',
#           numpoints=1,
#          fancybox=True)
plt.xlabel('time t')
plt.ylabel('E[I(t)]')
plt.title('the average number of infected nodes E[I(t)] together with its error bar')
plt.show()


#############################################
# (10) Ranking of the most influential nodes#
#############################################

infl =[] #[(node, time),..]
rankVar = (80*total_number_of_nodes)/100 #rank variable is 261.6
for i in range(0,len(nodesD)):
    infl.append((i, min(k for k in nodesD[i] if nodesD[i][k] >= rankVar))) 
#The shorter the time is, the more influential the seed node is.
R =  [int(tup[0]) for tup in sorted(infl, key=lambda x: x[1])] # sort in the ascending order 

save_obj(infl, "RTupleG2")
save_obj(R, "R10G2")

########################################
# (11) Degree & Clustering Coefficient #
########################################
degreeTupleList = []
clustTupleList = []

#TODO: not so many transformations! More efficiently! 
degreeList = list(aggr_graph.degree(list(aggr_graph.nodes()))) #[(node, degree),..]
clustD = nx.clustering(aggr_graph)

newnlist = sorted(degreeList, key = lambda x: x[1],reverse=True) 

# transform degree list into a sorted ordered list of nodes according to degree
D =[]
for row in newnlist:
    D.append(int(row[0]))

save_obj(newnlist, "Dtuple10")

print(D)#print degreeList
save_obj(D, "D10")

# transform clust. coef. list into a sorted ordered list of nodes according to clust. coef.
for i in clustD:
    clustTupleList.append((i, clustD[i]))
C = [int(tup[0]) for tup in sorted(clustTupleList, key=lambda x: x[1], reverse=True)]
#print clustList
save_obj(C, "C10")

#######################
# Computation of rRDf #
#######################
rRDf = {}
for f in np.arange(0.05,0.55,0.05):
    print (f)
    lastElement = int(f*len(D))
    rRDf[f] = float(len(intersect(R[:lastElement], D[:lastElement])))/float(len(R[:lastElement]))
    #print f, rRDf[f]
print(rRDf)

######################
# Computation of rRCf#
######################
rRCf = {}
for f in np.arange(0.05,0.55,0.05):
    print (f)
    lastElement = int(f*len(C))
    rRCf[f] = float(len(intersect(R[:lastElement], C[:lastElement])))/float(len(R[:lastElement]))

save_obj(rRDf, "rRDf11")
save_obj(rRCf, "rRCf11")


plot_Dict_2(rRDf, "lightblue")
plot_Dict_2(rRCf, "mediumturquoise")
plt.legend(['rRD(f)','rRC(f)'],
           loc='upper left',
           numpoints=1,
          fancybox=True)
plt.title('Comparison between rRC(f) and rRD(f)')
plt.show()


################################
# (12) Other centrality metrics#
################################

#########################
# Closeness of Temporal #
#########################
rRCLf = {}
clossAggTupleList = []

clossAgg = nx.closeness_centrality(aggr_graph) # or call function
for i in clossAgg:
    clossAggTupleList.append((i, clossAgg[i]))

Cl = [int(tup[0]) for tup in sorted(clossAggTupleList, key=lambda x: x[1], reverse=False)]

for f in np.arange(0.05,0.55,0.05):
    print (f)
    lastElement = int(f*len(Cl))
    rRCLf[f] = float(len(intersect(R[:lastElement], Cl[:lastElement])))/float(len(R[:lastElement]))


save_obj(rRCLf, "rRCLf11")
plot_Dict_2(rRDf, "lightblue")
plot_Dict_2(rRCf, "mediumturquoise")
plot_Dict_2(rRCLf, "mediumpurple")
plt.title('Closeness rRCL(f)')

#plt.show()

#########################
#betweenness of Temporal#
#########################

rRBef = {}
betweenAggTupleList = []

betweenAgg = nx.betweenness_centrality(aggr_graph) # or call function
for i in betweenAgg:
    betweenAggTupleList.append((i, betweenAgg[i]))

Be = [int(tup[0]) for tup in sorted(betweenAggTupleList, key=lambda x: x[1], reverse=False)]

for f in np.arange(0.05,0.55,0.05):
    print (f)
    lastElement = int(f*len(Cl))
    rRBef[f] = float(len(intersect(R[:lastElement], Be[:lastElement])))/float(len(R[:lastElement]))

rRClTf = load_obj("rRClTf")
save_obj(rRBef, "rRBef11")
plot_Dict_2(rRDf, "lightblue")
plot_Dict_2(rRCf, "mediumturquoise")
plot_Dict_2(rRCLf, "mediumpurple")
plot_Dict_2(rRClTf, "pink")


plt.legend(['rRD(f)','rRC(f)','rRCL(f)','rRClTf'],
           loc='upper left',
           numpoints=1,
          fancybox=True)
plt.title('Comparison among rRD(f),rRC(f), rRCL(f) and rRClT(f)')


plt.legend(['rRD(f)','rRC(f)','rRCL(f)','rRBe(f)'],
           loc='upper left',
           numpoints=1,
          fancybox=True)
plt.title('Comparison among rRD(f),rRC(f), rRCL(f) and rRBe(f)')

#plt.show()
