#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 21:46:57 2018

@author: Fenglu Xu
"""

"""
VARIABLES 
"""
import numpy as np
import igraph as ig
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

#print(ig.__version__)

data = pd.read_table('/Users/xuecho/Desktop/Data_Highschool.txt',header=0
                     ,delim_whitespace=True)
extCol = np.array(data.iloc[:,[0,1]]).tolist() #Extract len(extCol) = 188508 nodes columns 

remDup = list(); #list to store nodes

vertices = [i for i in range(3)]

"""
SEPARATE DUPLICATES
"""

for i in extCol:
    if i not in remDup:
        remDup.append(i)   # 5818 after remove duplicates
#list to store sorted unique nodes connections
newlist = sorted(remDup, key = lambda x: x[0]) 

remDupLink = list();
for i in extCol:
    inverlinks = [0,0]
    links = i;
    inverlinks[0] = inverlinks[1]
    inverlinks[1] = inverlinks[0]
    if i not in remDupLink and inverlinks not in remDupLink :
        remDupLink.append(i)

newList = sorted(remDupLink, key = lambda x: x[0])

#Create graph based on node connections
g = ig.Graph(edges=remDup, directed=False)
layout = g.layout("kk")

#QUESTION 1
#What is the number of nodes N, the number of links L, the link density p, the 
#average degree E[D] and the degree variance Var[D]?

#Number of Nodes
nodes = g.vcount()
#The total amount of nodes is 328. 

#Number of links 
links = g.ecount()
g#The total amount of links is 5818.

#Link density
grDen = g.density(g)
#According with the console grDen =  0.1078286010823634

#the average degree E[D]
avgDegVec = np.average(np.asarray(g.degree())) # asarray doesn't create a new copy
#According with the console the average degree vector is 35.475609756097562

#the degree variance Var[D]
p = g.degree()
p = p[1:-1]
varDeg = np.var(np.asarray(p))
#According with the console the variance is 182.68060521660581

#QUESTION 2
distr = g.degree_distribution()
print(distr) #N = 328, mean +- sd: 35.4756 +- 13.6397

xs, ys = zip(*[(left, count) for left, _, count in 
g.degree_distribution().bins()])
pylab.bar(xs, ys)
pylab.xlabel('degrees')
pylab.ylabel('#nodes')
pylab.title('degree distribution')
pylab.show()

#https://jeremykun.com/2013/08/22/the-erdos-renyi-random-graph/
#https://en.wikipedia.org/wiki/Barab%C3%A1si%E2%80%93Albert_model
# ER random graph
#why？
visual_style = {}
visual_style["layout"] = layout
visual_style["vertex_size"] = 20
visual_style["vertex_color"] = "lightblue"
visual_style["vertex_label"] = g.vs["name"]
visual_style["bbox"] = (1000, 1000)
visual_style["margin"] = 10
ig.plot(g, **visual_style)

#edges = np.array(data.iloc[:,[0,1]]).tolist()
#gx = ig.Graph(vertex_attrs={"label":vertices}, edges=edges, directed=False)
#gx.degree_distribution()
#layout = g.layout("kk")


#Question 3

#What is the degree correlation (assortativity)?  What is its physical meaning?
degCorrel = g.assortativity_degree(directed=False)
#The degree correlation is equal to 0.03317576993604877


#QUESTION 4

#What is the clustering coefficient C?
clustCoef = g.transitivity_undirected(mode="nan")
#According with the console result the clustCoef = 0.44442135612150807

#QUESTION 5

#What is the average hopcount E[H] of the shortest paths between all node pairs?  
#What is the diameter H_max?
avgHop = g.average_path_length()
#The average hopcount E[H] of the shortest paths between all node pairs is
#2.1594341569576554

diamG = g.diameter(directed=False,unconn=True,weights=None)
# diameter = 4


betnGr = g.closeness() #array with the closeness property for each link
avgBetGr = np.average(betnGr)#average result of the closeness array
maxBet = np.nanmax(betnGr) #max value of the closeness array.
#Below I compute the betweenness values as well
betwe = g.betweenness()
avgB = np.average(betwe)
maxBe = np.nanmax(avgB)


#QUESTION 6
#Has this network the small-world property? Justify your conclusion quantitatively.

#QUESTION 7
#What is the largest eigenvalue (spectral radius) lambda 1 of the adjacency matrix?
#create the adjacency matrix 
adjMat = g.get_adjacency()  
#another matrix to read the value
adjMatRead = np.array(adjMat.data)
eigVec = np.linalg.eigvals(adjMatRead)
maxEig = np.amax(eigVec)
#the max eigenvalue from the matrix is 41.231605032525721 

#QUESTION 8

#What is the second smallest eigenvalue of the Laplacian matrix (algebraic connectivity)?
lapVal = g.laplacian(normalized=True)
eigVector = np.linalg.eigvals(lapVal)
eigSort = np.sort(eigVector)
#create the laplacian matrix
#extract the eigenvalue vector and sort it to identify the second smallest value
#according with my console the second smallest eigenvalue is 1.38777878e-16. 


