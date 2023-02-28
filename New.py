'''
                        GANGA SINGH MANCHANDA
                        Imperial College London
                        February 2023
'''
#%%

import numpy as np
from collections import *
import random

#%%

class DAG:
    
    def __init__(self,N):
        # Initialise dictionary structures
        self.N = N
        self.nodes = defaultdict(list)
        self.adj = defaultdict(list)
        
        # Generate random nodes
        self.nodes[0].append([0,0])
        for i in range(1,self.N-1):
            x = np.random.uniform(0,1)
            y = np.random.uniform(0,1)
            self.nodes[i].append([x,y])
        self.nodes[self.N-1].append([1,1])
        
        # Generate edges using cube rule in adjacency dictionary
        for i in range(self.N):
            for j in range(self.N):
                source = self.nodes[i][0]
                sink = self.nodes[j][0]
                dx = sink[0] - source[0]
                dy = sink[1] - source[1]
                if dx > 0 and dy > 0:
                    dist = np.sqrt((dx * dx) + (dy * dy))
                    R = 1
                    if dist < R: 
                        self.adj[i].append([j,dist])
                    else:
                        continue
                else:
                    continue    
                
        # Reduce generated DAG to interval
        # while True:
        #     # Delete nodes with no outgoing edges
        #     count = 0
        #     for i in range(len(self.nodes)):
        #         if self.adj[i] == []:
        #             del self.nodes[i]
        #             del self.adj[i]
        #             count += 1
        #     # Delete nodes with no incoming edges
        #     edges = list(self.nodes.values())
        #     l = []
        #     for i in edges:
        #         if i != []:
        #             l.append(i[0][0])
        #     sinks = list(Counter(l).keys())
        #     allnodes = list(self.nodes.keys())
        #     missing = [x for x in allnodes if x not in sinks]
        #     for i in missing:
        #         del self.nodes[i]
        #         del self.adj[i]
        #     #Break loop
        #     if count == 0 and len(missing) == 0:
        #         break
        #     else:
        #         continue
            
    def tSort(self,v,visited,stack):        
        # Mark current node as visited
        visited[v] = True
 
        # Recur for all nodes adjacent to this node
        if v in self.adj.keys():
            for node,weight in self.adj[v]:
                if visited[node] == False:
                    self.tSort(node,visited,stack)
 
        # Push current node to stack which stores topological sort
        stack.append(v)

    def short(self):
        # Mark all nodes as not visited
        visited = [False]*self.N
        stack =[]
 
        # Call tSort to store Topological Sort starting from source
        for i in range(self.N):
            if visited[i] == False:
                self.tSort(0,visited,stack)
 
        # Initialize distances to all nodes as infinite and
        # distance to source as 0
        dist = [float("Inf")] * (self.N)
        dist[0] = 0
 
        # Process nodes in topological order
        while stack:
 
            # Get the next node from topological order
            i = stack.pop()
 
            # Update distances of all adjacent nodes
            for node,weight in self.adj[i]:
                if dist[node] > dist[i] + weight:
                    dist[node] = dist[i] + weight
         
        print (dist[self.N-1])
           
    def long(self):
        # Mark all nodes as not visited
        visited = [False]*self.N
        stack =[]
        dist = [-10**9 for i in range(self.N)]
 
        # Call tSort to store Topological Sort starting from all nodes
        for i in range(self.N):
            if (visited[i] == False):
                self.tSort(i,visited,stack)
 
        # Initialize distances to all nodes as infinite and
        # distance to source as 0
        dist[0] = 0
 
        # Process nodes in topological order
        while (len(stack) > 0):
       
            # Get the next node from topological order
            u = stack[-1]
            del stack[-1]
 
            # Update distances of all adjacent nodes
            if (dist[u] != 10**9):
                for i in self.adj[u]:
                    if (dist[i[0]] < dist[u] + i[1]):
                        dist[i[0]] = dist[u] + i[1]

        print (dist[self.N-1])
 

#%%

X = DAG(1000)
X.short()
X.long()
