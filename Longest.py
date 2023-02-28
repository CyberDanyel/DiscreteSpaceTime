'''
                        GANGA SINGH MANCHANDA
                        Imperial College London
                        January 2023
'''
#%%

import numpy as np
import matplotlib.pyplot as plt

#%%

def topologicalSortUtil(v):
    global Stack, visited, adj
    
    visited[v] = True
    for i in adj[v]:
        if (not visited[i[0]]):
            topologicalSortUtil(i[0])
    Stack.append(v)

def longestPath(s):
    global Stack, visited, adj, V
    
    dist = [-10**9 for i in range(V)]
    for i in range(V):
        if (visited[i] == False):
            topologicalSortUtil(i)

    dist[s] = 0
    while (len(Stack) > 0):
        u = Stack[-1]
        del Stack[-1]
 
        if (dist[u] != 10**9):
            for i in adj[u]:
                if (dist[i[0]] < dist[u] + i[1]):
                    dist[i[0]] = dist[u] + i[1]
 
    print("INF ",end="") if (dist[V-1] == -10**9) else print(dist[V-1],end=" ")

#%%

n = 5    
distances = np.zeros((n,n))

for i in range(n):
    for j in range(n):
        distances[i][j] = 1 + (i * j)
        distances[i][i] = 0       
distances = np.triu(distances,+1)

adj = [[] for i in range(n+1)]
for i in range(n):
    for j in range(n):
        if distances[i][j] != 0:
            adj[i].append([j,distances[i][j]])

#%%

V, Stack, visited = n, [], [False for i in range(n+1)]
adj = adj
longestPath(0)