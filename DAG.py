'''
                        GANGA SINGH MANCHANDA
                        DANIEL GAIVAO LOZANO
                        Imperial College London
                        February 2023
'''
#%%

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, ArrowStyle
from collections import *

#%%
    
class DAG:
    
    '''
    To Do:
        - Reduce the DAG to the interval
        - Carry out investigations for varying p
        - Implement CHI2 measure of path deviation from geodesic
    '''
    def __init__(self,N):
        # Initialise dictionary structures
        self.N = N
        self.nodes = defaultdict(list)
        self.adj = defaultdict(list)
        self.shortest_dic = defaultdict(list)
        self.longest_dic = defaultdict(list)
        self.geodesic_dic = defaultdict(list)
        self.shortestnodespdic = defaultdict(list)
        self.shortpaths = defaultdict(list)
        self.shortestnodespdic = [defaultdict(list) for i in range(self.N)]
        self.longpaths = defaultdict(list)
        self.longestnodespdic = [defaultdict(list) for i in range(self.N)]
        np.random.seed(1)
        
        # Generate random nodes
        self.t,self.s = [],[]
        self.nodes[0] = [0,0]
        for i in range(1,self.N-1):
            x,y = np.random.uniform(0,1,2)
            self.t.append(x)
            self.s.append(y)
            self.nodes[i] = [x,y]
        self.nodes[self.N-1] = [1,1]
        
        # Generate edges using cube rule in adjacency dictionary
        for i in range(self.N):
            for j in range(self.N):
                source = self.nodes[i]
                sink = self.nodes[j]
                dx = sink[0] - source[0]
                dy = sink[1] - source[1]
                if dx > 0 and dy > 0:
                    dist = np.sqrt((dx * dx) + (dy * dy))
                    R = 1
                    if dist < R: 
                        d_dict = defaultdict(list)
                        d_dict[2] = dist
                        self.geodesic_dic[2] = np.sqrt(2)
                        self.adj[i].append([j,d_dict])
                    else:
                        continue
                else:
                    continue    
        
        # # Reduce generated DAG to interval
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
        # self.N = len(self.nodes)
                
    def minkowski(self,ps):
        for i in range(len(ps)):
            p = ps[i]
            self.geodesic_dic[p] = ((self.nodes[self.N-1][0]-self.nodes[0][0])**p + (self.nodes[self.N-1][1]-self.nodes[0][1])**p)**(1/p)
            for i in range(self.N-1):
                for j in range(len(self.adj[i])):
                    dist = ((self.nodes[self.adj[i][j][0]][0]-self.nodes[i][0])**p+(self.nodes[self.adj[i][j][0]][1]-self.nodes[i][1])**p)**(1/p)
                    self.adj[i][j][1][p] = dist
                    
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

    def short(self,findpath=False,ps=None): 
        # Initialize distances to all nodes as infinite and
        # distance to source as 0
        for p in self.adj[0][0][1].keys():
            # Mark all nodes as not visited
            visited = [False]*self.N
            stack =[]
            # Call tSort to store Topological Sort starting from source
            for i in range(self.N):
                if visited[i] == False:
                    self.tSort(0,visited,stack)
            dist = [float("Inf")] * (self.N)
            dist[0] = 0
            # Process nodes in topological order
            while stack:
                # Get the next node from topological order
                i = stack.pop()
                # Update distances of all adjacent nodes
                for node,weights in self.adj[i]:
                    if dist[node] > dist[i] + weights[p]:
                        dist[node] = dist[i] + weights[p]
                        if findpath == True:
                            self.shortestnodespdic[node][p] = i            
                self.shortest_dic[p] = dist[self.N-1]
        if findpath == True:
            for p in ps:
                j = self.N - 1
                self.shortpaths[p].append(j)
                while True:
                    previous_node = self.shortestnodespdic[j][p]
                    self.shortpaths[p].append(previous_node)
                    j = previous_node
                    if previous_node == 0:
                        break
                self.shortpaths[p] = self.shortpaths[p][::-1]
                
    def long(self,findpath=False,ps=None):
        # Mark all nodes as not visited
        for p in self.adj[0][0][1].keys():
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
                        if (dist[i[0]] < dist[u] + i[1][p]):
                            dist[i[0]] = dist[u] + i[1][p]
                            if findpath == True:
                                self.longestnodespdic[i[0]][p] = u               
                    self.longest_dic[p] = dist[self.N-1]
            self.longest_dic[p] = dist[self.N-1]
        if findpath == True:
            for p in ps:
                j = self.N - 1
                self.longpaths[p].append(j)
                while True:
                    previous_node = self.longestnodespdic[j][p]
                    self.longpaths[p].append(previous_node)
                    j = previous_node
                    if previous_node == 0:
                        break
                self.longpaths[p] = self.longpaths[p][::-1]
                    
    def distance_comparison(self):
        keys = list(self.shortest_dic.keys())
        vals = [self.shortest_dic[j] for j in keys]
        geodesic_vals = [self.geodesic_dic[j] for j in keys]
        shortest_norm_vals = [i/j for i,j in zip(vals,geodesic_vals)]
        plt.figure()
        plt.plot(keys, shortest_norm_vals, 'x', color = 'green')
        plt.axvline(x=1, color='r', linestyle='-')
        keys = list(self.longest_dic.keys())
        vals = [self.longest_dic[j] for j in keys]
        longest_norm_vals = [i/j for i,j in zip(vals,geodesic_vals)]
        plt.plot(keys, longest_norm_vals, 'x', color = 'red')
        plt.axvline(x=1, color='r', linestyle='-')
    
    def chi2_comparison_y(self,ps):
        shortchi2 = defaultdict(list)
        longchi2 = defaultdict(list)
        for p in ps:
            shortchi2[p] = 0
            longchi2[p] = 0
            shortpath = self.shortpaths[p]
            longpath = self.longpaths[p]
            for i in range(len(shortpath)):
                shortchi2[p] = shortchi2[p] + ((self.nodes[shortpath[i]][1] - self.nodes[shortpath[i]][0])**2)
            for i in range(len(longpath)):    
                longchi2[p] = longchi2[p] + ((self.nodes[longpath[i]][1] - self.nodes[longpath[i]][0])**2) 
                
        shortchi2s = [shortchi2[j] for j in ps]
        longchi2s = [longchi2[j] for j in ps]
        plt.figure()
        plt.plot(ps,shortchi2s,'x',color='green')
        plt.plot(ps,longchi2s,'x',color='red')
        
    def show(self,ps):
        for p in ps:
            with plt.style.context('ggplot'):
                plt.figure(dpi=540)
                ax = plt.gca()
                ax.set_aspect('equal',adjustable='box')
                plt.grid(False)
                plt.xlabel(r'$x$')
                plt.ylabel(r'$ct$')
                plt.xlim(-0.01,1.01)
                plt.ylim(-0.01,1.01)
                plt.xticks([0,1])
                plt.yticks([0,1])
                plt.plot(self.t,self.s,'.',color='magenta')
                plt.arrow(0,0,1,1,width=0.005,length_includes_head=True,color='black')
                #for i in range(self.N-1):
                    #for j in range(len(self.adj[i])):
                        #plt.arrow(self.nodes[i][0],self.nodes[i][1],self.nodes[self.adj[i][j][0]][0]-self.nodes[i][0],self.nodes[self.adj[i][j][0]][1]-self.nodes[i][1],length_includes_head = True, color = 'black')
                for i in range(len(self.longpaths[p])-1):
                    plt.arrow(self.nodes[self.longpaths[p][i]][0],self.nodes[self.longpaths[p][i]][1],self.nodes[self.longpaths[p][i+1]][0]-self.nodes[self.longpaths[p][i]][0],self.nodes[self.longpaths[p][i+1]][1]-self.nodes[self.longpaths[p][i]][1],width=0.005,length_includes_head=True,color='darkorange')
                for i in range(len(self.shortpaths[p])-1):
                    plt.arrow(self.nodes[self.shortpaths[p][i]][0],self.nodes[self.shortpaths[p][i]][1],self.nodes[self.shortpaths[p][i+1]][0]-self.nodes[self.shortpaths[p][i]][0],self.nodes[self.shortpaths[p][i+1]][1]-self.nodes[self.shortpaths[p][i]][1],width=0.005,length_includes_head=True,color='dodgerblue')

#%%
ps = np.linspace(0.75,1.25,50)
X = DAG(100)
X.minkowski(ps)
X.short(True,ps)
X.long(True,ps)
X.distance_comparison()
X.chi2_comparison_y(ps)
X.show([0.75,1.25])
