'''
                        GANGA SINGH MANCHANDA
                        DANIEL GAIVAO LOZANO
                        Imperial College London
                        February 2023
'''
#%%

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch,ArrowStyle
from collections import *
from copy import deepcopy

#%%
    
class DAG:
    
    '''
    To Do:
        - Reduce DAG to the interval
        - Investigate path along greatest and least number of edges
        - Implement measure of jaggedness
    '''
    
    def __init__(self,N):
        # Initialise dictionary structures
        self.N = N
        self.nodes = defaultdict(list)
        self.adj = defaultdict(list)
        self.shortest_dic = defaultdict(list)
        self.shortestNone_dic = defaultdict(list)
        self.longest_dic = defaultdict(list)
        self.longestNone_dic = defaultdict(list)
        self.geodesic_dic = defaultdict(list)
        self.shortpaths = defaultdict(list)
        self.shortpathsNone = defaultdict(list)
        self.shortestnodespdic = [defaultdict(list) for i in range(self.N)]
        self.shortestnodesNonedic = defaultdict(list)
        self.longpaths = defaultdict(list)
        self.longpathsNone = defaultdict(list)
        self.longestnodespdic = [defaultdict(list) for i in range(self.N)]
        self.longestnodesNonedic = defaultdict(list)
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
                    R = 3 / np.sqrt(self.N)
                    if dist < R: 
                        d_dict = defaultdict(list)
                        #d_dict[2] = dist
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
        for p in ps:
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
            #self.longest_dic[p] = dist[self.N-1]
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


    def shortnum(self):
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
                    if dist[node] > dist[i] + 1:
                        dist[node] = dist[i] + 1                      
                        self.shortestnodesNonedic[node] = i            
                self.shortestNone_dic[None] = dist[self.N-1]

            j = self.N - 1
            self.shortpathsNone[None].append(j)
            while True:
                previous_node = self.shortestnodesNonedic[j]
                self.shortpathsNone[None].append(previous_node)
                j = previous_node
                if previous_node == 0:
                    break
            self.shortpathsNone[None] = self.shortpathsNone[None][::-1]

    def longnum(self):
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
                    if (dist[i[0]] < dist[u] + 1):
                        dist[i[0]] = dist[u] + 1
                        self.longestnodesNonedic[i[0]] = u               
                self.longestNone_dic[None] = dist[self.N-1]
        #self.longestNone_dic[None] = dist[self.N-1]

        j = self.N - 1
        self.longpathsNone[None].append(j)
        while True:
            previous_node = self.longestnodesNonedic[j]
            self.longpathsNone[None].append(previous_node)
            j = previous_node
            if previous_node == 0:
                break
        self.longpathsNone[None] = self.longpathsNone[None][::-1]

    def show(self,ps,showedges=False,showdistances=False):
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
                if showedges == True:
                    for i in range(self.N-1):
                        for j in range(len(self.adj[i])):
                            plt.arrow(self.nodes[i][0],self.nodes[i][1],self.nodes[self.adj[i][j][0]][0]-self.nodes[i][0],self.nodes[self.adj[i][j][0]][1]-self.nodes[i][1],length_includes_head = True, color = 'red')

                for i in range(len(self.longpaths[p])-1):
                    plt.arrow(self.nodes[self.longpaths[p][i]][0],self.nodes[self.longpaths[p][i]][1],self.nodes[self.longpaths[p][i+1]][0]-self.nodes[self.longpaths[p][i]][0],self.nodes[self.longpaths[p][i+1]][1]-self.nodes[self.longpaths[p][i]][1],width=0.005,length_includes_head=True,color='darkorange')
                for i in range(len(self.shortpaths[p])-1):
                    plt.arrow(self.nodes[self.shortpaths[p][i]][0],self.nodes[self.shortpaths[p][i]][1],self.nodes[self.shortpaths[p][i+1]][0]-self.nodes[self.shortpaths[p][i]][0],self.nodes[self.shortpaths[p][i+1]][1]-self.nodes[self.shortpaths[p][i]][1],width=0.005,length_includes_head=True,color='dodgerblue')
                if showdistances == True:                   
                    print('Shortest distance for p=' + str(p) + ': ' + str(self.shortest_dic[p]))            
                    print('Longest distance for p=' + str(p) + ': ' + str(self.longest_dic[p]))
                    print('')            

    def shownum(self,ps,showedges=False,showdistances=False):     
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
            if showedges == True:
                for i in range(self.N-1):
                    for j in range(len(self.adj[i])):
                        plt.arrow(self.nodes[i][0],self.nodes[i][1],self.nodes[self.adj[i][j][0]][0]-self.nodes[i][0],self.nodes[self.adj[i][j][0]][1]-self.nodes[i][1],length_includes_head = True, color = 'red')
      
            for i in range(len(self.longpathsNone[None])-1):
                plt.arrow(self.nodes[self.longpathsNone[None][i]][0],self.nodes[self.longpathsNone[None][i]][1],self.nodes[self.longpathsNone[None][i+1]][0]-self.nodes[self.longpathsNone[None][i]][0],self.nodes[self.longpathsNone[None][i+1]][1]-self.nodes[self.longpathsNone[None][i]][1],width=0.005,length_includes_head=True,color='darkorange')
            for i in range(len(self.shortpathsNone[None])-1):
                plt.arrow(self.nodes[self.shortpathsNone[None][i]][0],self.nodes[self.shortpathsNone[None][i]][1],self.nodes[self.shortpathsNone[None][i+1]][0]-self.nodes[self.shortpathsNone[None][i]][0],self.nodes[self.shortpathsNone[None][i+1]][1]-self.nodes[self.shortpathsNone[None][i]][1],width=0.005,length_includes_head=True,color='dodgerblue')
            if showdistances == True:
                print('Shortest distance: ' + str(self.shortestNone_dic[None]) + ' nodes')            
                print('Longest distance: ' + str(self.longestNone_dic[None]) + ' nodes')
                print('')
            
    def l_scaling(self):
        s_keys = list(self.shortest_dic.keys())
        s_vals = [self.shortest_dic[j] for j in s_keys]
        s_geodesic_vals = [self.geodesic_dic[j] for j in s_keys]
        s_norm_vals = [i/j for i,j in zip(s_vals,s_geodesic_vals)]
        l_keys = list(self.longest_dic.keys())
        l_vals = [self.longest_dic[j] for j in l_keys]
        l_geodesic_vals = [self.geodesic_dic[j] for j in l_keys]
        l_norm_vals = [i/j for i,j in zip(l_vals,l_geodesic_vals)]
        with plt.style.context('ggplot'):
            plt.figure(dpi=540)
            plt.grid(False)
            plt.xlabel(r'$p$')
            plt.ylabel(r'$\frac{\ell}{\sqrt{2}}$')
            plt.plot(s_keys,s_norm_vals,'.',color='dodgerblue',label='Shortest')
            plt.plot(l_keys,l_norm_vals,'.',color='darkorange',label='Longest')
            plt.axvline(x=1,color='black',linestyle='--')
            plt.legend()
            plt.show()
            
    def rss_scaling(self,ps):
        s_rss = defaultdict(list)
        l_rss = defaultdict(list)
        for p in ps:
            s_rss[p] = 0
            l_rss[p] = 0
            shortpath = self.shortpaths[p]
            longpath = self.longpaths[p]
            for i in range(len(shortpath)):
                s_rss[p] += ((self.nodes[shortpath[i]][1] - self.nodes[shortpath[i]][0])**2)
            for i in range(len(longpath)):    
                l_rss[p] += ((self.nodes[longpath[i]][1] - self.nodes[longpath[i]][0])**2)  
        s = [s_rss[j] for j in ps]
        l = [l_rss[j] for j in ps]
        with plt.style.context('ggplot'):
            plt.figure(dpi=540)
            plt.grid(False)
            plt.xlabel(r'$p$')
            plt.ylabel(r'$\sum(y_i-x_i)^2$')
            plt.plot(ps,s,'.',color='dodgerblue',label='Shortest')
            plt.plot(ps,l,'.',color='darkorange',label='Longest')
            plt.axvline(x=1,color='black',linestyle='--')
            plt.legend()
            plt.show()
            
    def investigate(self,ps):
        self.minkowski(ps)
        self.short(True,ps)
        self.long(True,ps)
        self.l_scaling()
        self.rss_scaling(ps)
 
#%%

ps = [0.5,2.5]
X = DAG(5000)
X.minkowski(ps)
X.short(True,ps)
X.long(True,ps)
X.shortnum()
X.longnum()
X.l_scaling()
X.rss_scaling(ps)
X.show(ps,showdistances=True)
X.shownum(ps,showdistances=True)