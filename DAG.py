'''
                        GANGA SINGH MANCHANDA
                        DANIEL GAIVAO LOZANO
                        Imperial College London
                        February 2023
'''
#%%

import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import time as time
from scipy.optimize import curve_fit

def linear(x,m,c):
    return (m * x) + c

plt.style.use('gangaplot')

#%%
    
class DAG:
    
    def __init__(self,N):
        # Initialise dictionary structures
        self.N = N
        self.nodes = defaultdict(list)
        self.adj = defaultdict(list)
        self.sources = defaultdict(list)
        self.sinks = defaultdict(list)
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
        self.greedy_short_path_length = defaultdict(list)
        self.greedy_short_path_dic = defaultdict(list)
        self.greedy_long_path_length = defaultdict(list)
        self.greedy_long_path_dic = defaultdict(list)

        np.random.seed(3)
        
        # Generate random nodes
        self.nodes[0] = [0,0]
        self.sources[0] = True
        self.sinks[0] = True
        for i in range(1,self.N-1):
            x,y = np.random.uniform(0,1,2)
            self.nodes[i] = [x,y]
            self.sources[i] = False
            self.sinks[i] = False
        self.nodes[self.N-1] = [1,1]
        self.sources[self.N-1] = True
        self.sinks[self.N-1] = True
        
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
                    #R = 0.5
                    if dist < R: 
                        self.sources[i] = True
                        self.sinks[j] = True
                        d_dict = defaultdict(list)
                        #d_dict[2] = dist
                        self.geodesic_dic[2] = np.sqrt(2)
                        self.adj[i].append([j,d_dict])
                    else:
                        continue
                else:
                    continue

        # deletion_number = 1
        # i = 0
        # while deletion_number != 0:
        #     delete_nodes = []
        #     deletion_number = 0
        #     for node in self.adj.keys():
        #         delete_list = []
        #         for i in range(len(self.adj[node])):
        #             current_node = self.adj[node][i][0]
        #             if self.sources[current_node] != True:
        #                 delete_list.append(i)
        #                 deletion_number = deletion_number + 1
        #         delete_list = delete_list[::-1]
        #         for delete in delete_list:
        #             del self.adj[node][delete]    
                    
        #         if self.adj[node] == []:
        #             delete_nodes.append(node)
        #             self.sources[node] = False
        #     for delete in delete_nodes:
        #         del self.adj[delete]
        #         del self.nodes[delete]
        #     i = i+1
        # while True:
        #     delete_list = []
        #     delete_num = 0
        #     for node in self.nodes:
        #         if self.sinks[node] != True:   
        #             delete_list.append(node)
        #             delete_num = delete_num + 1
        #     if delete_num == 0:
        #         break
        #     for delete in delete_list:
        #         del self.nodes[delete]
        #         del self.adj[delete]
                
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

        j = self.N - 1
        self.longpathsNone[None].append(j)
        while True:
            previous_node = self.longestnodesNonedic[j]
            self.longpathsNone[None].append(previous_node)
            j = previous_node
            if previous_node == 0:
                break
        self.longpathsNone[None] = self.longpathsNone[None][::-1]
        

    def greedy_short(self,ps):
        greedy_adj = {key: value for key, value in self.adj.items() if value != []}
        for p in ps:
            current_node = 0
            greedy_length = 0
            while current_node != self.N-1:
                edge = np.inf
                for i in range(len(greedy_adj[current_node])):
                    new_edge = greedy_adj[current_node][i][1][p]
                    if new_edge < edge:
                        edge = new_edge
                        next_node = greedy_adj[current_node][i][0]
                self.greedy_short_path_dic[p].append(next_node)
                greedy_length = greedy_length + edge
                current_node = next_node
            self.greedy_short_path_length[p] = greedy_length
            
    def greedy_long(self,ps):
        greedy_adj = {key: value for key, value in self.adj.items() if value != []}
        for p in ps:
            current_node = 0
            greedy_length = 0
            while current_node != self.N-1:
                edge = 0
                for i in range(len(greedy_adj[current_node])):
                    new_edge = greedy_adj[current_node][i][1][p]
                    if new_edge > edge:
                        edge = new_edge
                        next_node = greedy_adj[current_node][i][0]
                self.greedy_long_path_dic[p].append(next_node)
                greedy_length = greedy_length + edge
                current_node = next_node
            self.greedy_long_path_length[p] = greedy_length
            
    def show(self,ps,showedges=False,showdistances=False):
        for p in ps:
            ax = plt.gca()
            ax.set_aspect('equal',adjustable='box')
            plt.xlabel(r'$x$',fontsize=18)
            plt.ylabel(r'$ct$',rotation=0,ha='right',fontsize=18)
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.xticks([0,1])
            plt.yticks([0,1])
            for node in self.nodes:
                if self.sinks[node] and self.sources[node] == True:
                    plt.plot(self.nodes[node][0],self.nodes[node][1],'.',color='magenta',alpha=0.1)
            plt.arrow(0,0,1,1,width=0.005,length_includes_head=True,color='black',label='Geodesic')
            plt.arrow(0,0,0,0,color='dodgerblue',label='Shortest')
            plt.arrow(0,0,0,0,color='darkorange',label='Longest')
            if showedges == True:
                for i in range(self.N-1):
                    for j in range(len(self.adj[i])):
                        plt.arrow(self.nodes[i][0],self.nodes[i][1],self.nodes[self.adj[i][j][0]][0]-self.nodes[i][0],self.nodes[self.adj[i][j][0]][1]-self.nodes[i][1],length_includes_head = True, color = 'red')
            for i in range(len(self.longpaths[p])-1):
                plt.arrow(self.nodes[self.longpaths[p][i]][0],self.nodes[self.longpaths[p][i]][1],self.nodes[self.longpaths[p][i+1]][0]-self.nodes[self.longpaths[p][i]][0],self.nodes[self.longpaths[p][i+1]][1]-self.nodes[self.longpaths[p][i]][1],width=0.005,length_includes_head=True,color='darkorange')
            for i in range(len(self.shortpaths[p])-1):
                plt.arrow(self.nodes[self.shortpaths[p][i]][0],self.nodes[self.shortpaths[p][i]][1],self.nodes[self.shortpaths[p][i+1]][0]-self.nodes[self.shortpaths[p][i]][0],self.nodes[self.shortpaths[p][i+1]][1]-self.nodes[self.shortpaths[p][i]][1],width=0.005,length_includes_head=True,color='dodgerblue')
            plt.legend(bbox_to_anchor=(1,0),loc="lower left",borderaxespad=0)
            plt.show()     
            if showdistances == True:                   
                print('Shortest distance for p=' + str(p) + ': ' + str(self.shortest_dic[p]))            
                print('Longest distance for p=' + str(p) + ': ' + str(self.longest_dic[p]))
                print('') 

    def shownum(self,ps,showedges=False,showdistances=False):     
        ax = plt.gca()
        ax.set_aspect('equal',adjustable='box')
        plt.xlabel(r'$x$',fontsize=18)
        plt.ylabel(r'$ct$',rotation=0,ha='right',fontsize=18)
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.xticks([0,1])
        plt.yticks([0,1])
        for node in self.nodes:
            if self.sinks[node] and self.sources[node] == True:
                plt.plot(self.nodes[node][0],self.nodes[node][1],'.',color='magenta',alpha=0.1)
        plt.arrow(0,0,1,1,width=0.005,length_includes_head=True,color='black',label='Geodesic')
        plt.arrow(0,0,0,0,color='dodgerblue',label='Shortest')
        plt.arrow(0,0,0,0,color='darkorange',label='Longest')
        if showedges == True:
            for i in range(self.N-1):
                for j in range(len(self.adj[i])):
                    plt.arrow(self.nodes[i][0],self.nodes[i][1],self.nodes[self.adj[i][j][0]][0]-self.nodes[i][0],self.nodes[self.adj[i][j][0]][1]-self.nodes[i][1],length_includes_head = True, color = 'red')
        for i in range(len(self.longpathsNone[None])-1):
            plt.arrow(self.nodes[self.longpathsNone[None][i]][0],self.nodes[self.longpathsNone[None][i]][1],self.nodes[self.longpathsNone[None][i+1]][0]-self.nodes[self.longpathsNone[None][i]][0],self.nodes[self.longpathsNone[None][i+1]][1]-self.nodes[self.longpathsNone[None][i]][1],width=0.005,length_includes_head=True,color='darkorange')
        for i in range(len(self.shortpathsNone[None])-1):
            plt.arrow(self.nodes[self.shortpathsNone[None][i]][0],self.nodes[self.shortpathsNone[None][i]][1],self.nodes[self.shortpathsNone[None][i+1]][0]-self.nodes[self.shortpathsNone[None][i]][0],self.nodes[self.shortpathsNone[None][i+1]][1]-self.nodes[self.shortpathsNone[None][i]][1],width=0.005,length_includes_head=True,color='dodgerblue')
        plt.legend(bbox_to_anchor=(1,0),loc="lower left",borderaxespad=0)
        plt.show()
        if showdistances == True:
            print('Shortest distance: ' + str(self.shortestNone_dic[None]) + ' nodes')            
            print('Longest distance: ' + str(self.longestNone_dic[None]) + ' nodes')
            print('')
            
    def showgreedy(self,ps,showedges=False,showdistances=False):
        for p in ps:
            ax = plt.gca()
            ax.set_aspect('equal',adjustable='box')
            plt.xlabel(r'$x$',fontsize=18)
            plt.ylabel(r'$ct$',rotation=0,ha='right',fontsize=18)
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.xticks([0,1])
            plt.yticks([0,1])
            for node in self.nodes:
                if self.sinks[node] and self.sources[node] == True:
                    plt.plot(self.nodes[node][0],self.nodes[node][1],'.',color='magenta',alpha=0.1)
            plt.arrow(0,0,1,1,width=0.005,length_includes_head=True,color='black',label='Geodesic')
            plt.arrow(0,0,0,0,color='dodgerblue',label='Shortest')
            plt.arrow(0,0,0,0,color='darkorange',label='Longest')
            if showedges == True:
                for i in range(self.N-1):
                    for j in range(len(self.adj[i])):
                        plt.arrow(self.nodes[i][0],self.nodes[i][1],self.nodes[self.adj[i][j][0]][0]-self.nodes[i][0],self.nodes[self.adj[i][j][0]][1]-self.nodes[i][1],length_includes_head = True, color = 'red')
            for i in range(len(self.greedy_short_path_dic[p])-1):
                plt.arrow(self.nodes[self.greedy_short_path_dic[p][i]][0],self.nodes[self.greedy_short_path_dic[p][i]][1],self.nodes[self.greedy_short_path_dic[p][i+1]][0]-self.nodes[self.greedy_short_path_dic[p][i]][0],self.nodes[self.greedy_short_path_dic[p][i+1]][1]-self.nodes[self.greedy_short_path_dic[p][i]][1],width=0.005,length_includes_head=True,color='darkorange')
            for i in range(len(self.greedy_long_path_dic[p])-1):
                plt.arrow(self.nodes[self.greedy_long_path_dic[p][i]][0],self.nodes[self.greedy_long_path_dic[p][i]][1],self.nodes[self.greedy_long_path_dic[p][i+1]][0]-self.nodes[self.greedy_long_path_dic[p][i]][0],self.nodes[self.greedy_long_path_dic[p][i+1]][1]-self.nodes[self.greedy_long_path_dic[p][i]][1],width=0.005,length_includes_head=True,color='dodgerblue')
            plt.legend(bbox_to_anchor=(1,0),loc="lower left",borderaxespad=0)
            plt.show()
            if showdistances == True:                   
                print('Shortest distance for p=' + str(p) + ': ' + str(self.shortest_dic[p]))            
                print('Longest distance for p=' + str(p) + ': ' + str(self.longest_dic[p]))
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
        plt.xlabel(r'$p$',fontsize=18)
        plt.ylabel(r'$\frac{\ell}{\ell_g}$',rotation=0,ha='right',fontsize=18)
        plt.plot(s_keys,s_norm_vals,'.',color='dodgerblue',label='Shortest')
        plt.plot(l_keys,l_norm_vals,'.',color='darkorange',label='Longest')
        plt.axvline(x=1,color='black',linestyle='--')
        plt.legend(bbox_to_anchor=(1,0),loc="lower left",borderaxespad=0,fontsize=22)
        plt.show()
        
        s_keys = list(self.greedy_short_path_length.keys())
        s_vals = [self.greedy_short_path_length[j] for j in s_keys]
        s_geodesic_vals = [self.geodesic_dic[j] for j in s_keys]
        s_norm_vals = [i/j for i,j in zip(s_vals,s_geodesic_vals)]
        l_keys = list(self.greedy_long_path_length.keys())
        l_vals = [self.greedy_long_path_length[j] for j in l_keys]
        l_geodesic_vals = [self.geodesic_dic[j] for j in l_keys]
        l_norm_vals = [i/j for i,j in zip(l_vals,l_geodesic_vals)]
        plt.xlabel(r'$p$',fontsize=18)
        plt.ylabel(r'$\frac{\ell}{\ell_g}$',rotation=0,ha='right',fontsize=18)
        plt.plot(s_keys,s_norm_vals,'.',color='dodgerblue',label='Shortest')
        plt.plot(l_keys,l_norm_vals,'.',color='darkorange',label='Longest')
        plt.axvline(x=1,color='black',linestyle='--')
        plt.legend(bbox_to_anchor=(1,0),loc="lower left",borderaxespad=0,fontsize=22)
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
        plt.xlabel(r'$p$',fontsize=18)
        plt.ylabel(r'RSS',rotation=0,ha='right',fontsize=18)
        plt.plot(ps,s,'.',color='dodgerblue',label='Shortest')
        plt.plot(ps,l,'.',color='darkorange',label='Longest')
        plt.axvline(x=1,color='black',linestyle='--')
        plt.legend(bbox_to_anchor=(1,0),loc="lower left",borderaxespad=0,fontsize=22)
        plt.show()
        
        s_rss = defaultdict(list)
        l_rss = defaultdict(list)
        for p in ps:
            s_rss[p] = 0
            l_rss[p] = 0
            shortpath = self.greedy_short_path_dic[p]
            longpath = self.greedy_long_path_dic[p]
            for i in range(len(shortpath)):
                s_rss[p] += ((self.nodes[shortpath[i]][1] - self.nodes[shortpath[i]][0])**2)
            for i in range(len(longpath)):    
                l_rss[p] += ((self.nodes[longpath[i]][1] - self.nodes[longpath[i]][0])**2)  
        s = [s_rss[j] for j in ps]
        l = [l_rss[j] for j in ps]
        plt.xlabel(r'$p$',fontsize=18)
        plt.ylabel(r'RSS',rotation=0,ha='right',fontsize=18)
        plt.plot(ps,s,'.',color='dodgerblue',label='Shortest')
        plt.plot(ps,l,'.',color='darkorange',label='Longest')
        plt.axvline(x=1,color='black',linestyle='--')
        plt.legend(bbox_to_anchor=(1,0),loc="lower left",borderaxespad=0,fontsize=22)
        plt.show()
            
    def investigate(self,ps,showdistances=False):
        pshow = [ps[0],ps[-1]]
        self.minkowski(ps)
        self.short(True,ps)
        self.long(True,ps)
        self.greedy_short(ps)
        self.greedy_long(ps)
        self.shortnum()
        self.longnum()
        self.l_scaling()
        self.rss_scaling(ps)# print rss values for num too
        self.show(pshow,showdistances)
        self.shownum(pshow,showdistances)
        self.showgreedy(pshow,showdistances)
        
    def complexity(self,ps):        
        s = time.time()
        self.minkowski(ps)
        self.long(True,ps)
        e = time.time()
        t = e - s
        return t
                
#%%

ps = np.linspace(-0.5,2.5,100)
X = DAG(10000)
X.investigation(ps)

#%%

sizes = np.linspace(1000,100000,10,dtype=int)
times = []
num = 10

for i in sizes:
    t = []
    for j in range(num):
        t.append(DAG(i).complexity([1]))
    times.append(np.mean(t))
    
plt.xlabel(r'$N$',fontsize=18)
plt.ylabel(r'$t\,[\mathrm{s}]$',rotation=0,ha='right',fontsize=18)
plt.scatter(sizes,times,label='Data')
n = np.linspace(sizes[0],sizes[-1],100)
fit,cov = curve_fit(linear,sizes,times)
plt.plot(n,linear(n,*fit),color='dodgerblue',label='Curve Fit')
plt.legend(bbox_to_anchor=(1,0),loc="lower left",borderaxespad=0,fontsize=18)
plt.show()
print (fit)
