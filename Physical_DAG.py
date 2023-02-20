'''
                        GANGA SINGH MANCHANDA
                        Imperial College London
                        January 2023
'''
#%%

from Network_Elements import *
from math import inf
import copy
#%%

class DAG:
    
    '''
    
    This class defines a DAG object.
   
    Attributes
    ----------
    n : integer
        number of nodes
    n_edges : integer
              number of edges
    nodes : list
            list of node objects
    edges : list
            list of edge objects
    Adjacency : numpy array
                matrix of adjacent nodes
    Distance : numpy array
               matrix of adjacent nodes weights
               
    Methods
    -------
    show_shortest : plots and returns shortest path

    '''
    
    def __init__(self,n):
        self.n = n
        self.nodes = []
        self.edges = [[None for x in range(self.n)] for y in range(self.n)]
        self.distances = np.zeros((self.n,self.n))
        
        self.nodes.append(Node(position=np.array([0,0]),vertex=True, identifier = 0))
        for i in range(self.n - 2):
            x = np.random.uniform(0,1)
            y = np.random.uniform(0,1)
            self.nodes.append(Node(position=np.array([x,y]),vertex=False, identifier = i+1))
        self.nodes.append(Node(position=np.array([1,1]),vertex=True, identifier = self.n-1))
    
        for i in range(self.n):
            a = self.nodes[i]
            for j in range(self.n):
                b = self.nodes[j]
                if b.pos()[0] > a.pos()[0] and b.pos()[1] > a.pos()[1]:
                    edge = Edge(a,b)
                    f = edge.dist([2])[0] # Radial Inclusion (using p = 2)
                    #f = np.random.uniform(0,1) # Probabilistic Inclusion
                    P = 0.8
                    if f < P:
                        self.edges[i][j] = edge
                        self.nodes[i].connect(b)
                        self.distances[i][j] = edge.dist([2])[0]
                        a.issource(cond=True)
                        b.issink(cond=True)
                    else:
                        continue    
                else:
                    continue
        
        unfiltered_nodes = np.array(self.nodes)  
        self.nodes = np.array(self.nodes)        
        
        while True:
            lowdegree = [i for i in self.nodes if i.returnsource() == False or i.returnsink() == False]
            if len(lowdegree) == 0:
                self.nodes = [i for i in self.nodes]
                break
            else:
                index = []
                for i in lowdegree:
                    index.append(np.argwhere(i==self.nodes)[0][0])
                for i in index:
                    self.nodes[i].issource(cond=False)
                    self.nodes[i].issink(cond=False)
                self.nodes = [i for i in self.nodes if i not in lowdegree]
                index = []
                for i in lowdegree:
                    index.append(np.argwhere(i==unfiltered_nodes)[0][0])
                for i in index:
                    for j in range(self.n):
                        self.edges[i][j] = None
                        self.edges[j][i] = None
                        self.distances[i][j] = 0
                        self.distances[j][i] = 0
    
        self.edges = [i for x in self.edges for i in x]
        self.edges = [i for i in self.edges if i != None]

        self.n = len(self.nodes)    
        self.n_edges = len(self.edges)
        
        for i in range(len(self.nodes)):
            self.nodes[i].set_identifier(i)
    #def find_distances(p_values):
        
        #for i in range(len(self.edges)):
    def smallest_tent_dist_node(self, nodes):
        dist = []
        for i in range(len(nodes)):
            dist.append(nodes[i].tentative_distance)
        index_min = np.argmin(dist)
        node = nodes[index_min]
        return node
    
    def largest_tent_dist_node(self, nodes):
        dist = []
        for i in range(len(nodes)):
            dist.append(nodes[i].tentative_distance)
        index_max = np.argmax(dist)
        node = nodes[index_max]
        return node
    
    def dijkstra_hand(self, findpath = False, longest = False):
        #path = [self.nodes[0]] #Source node is always in shortest path
        unvisited_set = self.nodes
        self.nodes[0].set_tentative_distance(0)
        current = self.nodes[0]
        start = current.pos()
        neighbours = current.neighbours()
        
        for i in range(len(neighbours)):
            end = neighbours[i].pos()
            if longest == True:
                new_tentative_distance = -(((end[0] - start[0])**2+(end[1] - start[1])**2)**(1/2))
            else:
                new_tentative_distance = (((end[0] - start[0])**2+(end[1] - start[1])**2)**(1/2))                
            if new_tentative_distance < neighbours[i].tentative_distance:  
                neighbours[i].set_tentative_distance(new_tentative_distance)
                neighbours[i].set_previous_vertex(self.nodes[0])
                
        self.nodes[0].mark_visited()
        unvisited_set = [i for i in unvisited_set if i.returnvisited() == False]
        
        while True:
            if longest == True:              
                current = self.largest_tent_dist_node(unvisited_set)
            else:
                current = self.smallest_tent_dist_node(unvisited_set)
            if current.identifier == self.nodes[len(self.nodes)-1].identifier:
                break
            else:
                start = current.pos()
                neighbours = current.neighbours()
                for i in range(len(neighbours)):
                    end = neighbours[i].pos()
                    if longest == True:     
                        new_tentative_distance = current.tentative_distance - ((end[0] - start[0])**2+(end[1] - start[1])**2)**(1/2)
                    else:
                        new_tentative_distance = current.tentative_distance + ((end[0] - start[0])**2+(end[1] - start[1])**2)**(1/2)
                        
                    if new_tentative_distance < neighbours[i].tentative_distance:  
                        neighbours[i].set_tentative_distance(new_tentative_distance)
                        neighbours[i].set_previous_vertex(current)
                current.mark_visited()
                unvisited_set = [i for i in unvisited_set if i.returnvisited() == False]
        
        if findpath == True:
            self.path = [self.nodes[len(self.nodes)-1]]
            self.path_identities = [self.nodes[len(self.nodes)-1].identifier]
            while True:
                current = self.path[len(self.path)-1]
                nextnode = current.previous_vertex
                newidentifier = nextnode.identifier
                self.path.append(nextnode)
                self.path_identities.append(newidentifier)
                if nextnode == self.nodes[0]:
                    self.path = self.path[::-1]
                    self.path_identities = self.path_identities[::-1]
                    break
                else:
                    continue
    def reset_dijkstra(self):
        for i in range(len(self.nodes)):
            self.nodes[i].mark_unvisited()
            self.nodes[i].set_tentative_distance(inf)
    
    def _djikstra(self):
        
        """
        
        https://github.com/crixodia/python-dijkstra/blob/master/dijkstra.py
        
        """
        
        wmat = self.distances
        start = 0
        end = self.n - 1
        n = len(wmat)
        dist = [inf]*n
        dist[start] = wmat[start][start]
        spVertex = [False]*n
        parent = [-1]*n
        path = [{}]*n

        for count in range(n-1):
            minix = inf
            u = 0
            for v in range(len(spVertex)):
                if spVertex[v] == False and dist[v] <= minix:
                    minix = dist[v]
                    u = v
            spVertex[u] = True
            for v in range(n):
                if not(spVertex[v]) and wmat[u][v] != 0 and dist[u] + wmat[u][v] < dist[v]:
                    parent[v] = u
                    dist[v] = dist[u] + wmat[u][v]
        for i in range(n):
            j = i
            s = []
            while parent[j] != -1:
                s.append(j)
                j = parent[j]
            s.append(start)
            path[i] = s[::-1]
            
        self.path = path[end]
        self.path_distance = dist[end]
    
    def show_shortest(self):
        
        self._djikstra()
        
        print ('Shortest Path:',self.path)
        print ('Distance:',self.path_distance)
        
    def show_shortest_new(self, findpath = False):
        self.reset_dijkstra()
        self.dijkstra_hand(findpath)
        
        print ('Shortest Distance: ' + str(self.nodes[len(self.nodes)-1].tentative_distance))
        if findpath == True:
            print ('Shortest Path:' + str(self.path_identities))
            
    def show_longest_new(self, findpath = False):
        self.reset_dijkstra()
        self.dijkstra_hand(findpath, longest = True)
        
        print ('Longest Distance: ' + str(self.nodes[len(self.nodes)-1].tentative_distance))
        if findpath == True:
            print ('Longest Path:' + str(self.path_identities))
            
    def show_longest(self):
        
        self.distances = - self.distances
        self._djikstra()
        
        print ('Longest Path:',self.path)
        print ('Distance:',-self.path_distance)
        
    def draw(self,shortest=False,longest=False):
            
        with plt.style.context('dark_background'):
            plt.figure(dpi=540)
            ax = plt.axes(xlim=(0,1),ylim=(0,1))
            ax.set_aspect('equal',adjustable='box')
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            for i in range(self.n):
                ax.add_patch(self.nodes[i].patch())
            for i in range(self.n_edges):
                ax.add_patch(self.edges[i].edge_patch())
            t = np.linspace(0,1,2)
            plt.plot(t,t,color='white',ls='dashed')
        if shortest == True:
            self.reset_dijkstra()
            self.dijkstra_hand(findpath = True)
            for i in range(len(self.path_identities)-1):
                a = Edge(self.nodes[self.path_identities[i]],self.nodes[self.path_identities[i+1]])
                ax.add_patch(a.shortpath_patch())
        if longest == True:
            self.reset_dijkstra()
            self.dijkstra_hand(findpath = True, longest = True)
            for i in range(len(self.path)-1):
                a = Edge(self.nodes[self.path_identities[i]],self.nodes[self.path_identities[i+1]])
                ax.add_patch(a.longpath_patch())
    
#%%

X = DAG(15)
X.draw(shortest=True, longest = True)
X.show_shortest_new(findpath = True)
X.show_longest_new(findpath = True)


#%%