'''
                        GANGA SINGH MANCHANDA
                        Imperial College London
                        January 2023
'''
#%%

from Network_Elements import *
from math import inf

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
        
        self.nodes.append(Node(position=np.array([0,0]),vertex=True))
        for i in range(self.n - 2):
            x = np.random.uniform(0,1)
            y = np.random.uniform(0,1)
            self.nodes.append(Node(position=np.array([x,y]),vertex=False))
        self.nodes.append(Node(position=np.array([1,1]),vertex=True))
    
        for i in range(self.n):
            a = self.nodes[i]
            for j in range(self.n):
                b = self.nodes[j]
                edge = Edge(a,b)
                f = edge.dist() # Radial Inclusion
                #f = np.random.uniform(0,1) # Probabilistic Inclusion
                P = 0.6
                if b.pos()[0] > a.pos()[0] and b.pos()[1] > a.pos()[1] and f < P:
                    self.edges[i][j] = edge
                    self.distances[i][j] = edge.dist()
                    a.issource(cond=True)
                    b.issink(cond=True)
                else:
                    continue      
        
        unfiltered_nodes = np.array(self.nodes)  
        self.nodes = np.array(self.nodes)        
        u = True
        while u == True:
            lowdegree = [i for i in self.nodes if i.returnsource() == False or i.returnsink() == False]
            if len(lowdegree) == 0:
                u = False
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
            self._djikstra()
            for i in range(len(self.path)-1):
                a = Edge(self.nodes[self.path[i]],self.nodes[self.path[i+1]])
                ax.add_patch(a.shortpath_patch())
        if longest == True:
            self.distances = - self.distances
            self._djikstra()
            for i in range(len(self.path)-1):
                a = Edge(self.nodes[self.path[i]],self.nodes[self.path[i+1]])
                ax.add_patch(a.longpath_patch())
            self.distances = - self.distances
    
#%%

X = DAG(20)
X.draw(shortest=True,longest=True)
X.show_shortest()
X.show_longest()

#%%