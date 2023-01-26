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
        self.edges = []
        
        M = np.random.randint(2,size=(self.n,self.n))
        M[0][self.n-1] = 0
        #Test for correct distance:
        #M[0][self.n-1] = 1
        self.Adjacency = np.triu(M,+1)
        self.n_edges = np.sum(self.Adjacency)
        self.distances = np.zeros((self.n,self.n))
                
        self.nodes.append(Node(position=np.array([0,0]),vertex=True))
        for i in range(self.n - 2):
            x = np.random.uniform(0,1)
            y = np.random.uniform(0,1)
            self.nodes.append(Node(position=np.array([x,y]),vertex=False))
        self.nodes.append(Node(position=np.array([1,1]),vertex=True))

        x = np.asarray(np.where(self.Adjacency == 1))
        for i in range(self.n_edges):
            a = Edge(self.nodes[x[0][i]],self.nodes[x[1][i]])
            self.edges.append(a)
            self.distances[x[0][i]][x[1][i]] = a.dist()

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
        
    def draw(self, shortest = False, longest = False):
            
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
                    ax.add_patch(a.path_patch())
            if longest == True:
                self.distances = - self.distances
                self._djikstra()
                for i in range(len(self.path)-1):
                    a = Edge(self.nodes[self.path[i]],self.nodes[self.path[i+1]])
                    ax.add_patch(a.longpath_patch())
                self.distances = - self.distances
#%%

DAG1 = DAG(7)
DAG1.draw(shortest = True, longest = True)
DAG1.show_shortest()
DAG1.show_longest()
#%%