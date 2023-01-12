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
    
    Docstring
    
    '''
    
# We will generate random triangular adjacency matrices to represent DAGs.
# Further research must be done to determine whether this will excessively
# restrict the type of networks we can produce.
    
    def __init__(self,n):
        self._n = n
        self._nodes = []
        self._edges = []
        
        M = np.zeros((self._n,self._n),dtype = int)
        for i in range(0,self._n):
            M[i] = np.random.randint(2,size=self._n)
        self._Adjacency = np.triu(M,-1)
        print (self._Adjacency)
        n_edges = np.sum(self._Adjacency)
                
        self._nodes.append(Node(position=np.array([0,0]),vertex=True))
        for i in range(self._n - 2):
            x = np.random.uniform(0,1)
            y = np.random.uniform(0,1)
            self._nodes.append(Node(position=np.array([x,y]),vertex=False))
        self._nodes.append(Node(position=np.array([1,1]),vertex=True))

        x = np.asarray(np.where(self._Adjacency == 1))
        for i in range(n_edges):
            self._edges.append(Edge(self._nodes[x[0][i]],self._nodes[x[1][i]]))
            
        with plt.style.context('dark_background'):
            plt.figure(dpi=540)
            ax = plt.axes(xlim=(0,1),ylim=(0,1))
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            for i in range(self._n):
                ax.add_patch(self._nodes[i].get_patch())
            for i in range(n_edges):
                ax.add_patch(self._edges[i].get_patch())
            t = np.linspace(0,1,2)
            plt.plot(t,t,color='blue')
            plt.show()


    def djikstra(self):
        """
        https://github.com/crixodia/python-dijkstra/blob/master/dijkstra.py
        Returns a tuple with a distances' list and paths' list of
        all remaining vertices with the same indexing.
            (distances, paths)
        For example, distances[x] are the shortest distances from x
        vertex which shortest path is paths[x]. x is an element of
        {0, 1, ..., n-1} where n is the number of vertices
        Args:
        wmat    --  weighted graph's adjacency matrix
        start   --  paths' first vertex
        end     --  (optional) path's end vertex. Return just the 
                    distance and its path
        Exceptions:
        Index out of range, Be careful with start and end vertices
        """
        wmat = self._Adjacency
        start = 0
        end = self._n - 1
        
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

        return (dist[end], path[end])