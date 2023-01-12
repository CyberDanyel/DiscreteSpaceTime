'''
                        GANGA SINGH MANCHANDA
                        Imperial College London
                        January 2023
'''
#%%

from NetworkElements import *

#%%

class Centralised_Network:
    
    '''
    
    This class defines a centralised network object.
    
    Attributes
    ----------
    n : integer
        number of nodes
    nodes : list
            list of node objects
    edges : list
            list of edge objects
    
    '''

    def __init__(self,n):
        self._n = n
        self._nodes = []
        self._edges = []
        
        for i in range(self._n):
            x = np.random.uniform(0,1)
            y = np.random.uniform(0,1)
            self._nodes.append(Node(position=np.array([x,y])),vertex=False)
            
# Nodes are distributed by a Poisson Point Process
            
        for i in range(1,self._n):
            self._edges.append(Edge(self._nodes[0],self._nodes[i]))

        with plt.style.context('ggplot'):
            plt.figure(dpi=540)
            ax = plt.axes(xlim=(0,1),ylim=(0,1))
            for i in range(self._n):
                ax.add_patch(self._nodes[i].get_patch())
            for i in range(self._n - 1):
                ax.add_patch(self._edges[i].get_patch())
            plt.show()
            
# A complete graph contains N(N-1)/2 Edges

