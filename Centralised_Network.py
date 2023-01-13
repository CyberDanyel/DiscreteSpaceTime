'''
                        GANGA SINGH MANCHANDA
                        Imperial College London
                        January 2023
'''
#%%

from Network_Elements import *

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
        self.n = n
        self.nodes = []
        self.edges = []
        
        for i in range(self.n):
            x = np.random.uniform(0,1)
            y = np.random.uniform(0,1)
            self.nodes.append(Node(position=np.array([x,y]),vertex=False))
                        
        for i in range(1,self.n):
            self.edges.append(Edge(self.nodes[0],self.nodes[i]))

        with plt.style.context('dark_background'):
            plt.figure(dpi=540)
            ax = plt.axes(xlim=(0,1),ylim=(0,1))
            ax.set_aspect('equal',adjustable='box')
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            for i in range(self.n):
                ax.add_patch(self.nodes[i].patch())
            for i in range(self.n - 1):
                ax.add_patch(self.edges[i].edge_patch())
            plt.show()
            
#%%

Centralised_Network(50)

#%%
