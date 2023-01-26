'''
                        GANGA SINGH MANCHANDA
                        Imperial College London
                        January 2023
'''
#%%

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, ArrowStyle
import copy as cp

#%%

class Node:
    
    '''
    
    This class defines the node object.
    
    Attributes
    ----------
    r : numpy array
        position
    
    Methods
    -------
    pos : returns position
    patch : returns either a large(vertex) or small(node) patch
    
    '''
    
    def __init__(self,position=np.array([0.0,0.0]),vertex=False):
        self.r = position
        self.vertex = vertex
        self.nodepatch = plt.Circle(self.r,0.01,fc='white')
        self.vertexpatch = plt.Circle(self.r,0.05,fc='white')
        
        self.source = False
        self.sink = False
        
        if vertex == True:
            self.source = True
            self.sink = True

    def pos(self):
        return self.r

    def patch(self):
        if self.vertex == False:
            return self.nodepatch
        else:
            return self.vertexpatch
        
    def issource(self,cond=True):
        self.source = cond
        
    def issink(self,cond=True):
        self.sink = True
        
    def returnsource(self):
        return self.source
    
    def returnsink(self):
        return self.sink
    
#%%

class Edge:
    
    '''
    
    This class defines the edge object.
    
    Attributes
    ----------
    start : numpy array
            position of source node
    end : numpy array
          position of sink node
    dist : float
           distance between start and end or edge weight
    
    Methods
    -------
    dist : returns edge weight
    edge_patch : returns default edge patch
    path_patch : return special path patch
    
    '''
    
    def __init__(self,start=Node(position=np.array([0.0,0.0])),
                      end=Node(position=np.array([0.0,0.0]))):
        self.s = start
        self.e = end
        self.start = start.pos()
        self.end = end.pos()

        style = ArrowStyle("Fancy",head_length=5,head_width=5,tail_width=0.5)
        self.edgepatch = FancyArrowPatch(self.start,self.end,arrowstyle=style,color='dodgerblue')
        self.shortpathpatch = FancyArrowPatch(self.start,self.end,arrowstyle=style,color='green')
        self.longpathpatch = FancyArrowPatch(self.start,self.end,arrowstyle=style,color='red')
        
    def dist(self):
        self.distance = np.linalg.norm(self.start-self.end)
        return self.distance
        
    def edge_patch(self):
        return self.edgepatch
    
    def shortpath_patch(self):
        return self.shortpathpatch
    
    def longpath_patch(self):
        return self.longpathpatch
    
    def returnstart(self):
        return self.s
    
    def returnend(self):
        return self.e
    
#%%