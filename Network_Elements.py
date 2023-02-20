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
from math import inf

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
    
    def __init__(self,position=np.array([0.0,0.0]),vertex=False, identifier = None, visited = False):
        self.r = position
        self.vertex = vertex
        self.nodepatch = plt.Circle(self.r,0.01,fc='white')
        self.vertexpatch = plt.Circle(self.r,0.05,fc='white')
        
        self.source = False
        self.sink = False
        self.identifier = identifier
        self.connected_nodes = []
        self.tentative_distance = inf
        self.visited = visited
        self.previous_vertex = None
        
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
    
    def connect(self,connected_node):
        self.connected_nodes.append(connected_node)
        
    def neighbours(self):
        return self.connected_nodes
    
    def set_identifier(self,new_identifier):
        self.identifier = new_identifier
        
    def identifier(self):
        return self.identifier
    
    def tentative_distance(self):
        return self.tentative_distance
    
    def set_tentative_distance(self,new_tentative_distance):
        self.tentative_distance = new_tentative_distance
        
    def returnvisited(self):
        return self.visited
    
    def mark_visited(self):
        self.visited = True
        
    def mark_unvisited(self):
        self.visited = False
        
    def set_previous_vertex(self, new_previous_vertex):
        self.previous_vertex = new_previous_vertex
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
        
    def dist(self,p_values):
        self.distance = []
        for i in range(len(p_values)):
            if p_values[i] == 2:
                self.distance.append(np.linalg.norm(self.start-self.end))
            else:
                self.distance.append(((self.end[0] - self.start[0])**p_values[i]+(self.end[1] - self.start[1])**p_values[i])**(1/p_values[i]))
        
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