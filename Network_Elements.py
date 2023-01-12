'''
                        GANGA SINGH MANCHANDA
                        Imperial College London
                        January 2023
'''
#%%

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, ArrowStyle

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
    get_patch : returns patch
    
    '''
    
    def __init__(self,position=np.array([0.0,0.0]),vertex=False):
        self._r = position
        self._vertex = vertex
        self._patch1 = plt.Circle(self._r,0.01,fc='white')
        self._patch2 = plt.Circle(self._r,0.05,fc='red')

    def pos(self):
        return self._r

    def get_patch(self):
        if self._vertex == False:
            return self._patch1
        else:
            return self._patch2
    
#%%

class Edge:
    
    '''
    
    Docstring
    
    '''
    
    def __init__(self,start=Node(position=np.array([0.0,0.0])),
                      end=Node(position=np.array([0.0,0.0])),diag=False):
        self._start = start.pos()
        self._end = end.pos()
        self._diag = diag

        style = ArrowStyle("Fancy", head_length=5, head_width=5, tail_width=0.5)
        self._patch = FancyArrowPatch(self._start,self._end,arrowstyle=style,color='purple')
        
    def get_patch(self):
            return self._patch