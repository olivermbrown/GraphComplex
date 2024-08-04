"""
Code to solve the Laplacian on the distinguishable configuration space of two particles on a loop
"""

import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import squareFDM
import twocomplex

def looptwoparticles2cover(N):
      # Create the network with two particles on a single wire
      
      # Scaling factor
      h = (np.pi)/(N-1)
      
      # Create the domain
      D11_0 = squareFDM.Domain(N)
      D11_0.split_domain()
      D11_0.indices = (1,1,0)
      D11_1 = squareFDM.Domain(N)
      D11_1.split_domain()
      D11_1.indices = (1,1,1)
      
      cells = [D11_0, D11_1]
      
      Network = twocomplex.Complex(cells)

      gluingA = [D11_0.x0, D11_1.x1]
      gluingB = [D11_0.y0, D11_0.y1]
      gluingC = [D11_0.x1, D11_1.x0]
      gluingD = [D11_1.y0, D11_1.y1]

      Network.glue(gluingA)
      Network.glue(gluingB)
      Network.glue(gluingC)
      Network.glue(gluingD)

      Network.diagonal_bc("dirichlet")
      #Network.exterior_bc("dirichlet")
      
      Network.gen_lapl()
      
      return Network

if __name__ == "__main__":
    # Main
    
    N = 40 # The length of a wire
    
    Loop2Cover = looptwoparticles2cover(N)

    h = (np.pi)/(N-1)

    Loop2Cover.lapl_solve(h,2,50)
    print(Loop2Cover.spectrum)
    #spec = Loop2Cover.spectrum
    #scaling = 2/spec[0]
    #print(scaling*spec)

    # Plot the states
    Loop2Cover.plot_states(0)
    pass