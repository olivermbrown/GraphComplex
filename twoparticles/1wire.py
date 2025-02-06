"""
Code to solve the Laplacian on the distinguishable configuration space of two particles on a single wire
"""

import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import squareFDM
import twocomplex

def onewiretwoparticles(N):
      # Create the network with two particles on a single wire
      
      # Scaling factor
      h = (np.pi)/(N-1)
      
      # Create the domain
      D11 = squareFDM.Domain(N)
      D11.split_domain()

      D11.indices = (1,1)
      
      cells = [D11]
      
      Network = twocomplex.Complex(cells)
      
      Network.diagonal_bc("dirichlet")
      Network.exterior_bc("dirichlet")
      
      Network.gen_lapl()
      
      return Network

if __name__ == "__main__":
    # Main
    
    N = 50 # The length of a wire
    
    OneWire = onewiretwoparticles(N)

    h = (np.pi)/(N-1)

    OneWire.lapl_solve(h,2,50)
    print(sorted(OneWire.spectrum))

    # Plot the states
    OneWire.plot_states(0)
    pass