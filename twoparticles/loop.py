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

def looptwoparticles(N):
      # Create the network with two particles on a single wire
      
      # Scaling factor
      h = (np.pi)/(N-1)
      
      # Create the domain
      D11 = squareFDM.Domain(N)
      D11.split_domain()

      D11.indices = (1,1)
      
      cells = [D11]
      
      Network = twocomplex.Complex(cells)

      gluingA = [D11.x0, D11.x1]
      gluingB = [D11.y0, D11.y1]

      Network.glue(gluingA)
      Network.glue(gluingB)
      
      Network.diagonal_bc("dirichlet")
      #Network.exterior_bc("dirichlet")
      
      Network.gen_lapl()
      
      return Network

if __name__ == "__main__":
    # Main
    
    N = 100 # The length of a wire
    
    OneWire = looptwoparticles(N)

    h = (np.pi)/(N-1)

    #OneWire.print_eqs()

    OneWire.lapl_solve(h,2,20)
    print(OneWire.spectrum)
    s = sorted(OneWire.spectrum)
    print(s)

    # Plot the states
    OneWire.plot_states(0)
    OneWire.plot_states(1)
    OneWire.plot_states(2)
    pass