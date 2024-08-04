"""
Code to solve the Laplacian on the distinguishable configuration space of two particles on the lasso graph.
"""

import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import squareFDM
import twocomplex

def lassotwoparticles(N):
    # Create the network with two particles on the lasso graph

    # Scaling factor
    h = (np.pi)/(N-1)

    D11 = squareFDM.Domain(N)
    D12 = squareFDM.Domain(N)
    D21 = squareFDM.Domain(N)
    D22 = squareFDM.Domain(N)
    D11.split_domain()
    D22.split_domain()
    D11.indices = (1,1)
    D12.indices = (1,2)
    D21.indices = (2,1)
    D22.indices = (2,2)

    cells = [D11,D12,D21,D22]

    Network = twocomplex.Complex(cells)

    gluingA = [D11.x1, D12.x0, D12.x1]
    gluingB = [D12.y1, D22.y0, D22.y1]
    gluingC = [D22.x0, D21.x1, D22.x1]
    gluingD = [D21.y0, D11.y1, D21.y1]

    Network.glue(gluingA)
    Network.glue(gluingB)
    Network.glue(gluingC)
    Network.glue(gluingD)

    Network.diagonal_bc("dirichlet")
    Network.exterior_bc("dirichlet")

    Network.gen_lapl()
    
    return Network

if __name__ == "__main__":
    # Main
    
    N = 50 # The length of a wire
    
    # Scaling factor
    h = (np.pi)/(N-1)

    Lasso = lassotwoparticles(N)

    #TwoWires.print_eqs()

    #spectrum = TwoWires.lapl_spectrum(h,2,20)
    #print(spectrum)
    #print("Scaling factor:")
    #print(5/spectrum[0])
    #print(4*spectrum)

    Lasso.lapl_solve(h,2,50)
    print(Lasso.spectrum)

    spectrum = Lasso.spectrum
    print("Scaling factor:")
    print(18/spectrum[1])
    print(18/spectrum[1]*spectrum)

    # Plot the states
    Lasso.plot_states(1)
    
    pass