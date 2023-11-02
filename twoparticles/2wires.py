"""
Code to solve the Laplacian on the distinguishable configuration space of two particles on two wires
"""

import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import squareFDM
import twocomplex

def twowirestwoparticles(N):
    # Create the network with two particles on two connected wires

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

    gluingA = [D11.x1, D12.x0]
    gluingB = [D12.y1, D22.y0]
    gluingC = [D22.x0, D21.x1]
    gluingD = [D21.y0, D11.y1]

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
    
    N = 20 # The length of a wire
    
    # Scaling factor
    h = (np.pi)/(N-1)

    TwoWires = twowirestwoparticles(N)

    #TwoWires.print_eqs()

    #spectrum = TwoWires.lapl_spectrum(h,2,20)
    #print(spectrum)
    #print("Scaling factor:")
    #print(5/spectrum[0])
    #print(4*spectrum)

    TwoWires.lapl_solve(h,2,50)
    print(TwoWires.spectrum)

    # Plot the states
    TwoWires.plot_states(0)
    
    pass