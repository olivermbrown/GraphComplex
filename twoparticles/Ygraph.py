"""
Code to solve the Laplacian on the distinguishable configuration space of two particles on the Y graph
"""

import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import squareFDM
import twocomplex
import configspace

def twoparticlesYgraph(N):
    # Evaluate the spectrum of the Laplacian on the distinguishable configuration space of two particles on the Y graph

    # Scaling factor
    h = (np.pi)/(N-1)

    # Create domains
    D11 = squareFDM.Domain(N)
    D22 = squareFDM.Domain(N)
    D33 = squareFDM.Domain(N)
    D12 = squareFDM.Domain(N)
    D13 = squareFDM.Domain(N)
    D23 = squareFDM.Domain(N)
    D21 = squareFDM.Domain(N)
    D31 = squareFDM.Domain(N)
    D32 = squareFDM.Domain(N)
    D11.split_domain()
    D22.split_domain()
    D33.split_domain()
    D11.indices = (1,1,0)
    D22.indices = (2,2,0)
    D33.indices = (3,3,0)
    D12.indices = (1,2,0)
    D13.indices = (1,3,0)
    D23.indices = (2,3,0)
    D21.indices = (2,1,0)
    D31.indices = (3,1,0)
    D32.indices = (3,2,0)
    
    cells = [D11,D22,D33,D12,D13,D23,D21,D31,D32]
    
    # Create two complex from domains
    Network = twocomplex.Complex(cells)
    
    # Glue the domains together
    gluingA = [D11.x0, D13.x0, D12.x0]
    gluingB = [D22.y0, D12.y0, D32.y0]
    gluingC = [D33.x0, D32.x0, D31.x0]
    gluingD = [D11.y0, D31.y0, D21.y0]
    gluingE = [D22.x0, D21.x0, D23.x0]
    gluingF = [D33.y0, D23.y0, D13.y0]
    
    Network.glue(gluingA)
    Network.glue(gluingB)
    Network.glue(gluingC)
    Network.glue(gluingD)
    Network.glue(gluingE)
    Network.glue(gluingF)

    # Apply boundary conditions to the exterior and diagonal of the complex
    Network.diagonal_bc("dirichlet")
    Network.exterior_bc("dirichlet")
    
    # Generate the discrete Laplacian
    Network.gen_lapl()
    
    return Network

if __name__ == "__main__":
    # Main
    
    N = 50 # The length of a wire

    # Scaling factor
    h = (np.pi)/(N-1)

    Ygraph = twoparticlesYgraph(N)

    #Ygraph.print_eqs()

    #spectrum = Ygraph.lapl_spectrum(h,2,50)
    #print(spectrum)

    Ygraph.lapl_solve(h,2,100)
    print("Unsorted eigenvalues = " + str(Ygraph.spectrum))
    print("Sorted eigenvalues = " + str(np.sort(Ygraph.spectrum)))

    # Plot the states
    Ygraph.plot_states(3)

    Configs = configspace.ConfigurationSpace(Ygraph)

    Configs.bosonic_projection()
    Configs.plot_projected_wavefunction(3)
    #Configs.save_projected_wavefunction(0)
    
    pass
