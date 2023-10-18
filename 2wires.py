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

if __name__ == "__main__":
    # Main
    
    N = 50 # The length of a wire
    
    # Scaling factor
    h = (np.pi)/(N-1)
    
    D11 = squareFDM.Domain(N)
    D12 = squareFDM.Domain(N)
    D21 = squareFDM.Domain(N)
    D22 = squareFDM.Domain(N)
    D11.split_domain("dirichlet")
    D22.split_domain("dirichlet")
    
    cells = [D11,D12,D21,D22]
    
    Network = twocomplex.Complex(cells)
    
    gluingA = [(D11, D11.x1),(D12, D12.x0)]
    gluingB = [(D12, D12.y1),(D22, D22.y0)]
    gluingC = [(D22, D22.x0),(D21, D21.x1)]
    gluingD = [(D21, D21.y0),(D11, D11.y1)]
    
    Network.diagonal_bc("dirichlet")
    
    Network.glue(gluingA)
    Network.glue(gluingB)
    Network.glue(gluingC)
    Network.glue(gluingD)
    
    Network.exterior_bc("dirichlet")
    
    #print(Network.L)
    M = Network.L.toarray()
    #print(M)
    N = Network.simplify_lapl()
    new = Network.sL
    #print(new.toarray())
    spectrum = Network.lapl_spectrum(h,2,20)
    print(spectrum)
    print("Scaling factor:")
    print(5/spectrum[0])
    print(4*spectrum)
    
    spectrum, states = Network.lapl_solve(h,2)
    print(spectrum)
    n = 4
    for cell in cells:
            plot = Network.plot_state(cell, states[:,n])
            plot.show()
            pass
    
    pass