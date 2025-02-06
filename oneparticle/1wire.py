"""
Code to solve the Laplacian on a single wire using the finite difference method.
"""

import matplotlib.pyplot as plt
from matplotlib import cm
import lineFDM
import onecomplex

def oneparticleonewire(N):
    # Create a 1 complex describing one particle on a single wire

    # Create one particle domain
    l1 = lineFDM.Line(N)

    # Create one particle complex
    cells = [l1]
    Network = onecomplex.Complex(cells)

    # Set boundary conditions
    Network.exterior_bc("dirichlet")

    Network.gen_lapl()

    return Network

if __name__=="__main__":
    # Main
    
    # Define length of wire
    N = 50

    onewire = oneparticleonewire(N)

    spectrum, states = onewire.lapl_solve(2)

    spectrum.sort()
    print(spectrum)

    pass