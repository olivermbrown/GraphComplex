"""
Code to solve the Laplacian on a single wire using the finite difference method.
"""

import matplotlib.pyplot as plt
from matplotlib import cm
import lineFDM
import onecomplex

def oneparticleloop(N):
    # Create a 1 complex describing one particle on a single wire

    # Create one particle domain
    l1 = lineFDM.Line(N)

    # Create one particle complex
    cells = [l1]
    Network = onecomplex.Complex(cells)

    # Set boundary conditions
    Network.exterior_bc("dirichlet")

    # Glue the ends of the wire together
    gluing = [l1.start, l1.end]
    Network.glue(gluing)

    Network.gen_lapl()

    return Network

if __name__=="__main__":
    # Main
    
    # Define length of wire
    N = 30

    loop = oneparticleloop(N)

    spectrum, states = loop.lapl_solve(2,n_eigs=10)

    print("Unsorted spectrum:")
    print(spectrum)
    spectrum.sort()
    print("Sorted spectrum:")
    print(spectrum)

    #loop.print_eqs()
    loop.plot_states(1)

    pass