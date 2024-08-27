"""
Code to solve the Laplacian on a the lasso graph using the finite difference method.
"""

import matplotlib.pyplot as plt
from matplotlib import cm
import lineFDM
import onecomplex

def oneparticletwowires(N):
    # Create a 1 complex describing one particle on a single wire

    # Create one particle domain
    l1 = lineFDM.Line(N)
    l2 = lineFDM.Line(N)
    l3 = lineFDM.Line(N)

    # Create one particle complex
    cells = [l1, l2, l3]
    Network = onecomplex.Complex(cells)

    # Set boundary conditions
    Network.exterior_bc("dirichlet")

    # Glue the two wires together
    gluing1 = [l1.start, l1.end, l2.start]
    gluing2 = [l2.end, l3.start, l3.end]
    Network.glue(gluing1)
    Network.glue(gluing2)

    Network.gen_lapl()

    return Network

if __name__=="__main__":
    # Main
    
    # Define length of wire
    N = 100

    twowires = oneparticletwowires(N)

    spectrum, states = twowires.lapl_solve(2,n_eigs=100)

    #spectrum.sort()
    print(spectrum)

    twowires.plot_states(3)

    pass