"""
Code to solve the Laplacian on a two glued wires using the finite difference method.
"""

import matplotlib.pyplot as plt
from matplotlib import cm
import lineFDM
import onecomplex

def oneparticleHgraph(N):
    # Create a 1 complex describing one particle on a single wire

    # Create one particle domains
    l1 = lineFDM.Line(N)
    l2 = lineFDM.Line(N)
    l3 = lineFDM.Line(N)
    l4 = lineFDM.Line(N)
    l5 = lineFDM.Line(N)

    # Create one particle complex
    cells = [l1, l2, l3, l4, l5]
    Network = onecomplex.Complex(cells)

    # Set boundary conditions
    Network.exterior_bc("dirichlet")

    # Glue the two wires together
    gluing1 = [l1.start, l2.start, l3.start]
    gluing2 = [l3.end, l4.start, l5.start]
    Network.glue(gluing1)
    Network.glue(gluing2)

    Network.gen_lapl()

    return Network

if __name__=="__main__":
    # Main
    
    # Define length of wire
    N = 50

    Y = oneparticleHgraph(N)

    Y.print_eqs()

    spectrum, states = Y.lapl_solve(2,30)

    spectrum.sort()
    print(spectrum)

    Y.plot_states(0)

    pass