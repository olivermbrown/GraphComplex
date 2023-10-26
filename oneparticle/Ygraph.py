"""
Code to solve the Laplacian on a two glued wires using the finite difference method.
"""

import matplotlib.pyplot as plt
from matplotlib import cm
import lineFDM
import onecomplex

def oneparticleYgraph(N):
    # Create a 1 complex describing one particle on a single wire

    # Create one particle domains
    l1 = lineFDM.Line(N)
    l2 = lineFDM.Line(N)
    l3 = lineFDM.Line(N)

    # Create one particle complex
    cells = [l1, l2, l3]
    Network = onecomplex.Complex(cells)

    # Set boundary conditions
    Network.exterior_bc("dirichlet")

    # Glue the two wires together
    gluing = [l1.start, l2.start, l3.start]
    Network.glue(gluing)

    Network.gen_lapl()

    return Network

if __name__=="__main__":
    # Main
    
    # Define length of wire
    N = 50

    Y = oneparticleYgraph(N)

    Y.print_eqs()

    spectrum, states = Y.lapl_solve(2)

    spectrum.sort()
    print(spectrum)

    onecomplex.plot_state(states[0])

    pass