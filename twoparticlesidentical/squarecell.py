import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs

def Square(N, alpha):

    D11 = cls.SquareCell(N)

    D11.indices = (1,1)
    
    CW = configs.ConfigurationSpace([D11])

    gluing = [D11.x0, D11.y0]

    phase = np.exp(1j*np.pi*alpha)

    CW.glue_with_branch_cut(gluing, phase)

    CW.exterior_bc("dirichlet")
    CW.diagonal_bc("dirichlet")

    return CW

def ArrangeYgraphPlots(C):

    C.plot_dim = (1,1)

    C.figuresize = (19,9)

    for cell in C.cells:
        if cell.indices == (1,1):
            cell.plot_loc = (0,0)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        else:
            pass
        pass

    return None

if __name__ == "__main__":
    # Main

    N = 50
    N2 = 50
    h = (np.pi)/(N-1)

    CW = Square(N)

    #CW.print_eqs()

    CW.lapl_solve(h, N_eigs=20)
    spec = CW.spectrum
    spec.sort()
    # Round spectrum to 2 decimal places
    spec = np.round(spec, 2)
    print(spec)
    CW.plot_states(0)