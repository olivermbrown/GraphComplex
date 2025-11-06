import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs

def TwoWires(N):

    D11 = cls.TriangleCell(N)
    D21 = cls.SquareCell(N)
    D22 = cls.TriangleCell(N)

    D11.indices = (1,1)
    D21.indices = (2,1)
    D22.indices = (2,2)

    C = configs.ConfigurationSpace([D11,D21,D22])

    gluing1 = [D11.x1,D21.x0]
    gluing2 = [D21.y1,D22.y0]
    C.glue(gluing1)
    C.glue(gluing2)

    C.exterior_bc("dirichlet")
    C.diagonal_bc("dirichlet")

    C.gen_lapl()

    return C

def TwoWiresDifferentLengths(N,N2):

    D11 = cls.TriangleCell(N)
    D21 = cls.SquareCell(N2,N)
    D22 = cls.TriangleCell(N2)

    D11.indices = (1,1)
    D21.indices = (2,1)
    D22.indices = (2,2)

    C = configs.ConfigurationSpace([D11,D21,D22])

    gluing1 = [D11.y1,D21.y0]
    gluing2 = [D21.x1,D22.x0]
    C.glue(gluing1)
    C.glue(gluing2)

    C.exterior_bc("dirichlet")
    C.diagonal_bc("dirichlet")

    C.gen_lapl()

    return C

def TwoWiresWithPhases(N, alpha1, alpha2):

    D11 = cls.TriangleCell(N)
    D12 = cls.SquareCell(N)
    D22 = cls.TriangleCell(N)

    D11.indices = (1,1)
    D12.indices = (2,1)
    D22.indices = (2,2)

    C = configs.ConfigurationSpace([D11,D12,D22])

    gluing1 = [D11.y1,D12.y0]
    gluing2 = [D12.x1,D22.x0]

    phase1 = np.exp(1j*np.pi*alpha1)
    phase2 = np.exp(1j*np.pi*alpha2)

    C.glue_with_branch_cut(gluing1, phase1)
    C.glue_with_branch_cut(gluing2, phase2)

    C.exterior_bc("dirichlet")
    C.diagonal_bc("dirichlet")

    return C

def ArrangeTwoWiresPlots(C):

    C.plot_dim = (2,2)

    C.figuresize = (8,6)

    for cell in C.cells:
        if cell.indices == (1,1):
            cell.plot_loc = (1,0)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (1,2):
            cell.plot_loc = (0,0)
            cell.use_x_labels = False
            cell.use_y_labels = True
            pass
        elif cell.indices == (2,2):
            cell.plot_loc = (0,1)
            cell.use_x_labels = True
            cell.use_y_labels = False
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
    alpha1 = 0.0
    alpha2 = 0.0

    C = TwoWiresWithPhases(N, alpha1, alpha2)

    #C = TwoWiresDifferentLengths(N,N2)

    #CY.print_eqs()

    C.lapl_solve(h,2)
    spec = C.spectrum
    spec.sort()
    print(spec)
    ArrangeTwoWiresPlots(C)
    C.plot_states(0, plotting_method="contour", realimag="abs", N_levels=20)