import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs

def LoopAnyons(N, alpha):

    D11 = cls.TriangleCell(N)

    D11.indices = (1,1)
    
    CL = configs.ConfigurationSpace([D11])

    gluing = [D11.x0, D11.y1]

    phase = np.exp(1j*np.pi*alpha)

    CL.glue_with_branch_cut(gluing, phase)

    CL.exterior_bc("dirichlet")
    CL.diagonal_bc("dirichlet")

    CL.gen_lapl()

    return CL

def ArrangeLoopPlots(C):

    C.plot_dim = (1,1)

    C.figuresize = (8,6)

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

    N = 30
    h = (np.pi)/(N-1)
    a = 0.3

    TG = LoopAnyons(N,a)

    #TG.print_eqs()

    TG.lapl_solve(h,2)
    spec = TG.spectrum
    spec.sort()
    print(spec)
    ArrangeLoopPlots(TG)
    TG.plot_states(0, plotting_method="contour", realimag="phase", N_levels=40)
    pass