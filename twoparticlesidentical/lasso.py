import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs

def LassoAnyonsDirichlet(N, alpha):

    D11 = cls.TriangleCell(N)
    D12 = cls.SquareCell(N)
    D22 = cls.TriangleCell(N)

    D11.indices = (1,1)
    D12.indices = (1,2)
    D22.indices = (2,2)

    C = configs.ConfigurationSpace([D11,D12,D22])

    gluing1 = [D11.y1,D12.y0,D12.y1]
    gluing2 = [D12.x1,D22.x0,D22.y1]

    phase = np.exp(1j*np.pi*alpha)

    C.glue(gluing1)
    C.glue_with_branch_cut(gluing2, phase)

    C.exterior_bc("dirichlet")
    C.diagonal_bc("dirichlet")

    return C

def LassoAnyonsRobin(N, alpha):

    D11 = cls.TriangleCell(N)
    D12 = cls.SquareCell(N)
    D22 = cls.TriangleCell(N)

    D11.indices = (1,1)
    D12.indices = (1,2)
    D22.indices = (2,2)

    C = configs.ConfigurationSpace([D11,D12,D22])

    gluing1 = [D11.y1,D12.y0,D12.y1]
    gluing2 = [D12.x1,D22.x0,D22.y1]

    phase = np.exp(1j*np.pi*alpha)

    C.glue(gluing1)
    C.glue_with_branch_cut(gluing2, phase)

    C.exterior_bc("dirichlet")
    C.diagonal_bc("robin")

    return C

if __name__ == "__main__":
    # Main

    N = 80
    h = (np.pi)/(N-1)
    a = 0.999

    C = LassoAnyonsRobin(N,a)

    C.robin_constant = np.tan(np.pi*a/2)

    C.gen_lapl()

    #C.print_eqs()

    C.lapl_solve(h, N_eigs=10)
    spec = C.spectrum
    spec.sort()
    print(spec)
    C.plot_states(0, plotting_method="contour",realimag="real")

    pass