import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs

def CrossgraphBosons(N):

    D11 = cls.TriangleCell(N)
    D22 = cls.TriangleCell(N)
    D33 = cls.TriangleCell(N)
    D44 = cls.TriangleCell(N)
    D12 = cls.SquareCell(N)
    D13 = cls.SquareCell(N)
    D14 = cls.SquareCell(N)
    D23 = cls.SquareCell(N)
    D24 = cls.SquareCell(N)
    D34 = cls.SquareCell(N)

    D11.indices = (1,1)
    D22.indices = (2,2)
    D33.indices = (3,3)
    D44.indices = (4,4)
    D12.indices = (1,2)
    D13.indices = (1,3)
    D14.indices = (1,4)
    D23.indices = (2,3)
    D24.indices = (2,4)
    D34.indices = (3,4)

    CY = configs.ConfigurationSpace([D11,D22,D33,D44,D12,D13,D14,D23,D24,D34])

    gluing1 = [D11.x0, D12.y0, D13.y0, D14.y0]
    gluing2 = [D22.x0, D12.x0, D23.y0, D24.y0]
    gluing3 = [D33.x0, D13.x0, D23.x0, D34.y0]
    gluing4 = [D44.x0, D14.x0, D24.x0, D34.x0]

    CY.glue(gluing1)
    CY.glue(gluing2)
    CY.glue(gluing3)
    CY.glue(gluing4)

    CY.exterior_bc("dirichlet")
    CY.diagonal_bc("neumann")

    CY.gen_lapl()

    return CY

def CrossgraphAnyonsHardcore(N, alpha):

    D11 = cls.TriangleCell(N)
    D22 = cls.TriangleCell(N)
    D33 = cls.TriangleCell(N)
    D44 = cls.TriangleCell(N)
    D12 = cls.SquareCell(N)
    D13 = cls.SquareCell(N)
    D14 = cls.SquareCell(N)
    D23 = cls.SquareCell(N)
    D24 = cls.SquareCell(N)
    D34 = cls.SquareCell(N)

    D11.indices = (1,1)
    D22.indices = (2,2)
    D33.indices = (3,3)
    D44.indices = (4,4)
    D12.indices = (1,2)
    D13.indices = (1,3)
    D14.indices = (1,4)
    D23.indices = (2,3)
    D24.indices = (2,4)
    D34.indices = (3,4)

    CY = configs.ConfigurationSpace([D11,D22,D33,D44,D12,D13,D14,D23,D24,D34])

    gluing1 = [D11.x0, D12.y0, D13.y0, D14.y0]
    gluing2 = [D22.x0, D12.x0, D24.y0, D23.y0]
    gluing3 = [D33.x0, D13.x0, D23.x0, D34.y0]
    gluing4 = [D44.x0, D14.x0, D34.x0, D24.x0]

    phase = np.exp(1j*np.pi*alpha)

    CY.glue(gluing1)
    CY.glue_with_branch_cut(gluing2, phase)
    CY.glue_with_branch_cut(gluing3, phase)
    CY.glue_with_branch_cut(gluing4, phase)

    CY.exterior_bc("dirichlet")
    CY.diagonal_bc("dirichlet")

    CY.gen_lapl()

    return CY

if __name__ == "__main__":
    # Main

    N = 100
    h = (np.pi)/(N-1)
    alpha = 1.0

    C = CrossgraphBosons(N)

    #C = CrossgraphAnyonsHardcore(N,alpha)

    C.lapl_solve(h,20)
    spec = C.spectrum
    spec.sort()
    print(spec)
    C.plot_states(0)