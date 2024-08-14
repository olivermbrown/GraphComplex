import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs

def YgraphBosons(N):

    D11 = cls.TriangleCell(N)
    D22 = cls.TriangleCell(N)
    D33 = cls.TriangleCell(N)
    D12 = cls.SquareCell(N)
    D13 = cls.SquareCell(N)
    D23 = cls.SquareCell(N)

    D11.indices = (1,1)
    D22.indices = (2,2)
    D33.indices = (3,3)
    D12.indices = (1,2)
    D13.indices = (1,3)
    D23.indices = (2,3)

    CY = configs.ConfigurationSpace([D11,D22,D33,D12,D13,D23])

    gluing1 = [D11.y0, D13.x0, D12.x0]
    gluing2 = [D22.y0, D12.y0, D23.x0]
    gluing3 = [D33.y0, D23.y0, D13.y0]

    CY.glue(gluing1)
    CY.glue(gluing2)
    CY.glue(gluing3)

    CY.exterior_bc("dirichlet")
    CY.diagonal_bc("dirichlet")

    CY.gen_lapl()

    return CY

if __name__ == "__main__":
    # Main

    N = 50
    h = (np.pi)/(N-1)

    CY = YgraphBosons(N)

    CY.lapl_solve(h,2,20)
    spec = CY.spectrum
    spec.sort()
    print(spec)
    CY.plot_states(0)