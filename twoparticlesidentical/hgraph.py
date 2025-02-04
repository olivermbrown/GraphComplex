import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs

def HgraphBosons(N):

    D11 = cls.TriangleCell(N)
    D22 = cls.TriangleCell(N)
    D33 = cls.TriangleCell(N)
    D44 = cls.TriangleCell(N)
    D55 = cls.TriangleCell(N)
    D12 = cls.SquareCell(N)
    D13 = cls.SquareCell(N)
    D14 = cls.SquareCell(N)
    D15 = cls.SquareCell(N)
    D23 = cls.SquareCell(N)
    D24 = cls.SquareCell(N)
    D25 = cls.SquareCell(N)
    D34 = cls.SquareCell(N)
    D35 = cls.SquareCell(N)
    D45 = cls.SquareCell(N)

    D11.indices = (1,1)
    D22.indices = (2,2)
    D33.indices = (3,3)
    D44.indices = (4,4)
    D55.indices = (5,5)
    D12.indices = (1,2)
    D13.indices = (1,3)
    D14.indices = (1,4)
    D15.indices = (1,5)
    D23.indices = (2,3)
    D24.indices = (2,4)
    D25.indices = (2,5)
    D34.indices = (3,4)
    D35.indices = (3,5)
    D45.indices = (4,5)

    CY = configs.ConfigurationSpace([D11,D22,D33,D44,D55,D12,D13,D14,D15,D23,D24,D25,D34,D35,D45])

    gluing1 = [D11.x0, D13.y0, D12.x0]
    gluing2 = [D22.x0, D12.y0, D23.y0]
    gluing3 = [D33.x0, D23.x0, D13.x0]

    gluing4 = [D33.y1, D34.y0, D35.y0]
    gluing5 = [D44.x0, D34.x1, D45.y0]
    gluing6 = [D55.x0, D35.x1, D45.x0]

    gluing7 = [D14.x0, D24.x0, D34.x0]
    gluing8 = [D15.x0, D25.x0, D35.x0]

    gluing9 = [D13.y1, D14.y0, D15.y0]
    gluing10 = [D23.y1, D24.y0, D25.y0]

    CY.glue(gluing1)
    CY.glue(gluing2)
    CY.glue(gluing3)
    CY.glue(gluing4)
    CY.glue(gluing5)
    CY.glue(gluing6)
    CY.glue(gluing7)
    CY.glue(gluing8)
    CY.glue(gluing9)
    CY.glue(gluing10)

    CY.exterior_bc("dirichlet")
    CY.diagonal_bc("neumann")

    CY.gen_lapl()

    return CY

def HgraphAnyons(N, alpha):

    D11 = cls.TriangleCell(N)
    D22 = cls.TriangleCell(N)
    D33 = cls.TriangleCell(N)
    D44 = cls.TriangleCell(N)
    D55 = cls.TriangleCell(N)
    D12 = cls.SquareCell(N)
    D13 = cls.SquareCell(N)
    D14 = cls.SquareCell(N)
    D15 = cls.SquareCell(N)
    D23 = cls.SquareCell(N)
    D24 = cls.SquareCell(N)
    D25 = cls.SquareCell(N)
    D34 = cls.SquareCell(N)
    D35 = cls.SquareCell(N)
    D45 = cls.SquareCell(N)

    D11.indices = (1,1)
    D22.indices = (2,2)
    D33.indices = (3,3)
    D44.indices = (4,4)
    D55.indices = (5,5)
    D12.indices = (1,2)
    D13.indices = (1,3)
    D14.indices = (1,4)
    D15.indices = (1,5)
    D23.indices = (2,3)
    D24.indices = (2,4)
    D25.indices = (2,5)
    D34.indices = (3,4)
    D35.indices = (3,5)
    D45.indices = (4,5)

    CY = configs.ConfigurationSpace([D11,D22,D33,D44,D55,D12,D13,D14,D15,D23,D24,D25,D34,D35,D45])

    gluing1 = [D11.x0, D12.y0, D13.y0]
    gluing2 = [D33.x0, D13.x0, D23.x0]
    gluing3 = [D22.x0, D23.y0, D12.x0]

    gluing4 = [D44.x0, D45.y0, D34.x1]
    gluing5 = [D33.y1, D34.y0, D35.y0]
    gluing6 = [D55.x0, D35.x1, D45.x0]

    gluing7 = [D14.x0, D24.x0, D34.x0]
    gluing8 = [D15.x0, D25.x0, D35.x0]

    gluing9 = [D13.y1, D14.y0, D15.y0]
    gluing10 = [D23.y1, D24.y0, D25.y0]

    phase = np.exp(1j*np.pi*alpha)

    CY.glue(gluing1)
    CY.glue(gluing2)
    CY.glue_with_branch_cut(gluing3, phase)
    CY.glue_with_branch_cut(gluing4, phase)
    CY.glue(gluing5)
    CY.glue(gluing6)
    CY.glue(gluing7)
    CY.glue(gluing8)
    CY.glue(gluing9)
    CY.glue(gluing10)

    CY.exterior_bc("dirichlet")
    CY.diagonal_bc("dirichlet")

    CY.gen_lapl()

    return CY

if __name__ == "__main__":
    # Main

    N = 70
    h = (np.pi)/(N-1)
    alpha = 1.0

    #CY = HgraphBosons(N)

    CY = HgraphAnyons(N,alpha)

    CY.lapl_solve(h,2)
    spec = CY.spectrum
    spec.sort()
    print(spec)
    CY.plot_states(0)