import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs

def DumbbellBosons(N):

    D11 = cls.TriangleCell(N)
    D22 = cls.TriangleCell(N)
    D33 = cls.TriangleCell(N)
    D21 = cls.SquareCell(N)
    D31 = cls.SquareCell(N)
    D32 = cls.SquareCell(N)

    D11.indices = (1,1)
    D22.indices = (2,2)
    D33.indices = (3,3)
    D21.indices = (2,1)
    D31.indices = (3,1)
    D32.indices = (3,2)

    C = configs.ConfigurationSpace([D11,D22,D33,D21,D31,D32])

    gluing1 = [D11.y0,D11.x1,D21.x0]
    gluing2 = [D21.y0,D21.y1,D22.y0]
    gluing3 = [D21.x1,D31.x0,D31.x1]
    gluing4 = [D31.y0,D31.y1,D32.y0]
    gluing5 = [D22.x1,D32.x0,D32.x1]
    gluing6 = [D32.y1,D33.y0,D33.x1]
    C.glue(gluing1)
    C.glue(gluing2)
    C.glue(gluing3)
    C.glue(gluing4)
    C.glue(gluing5)
    C.glue(gluing6)

    C.exterior_bc("dirichlet")
    C.diagonal_bc("dirichlet")

    C.gen_lapl()

    return C

if __name__ == "__main__":
    # Main

    N = 30
    h = (np.pi)/(N-1)

    C = DumbbellBosons(N)

    #CY.print_eqs()

    C.lapl_solve(h,2,50)
    spec = C.spectrum
    spec.sort()
    print(spec)
    C.plot_states(0)
    #C.plot_states(1)
    #C.plot_states(2)