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

if __name__ == "__main__":
    # Main

    N = 20
    h = (np.pi)/(N-1)

    C = TwoWires(N)

    #CY.print_eqs()

    C.lapl_solve(h,2,50)
    spec = C.spectrum
    spec.sort()
    print(spec)
    C.plot_states(0)