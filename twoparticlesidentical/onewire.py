import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs

def Wire(N):

    D11 = cls.TriangleCell(N)

    D11.indices = (1,1)
    

    CW = configs.ConfigurationSpace([D11])

    CW.exterior_bc("dirichlet")
    CW.diagonal_bc("dirichlet")

    CW.gen_lapl()

    return CW

if __name__ == "__main__":
    # Main

    N = 50
    h = (np.pi)/(N-1)

    CW = Wire(N)

    #CY.print_eqs()

    CW.lapl_solve(h,2,30)
    spec = CW.spectrum
    spec.sort()
    print(spec)
    CW.plot_states(0)