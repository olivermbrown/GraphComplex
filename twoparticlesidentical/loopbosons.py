import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs

def LoopBosons(N):

    D11 = cls.TriangleCell(N)

    D11.indices = (1,1)
    

    CL = configs.ConfigurationSpace([D11])

    gluing = [D11.y0, D11.x1]
    #gluing2 = [D11.y1, D11.x0]

    CL.glue(gluing)

    CL.exterior_bc("dirichlet")
    CL.diagonal_bc("dirichlet")

    CL.gen_lapl()

    return CL

if __name__ == "__main__":
    # Main

    N = 10
    h = (np.pi)/(N-1)

    TG = LoopBosons(N)

    TG.print_eqs()

    TG.lapl_solve(h,2,10)
    spec = TG.spectrum
    spec.sort()
    print(spec)
    TG.plot_states(0)
    TG.plot_states(1)
    TG.plot_states(2)
    pass