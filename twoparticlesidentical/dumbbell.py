import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs

def DumbbellAnyons(N, alpha):

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

    gluing1 = [D11.y1,D21.x0,D11.x0]
    gluing2 = [D21.y0,D21.y1,D22.x0]
    gluing3 = [D21.x1,D31.x0,D31.x1]
    gluing4 = [D31.y0,D31.y1,D32.y0]
    gluing5 = [D22.y1,D32.x0,D32.x1]
    gluing6 = [D32.y1,D33.x0,D33.y1]

    phase = np.exp(1j*np.pi*alpha)

    C.glue_with_branch_cut(gluing1, phase)
    C.glue(gluing2)
    C.glue(gluing3)
    C.glue(gluing4)
    C.glue(gluing5)
    C.glue_with_branch_cut(gluing6, phase)

    C.exterior_bc("dirichlet")
    C.diagonal_bc("dirichlet")

    return C

if __name__ == "__main__":
    # Main

    N = 70
    h = (np.pi)/(N-1)
    alpha = 0.0

    C = DumbbellAnyons(N,alpha)

    #CY.print_eqs()

    #C.gen_lapl()

    #C.lapl_solve(h,50)

    # Load the states
    # Eigenvalues filepath
    states_path = "dumbbell_states/dumbbell_N" + str(N) + "_alpha" + str(alpha) + "_"
    # Eigenstates filepath
    eigs_path = "dumbbell_eigenvalues/dumbbell_N" + str(N) + "_alpha" + str(alpha)

    # Load the eigenvalues
    C.load_eigenvalues(eigs_path)
    # Load the eigenstates
    C.load_states(states_path)

    spec = C.spectrum
    spec.sort()
    print(spec)
    C.plot_states(0)
    #C.plot_states(1)
    #C.plot_states(2)