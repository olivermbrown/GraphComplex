import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs

def DumbbellAnyonsHardcore(N, alpha1, alpha2):

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

    C = configs.ConfigurationSpace([D11,D22,D33,D12,D13,D23])

    gluing1 = [D11.y1,D12.y0,D11.x0]
    gluing2 = [D12.x0,D12.x1,D22.x0]
    gluing3 = [D12.y1,D13.y0,D13.y1]
    gluing4 = [D13.x0,D13.x1,D23.x0]
    gluing5 = [D22.y1,D23.y0,D23.y1]
    gluing6 = [D23.x1,D33.x0,D33.y1]

    phase1 = np.exp(1j*np.pi*alpha1)
    phase2 = np.exp(1j*np.pi*alpha2)

    C.glue_with_branch_cut(gluing1, phase1)
    C.glue(gluing2)
    C.glue(gluing3)
    C.glue(gluing4)
    C.glue(gluing5)
    C.glue_with_branch_cut(gluing6, phase2)

    C.exterior_bc("dirichlet")
    C.diagonal_bc("dirichlet")

    return C

def DumbbellAnyonsNeumann(N, alpha1, alpha2):

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

    C = configs.ConfigurationSpace([D11,D22,D33,D12,D13,D23])

    gluing1 = [D11.y1,D12.y0,D11.x0]
    gluing2 = [D12.x0,D12.x1,D22.x0]
    gluing3 = [D12.y1,D13.y0,D13.y1]
    gluing4 = [D13.x0,D13.x1,D23.x0]
    gluing5 = [D22.y1,D23.y0,D23.y1]
    gluing6 = [D23.x1,D33.x0,D33.y1]

    phase1 = np.exp(1j*np.pi*alpha1)
    phase2 = np.exp(1j*np.pi*alpha2)

    C.glue_with_branch_cut(gluing1, phase1)
    C.glue(gluing2)
    C.glue(gluing3)
    C.glue(gluing4)
    C.glue(gluing5)
    C.glue_with_branch_cut(gluing6, phase2)

    C.exterior_bc("dirichlet")
    C.diagonal_bc("neumann")

    return C

def ArrangeDumbbellPlots(C):

    C.plot_dim = (3,3)

    C.figuresize = (8,6)

    for cell in C.cells:
        if cell.indices == (1,1):
            cell.plot_loc = (2,0)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (1,2):
            cell.plot_loc = (1,0)
            cell.use_x_labels = False
            cell.use_y_labels = True
            pass
        elif cell.indices == (1,3):
            cell.plot_loc = (0,0)
            cell.use_x_labels = False
            cell.use_y_labels = True
            pass
        elif cell.indices == (2,2):
            cell.plot_loc = (1,1)
            cell.use_x_labels = True
            cell.use_y_labels = False
            pass
        elif cell.indices == (2,3):
            cell.plot_loc = (0,1)
            cell.use_x_labels = False
            cell.use_y_labels = False
            pass
        elif cell.indices == (3,3):
            cell.plot_loc = (0,2)
            cell.use_x_labels = True
            cell.use_y_labels = False
            pass
        else:
            pass
        pass

    C.plots_off = [(2,1),(2,2),(1,2)]

    return None

if __name__ == "__main__":
    # Main

    N = 120
    h = (np.pi)/(N-1)
    alpha1 = 0.0
    alpha2 = 0.0
    
    C = DumbbellAnyonsNeumann(N,alpha1,alpha2)

    C.gen_lapl()

    C.lapl_solve(h,50)

    # Load the states
    # Eigenvalues filepath
    states_path = "dumbbell_hardcore_states/dumbbell_N" + str(N) + "_alpha1" + str(alpha1) + "_alpha2" + str(alpha2) + "_"
    # Eigenstates filepath
    eigs_path = "dumbbell_hardcore_eigenvalues/dumbbell_N" + str(N) + "_alpha1" + str(alpha1) + "_alpha2" + str(alpha2)

    # Load the eigenvalues
    #C.load_eigenvalues(eigs_path)
    # Load the eigenstates
    #C.load_states(states_path)

    spec = C.spectrum
    spec.sort()
    print(spec)
    ArrangeDumbbellPlots(C)
    C.plot_states(0, plotting_method="contour", realimag="phase", N_levels=20)
    #C.plot_states(1)
    #C.plot_states(2)