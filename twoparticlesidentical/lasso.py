import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs

def LassoAnyonsHardcore(N, alpha):

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

def LassoAnyonsContact(N, alpha):

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

    C.alpha = alpha

    C.glue(gluing1)
    C.glue_with_branch_cut(gluing2, phase)

    C.exterior_bc("dirichlet")
    C.diagonal_bc("robin")

    return C

def ArrangeLassoPlots(C):

    C.plot_dim = (2,2)

    C.figuresize = (8,6)

    for cell in C.cells:
        if cell.indices == (1,1):
            cell.plot_loc = (1,0)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (1,2):
            cell.plot_loc = (0,0)
            cell.use_x_labels = False
            cell.use_y_labels = True
            pass
        elif cell.indices == (2,2):
            cell.plot_loc = (0,1)
            cell.use_x_labels = True
            cell.use_y_labels = False
            pass
        else:
            pass
        pass

    C.plots_off = [(1,1)]

    return None

if __name__ == "__main__":
    # Main

    N = 20
    h = (np.pi)/(N-1)
    alpha = 0.0

    #C = LassoAnyonsContact(N,alpha)

    #C.robin_constant = np.tan(np.pi*alpha/2)

    #C.robin_constant = 0



    C = LassoAnyonsHardcore(N,alpha)
    #C.robin_constant = 0

    # Eigenvalues filepath
    eigs_path = "lasso_hardcore_eigenvalues/lasso_N"+str(N)+"_alpha"+str(alpha)
    # Eigenstates filepath
    states_path = "lasso_hardcore_states/lasso_N"+str(N)+"_alpha"+str(alpha)+"_"

    # # Neumann Eigenvalues filepath
    # eigs_path = "lasso_neumann_eigenvalues/lasso_N"+str(N)+"_alpha"+str(alpha)
    # # Neumann Eigenstates filepath
    # states_path = "lasso_neumann_states/lasso_N"+str(N)+"_alpha"+str(alpha)+"_"

    # Load the eigenvalues
    #C.load_eigenvalues(eigs_path,override_directory_path=True)
    # Load the eigenstates
    #C.load_states(states_path,override_directory_path=True)

    C.gen_lapl()

    C.lapl_solve(h, N_eigs=50)

    # Save the eigenvalues
    #C.save_eigenvalues(eigs_path)
    # Save the eigenstates
    #C.save_states(states_path)

    spec = C.spectrum
    #spec.sort()
    print(spec)
    ArrangeLassoPlots(C)
    C.plot_states(0, plotting_method="surface", realimag="real", N_levels=20)

    pass