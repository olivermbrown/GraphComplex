import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs

def YgraphAnyonsHardcore(N, alpha):

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

    gluing1 = [D11.x0, D12.y0, D13.y0]
    gluing2 = [D33.x0, D13.x0, D23.x0]
    gluing3 = [D22.x0, D23.y0, D12.x0]

    phase = np.exp(1j*np.pi*alpha)

    CY.glue(gluing1)
    CY.glue(gluing2)
    CY.glue_with_branch_cut(gluing3, phase)

    CY.exterior_bc("dirichlet")
    CY.diagonal_bc("dirichlet")

    return CY

def YgraphAnyonsRobin(N, alpha):

    # Create domains
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

    gluing1 = [D11.x0, D12.y0, D13.y0]
    gluing2 = [D33.x0, D13.x0, D23.x0]
    gluing3 = [D22.x0, D23.y0, D12.x0]

    phase = np.exp(1j*np.pi*alpha)

    CY.glue(gluing1)
    CY.glue(gluing2)
    CY.glue_with_branch_cut(gluing3, phase)

    # Apply boundary conditions to the exterior and diagonal of the complex
    CY.diagonal_bc("robin")
    CY.exterior_bc("dirichlet")

    return CY

def ArrangeYgraphPlots(C):

    C.plot_dim = (2,3)

    C.figuresize = (19,9)

    for cell in C.cells:
        if cell.indices == (1,1):
            cell.plot_loc = (0,0)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (1,2):
            cell.plot_loc = (1,0)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (2,2):
            cell.plot_loc = (0,1)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (1,3):
            cell.plot_loc = (1,1)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (2,3):
            cell.plot_loc = (1,2)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (3,3):
            cell.plot_loc = (0,2)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        else:
            pass
        pass

    return None

if __name__ == "__main__":
    # Main

    N = 100
    h = (np.pi)/(N-1)
    alpha = 0.0

    CY = YgraphAnyonsHardcore(N,alpha)

    #CY.robin_constant = np.tan(np.pi*alpha/2)

    #CY.print_eqs()

    # Eigenvalues filepath
    eigs_path = "ygraph_hardcore_eigenvalues/ygraph_N"+str(N)+"_alpha"+str(alpha)
    # Eigenstates filepath
    states_path = "ygraph_hardcore_states/ygraph_N"+str(N)+"_alpha"+str(alpha)+"_"

    #CY = YgraphAnyonsHardcore(N,alpha)

    # Eigenvalues filepath
    #eigs_path = "ygraph_neumann_eigenvalues/ygraph_N"+str(N)+"_alpha"+str(alpha)
    # Eigenstates filepath
    #states_path = "ygraph_neumann_states/ygraph_N"+str(N)+"_alpha"+str(alpha)+"_"

    #CY = YgraphAnyonsRobin(N,alpha)
    #CY.robin_constant = 0

    #CY.gen_lapl()

    #CY.lapl_solve(h, N_eigs=50)

    # Save the eigenvalues
    #CY.save_eigenvalues(eigs_path)
    # Save the eigenstates
    #CY.save_states(states_path)

    # Load the eigenvalues
    CY.load_eigenvalues(eigs_path, override_directory_path=False)
    # Load the eigenstates
    CY.load_states(states_path, override_directory_path=False)
    spec = CY.spectrum
    print(spec)
    print(spec[6])
    spec.sort()
    print(spec)
    ArrangeYgraphPlots(CY)
    CY.plot_states(0, plotting_method="surface", realimag="real", N_levels=20)
    #CY.plot_states(1)
    #CY.plot_states(2)
    #CY.plot_states(3)
    #CY.plot_states(4)