import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs

def StargraphBosonsHardcore(N):

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

    gluing1 = [D11.x0, D12.y0, D13.y0, D14.y0, D15.y0]
    gluing2 = [D22.x0, D12.x0, D23.y0, D24.y0, D25.y0]
    gluing3 = [D33.x0, D13.x0, D23.x0, D34.y0, D35.y0]
    gluing4 = [D44.x0, D14.x0, D24.x0, D34.x0, D45.y0]
    gluing5 = [D55.x0, D15.x0, D25.x0, D35.x0, D45.x0]

    CY.glue(gluing1)
    CY.glue(gluing2)
    CY.glue(gluing3)
    CY.glue(gluing4)
    CY.glue(gluing5)

    CY.exterior_bc("dirichlet")
    CY.diagonal_bc("dirichlet")

    return CY

def ArrangeStargraphPlots(C):

    C.plot_dim = (2,6)

    C.figuresize = (40,10)

    for cell in C.cells:
        if cell.indices == (1,1):
            cell.plot_loc = (1,0)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (1,2):
            cell.plot_loc = (0,0)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (2,2):
            cell.plot_loc = (1,1)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (1,3):
            cell.plot_loc = (0,1)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (2,3):
            cell.plot_loc = (0,3)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (3,3):
            cell.plot_loc = (1,2)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (1,4):
            cell.plot_loc = (0,2)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (2,4):
            cell.plot_loc = (0,4)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (3,4):
            cell.plot_loc = (0,5)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (4,4):
            cell.plot_loc = (1,3)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        else:
            pass
        pass

    C.plots_off = [(0,4), (0,5)]

    return None

if __name__ == "__main__":
    # Main

    N = 100
    h = (np.pi)/(N-1)
    alpha1 = 0.0
    alpha2 = 0.0
    alpha3 = 0.0

    #C = CrossgraphBosons(N)

    C = StargraphBosonsHardcore(N)

    #C.gen_lapl()

    #C.lapl_solve(h,2)

    # Eigenvalues filepath
    eigs_path = "star_hardcore_eigenvalues/star_N"+str(N)+"_alpha"+str(alpha1)
    # Eigenstates filepath
    states_path = "star_hardcore_states/star_N"+str(N)+"_alpha"+str(alpha1)+"_"

    # Load the eigenvalues
    C.load_eigenvalues(eigs_path, override_directory_path=False)
    # Load the eigenstates
    C.load_states(states_path, override_directory_path=False)

    spec = C.spectrum
    #spec.sort()
    print(spec)
    ArrangeStargraphPlots(C)
    C.plot_states(0,plotting_method="surface",realimag="real")