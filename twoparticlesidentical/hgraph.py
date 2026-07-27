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

def HgraphAnyonsHardcore(N, alpha1, alpha2):

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

    CH = configs.ConfigurationSpace([D11,D22,D33,D44,D55,D12,D13,D14,D15,D23,D24,D25,D34,D35,D45])

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

    phase1 = np.exp(1j*np.pi*alpha1)
    phase2 = np.exp(1j*np.pi*alpha2)

    CH.glue(gluing1)
    CH.glue(gluing2)
    CH.glue_with_branch_cut(gluing3, phase1)
    CH.glue_with_branch_cut(gluing4, phase2)
    CH.glue(gluing5)
    CH.glue(gluing6)
    CH.glue(gluing7)
    CH.glue(gluing8)
    CH.glue(gluing9)
    CH.glue(gluing10)

    CH.exterior_bc("dirichlet")
    CH.diagonal_bc("dirichlet")

    return CH

def HgraphAnyonsDelta(N, alpha1, alpha2, robin_constant):

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

    CH = configs.ConfigurationSpace([D11,D22,D33,D44,D55,D12,D13,D14,D15,D23,D24,D25,D34,D35,D45])

    CH.robin_constant = robin_constant

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

    phase1 = np.exp(1j*np.pi*alpha1)
    phase2 = np.exp(1j*np.pi*alpha2)

    CH.glue(gluing1)
    CH.glue(gluing2)
    CH.glue_with_branch_cut(gluing3, phase1)
    CH.glue_with_branch_cut(gluing4, phase2)
    CH.glue(gluing5)
    CH.glue(gluing6)
    CH.glue(gluing7)
    CH.glue(gluing8)
    CH.glue(gluing9)
    CH.glue(gluing10)

    CH.diagonal_bc("robin")
    CH.exterior_bc("dirichlet")

    return CH

def ArrangeHgraphPlots(C):

    C.plot_dim = (3,5)

    C.figuresize = (18,10)

    for cell in C.cells:
        if cell.indices == (1,1):
            cell.plot_loc = (2,0)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (1,2):
            cell.plot_loc = (0,0)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (1,3):
            cell.plot_loc = (0,1)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (1,4):
            cell.plot_loc = (0,2)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (1,5):
            cell.plot_loc = (0,3)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (2,2):
            cell.plot_loc = (2,1)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (2,3):
            cell.plot_loc = (0,4)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (2,4):
            cell.plot_loc = (1,0)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (2,5):
            cell.plot_loc = (1,1)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (3,3):
            cell.plot_loc = (2,2)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (3,4):
            cell.plot_loc = (1,2)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (3,5):
            cell.plot_loc = (1,3)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (4,4):
            cell.plot_loc = (2,3)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (4,5):
            cell.plot_loc = (1,4)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (5,5):
            cell.plot_loc = (2,4)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        else:
            pass
        pass
    
    C.panel_wspace = -0.1
    C.panel_hspace = 0.6
    C.title_size = 30
    C.title_pad = 10.0
    C.label_size = 24
    C.tick_size = 16
    C.legend_size = 24

    C.y_label_coords = (-0.35, 0.58)

    return None

if __name__ == "__main__":
    # Main

    N = 20
    h = (np.pi)/(N-1)
    alpha1 = 1.0
    alpha2 = 1.0

    CH = HgraphAnyonsHardcore(N,alpha1,alpha2)

    CH.gen_lapl()

    CH.lapl_solve(h,10)

    spec = CH.spectrum
    spec.sort()
    print(spec)
    ArrangeHgraphPlots(CH)
    ri = "abs"
    save_name = "Hgraph_N"+str(N)+"_alpha1="+str(alpha1)+"_alpha2="+str(alpha2)+"_"+ri+"_contour.pdf"
    CH.plot_states(0,plotting_method="contour",realimag=ri,N_levels=20, save_as=save_name)

    # for p in [0.25, 0.5, 0.75, 1.0]:

    #     alpha1 = p
    #     alpha2 = 0.0

    #     CH = HgraphAnyonsHardcore(N,alpha1,alpha2)

    #     CH.gen_lapl()

    #     CH.lapl_solve(h,1)

    #     ArrangeHgraphPlots(CH)

    #     ri = "abs"
    #     save_name = "Hgraph_N"+str(N)+"_alpha1="+str(alpha1)+"_alpha2="+str(alpha2)+"_"+ri+"_contour.pdf"
    #     CH.plot_states(0,plotting_method="contour",realimag=ri,N_levels=20, save_as=save_name)

    #     ri = "phase"
    #     save_name = "Hgraph_N"+str(N)+"_alpha1="+str(alpha1)+"_alpha2="+str(alpha2)+"_"+ri+"_contour.pdf"
    #     CH.plot_states(0,plotting_method="contour",realimag=ri,N_levels=20, save_as=save_name)

    #     pass

    # for p in [0.25, 0.5, 0.75, 1.0]:

    #     alpha1 = 0.0
    #     alpha2 = p

    #     CH = HgraphAnyonsHardcore(N,alpha1,alpha2)

    #     CH.gen_lapl()

    #     CH.lapl_solve(h,1)

    #     ArrangeHgraphPlots(CH)

    #     ri = "abs"
    #     save_name = "Hgraph_N"+str(N)+"_alpha1="+str(alpha1)+"_alpha2="+str(alpha2)+"_"+ri+"_contour.pdf"
    #     CH.plot_states(0,plotting_method="contour",realimag=ri,N_levels=20, save_as=save_name)

    #     ri = "phase"
    #     save_name = "Hgraph_N"+str(N)+"_alpha1="+str(alpha1)+"_alpha2="+str(alpha2)+"_"+ri+"_contour.pdf"
    #     CH.plot_states(0,plotting_method="contour",realimag=ri,N_levels=20, save_as=save_name)

    #     pass