import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs

def DumbbellSubdividedAnyonsHardcore(N, N2, N3, alpha1, alpha2, aharonovbohm1, aharonovbohm2):

    D11 = cls.TriangleCell(N,N)
    D22 = cls.TriangleCell(N,N)
    D33 = cls.TriangleCell(N2,N2)
    D44 = cls.TriangleCell(N3,N3)
    D55 = cls.TriangleCell(N3,N3)
    D12 = cls.SquareCell(N,N)
    D13 = cls.SquareCell(N2,N)
    D14 = cls.SquareCell(N3,N)
    D15 = cls.SquareCell(N3,N)
    D23 = cls.SquareCell(N2,N)
    D24 = cls.SquareCell(N3,N)
    D25 = cls.SquareCell(N3,N)
    D34 = cls.SquareCell(N3,N2)
    D35 = cls.SquareCell(N3,N2)
    D45 = cls.SquareCell(N3,N3)

    D11.indices = (1,1)
    D22.indices = (2,2)
    D33.indices = (3,3)
    D44.indices = (4,4)
    D55.indices = (5,5)
    D12.indices = (1,2)
    D13.indices = (1,3)
    D23.indices = (2,3)
    D14.indices = (1,4)
    D15.indices = (1,5)
    D24.indices = (2,4)
    D25.indices = (2,5)
    D34.indices = (3,4)
    D35.indices = (3,5)
    D45.indices = (4,5)

    C = configs.ConfigurationSpace([D11,D22,D33,D44,D55,D12,D13,D14,D15,D23,D24,D25,D34,D35,D45])

    gluing1 = [D12.y0,D11.y1]
    gluing2 = [D11.x0,D12.y1,D13.y0]
    gluing3 = [D12.x1,D22.x0]
    gluing4 = [D22.y1,D23.y0,D12.x0]
    gluing5 = [D13.x0,D23.x1,D33.x0]
    gluing6 = [D13.x1,D23.x0]
    gluing7 = [D13.y1,D14.y0,D15.y1]
    gluing8 = [D23.y1,D24.y0,D25.y1]
    gluing9 = [D14.x0,D24.x1,D34.x0]
    gluing10 = [D14.x1,D24.x0]
    gluing11 = [D15.y0,D14.y1]
    gluing12 = [D25.y0,D24.y1]
    gluing13 = [D15.x0,D25.x1,D35.x0]
    gluing14 = [D15.x1,D25.x0]
    gluing15 = [D35.y0,D34.y1]
    gluing16 = [D33.y1,D34.y0,D35.y1]
    gluing17 = [D34.x1,D44.x0,D45.y1]
    gluing18 = [D35.x1,D45.x0,D55.y1]
    gluing19 = [D45.y0,D44.y1]
    gluing20 = [D45.x1,D55.x0]

    phase1a = np.exp(1j*np.pi*alpha1)
    phase1b = np.exp(1j*np.pi*aharonovbohm1)
    phase1c = np.exp(1j*np.pi*(aharonovbohm1-alpha1))
    phase2a = np.exp(1j*np.pi*alpha2)
    phase2b = np.exp(1j*np.pi*aharonovbohm2)
    phase2c = np.exp(1j*np.pi*(aharonovbohm2-alpha2))

    C.glue_with_branch_cut(gluing1, phase1c)
    C.glue(gluing2)
    C.glue_with_branch_cut(gluing3, phase1b)
    C.glue_with_branch_cut(gluing4, phase1a)
    C.glue(gluing5)
    C.glue_with_branch_cut(gluing6, phase1b)
    C.glue(gluing7)
    C.glue(gluing8)
    C.glue(gluing9)
    C.glue_with_branch_cut(gluing10, phase1b)
    C.glue_with_branch_cut(gluing11, phase2b)
    C.glue_with_branch_cut(gluing12, phase2b)
    C.glue(gluing13)
    C.glue_with_branch_cut(gluing14, phase1b)
    C.glue_with_branch_cut(gluing15, phase2b)
    C.glue(gluing16)
    C.glue_with_branch_cut(gluing17, phase2a)
    C.glue(gluing18)
    C.glue_with_branch_cut(gluing19, phase2b)
    C.glue_with_branch_cut(gluing20, phase2c)

    C.exterior_bc("dirichlet")
    C.diagonal_bc("dirichlet")

    return C

def DumbbellSubdividedAnyonsNeumann(N, alpha1, alpha2):

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

def ArrangeDumbbellSubdividedPlots(C):

    C.plot_dim = (5,5)

    C.figuresize = (11,10)

    for cell in C.cells:
        if cell.indices == (1,1):
            cell.plot_loc = (4,0)
            cell.use_x_labels = True
            cell.use_y_labels = True
            pass
        elif cell.indices == (1,2):
            cell.plot_loc = (3,0)
            cell.use_x_labels = False
            cell.use_y_labels = True
            pass
        elif cell.indices == (1,3):
            cell.plot_loc = (2,0)
            cell.use_x_labels = False
            cell.use_y_labels = True
            pass
        elif cell.indices == (1,4):
            cell.plot_loc = (1,0)
            cell.use_x_labels = False
            cell.use_y_labels = True
            pass
        elif cell.indices == (1,5):
            cell.plot_loc = (0,0)
            cell.use_x_labels = False
            cell.use_y_labels = True
            pass
        elif cell.indices == (2,2):
            cell.plot_loc = (3,1)
            cell.use_x_labels = True
            cell.use_y_labels = False
            pass
        elif cell.indices == (2,3):
            cell.plot_loc = (2,1)
            cell.use_x_labels = False
            cell.use_y_labels = False
            pass
        elif cell.indices == (2,4):
            cell.plot_loc = (1,1)
            cell.use_x_labels = False
            cell.use_y_labels = False
            pass
        elif cell.indices == (2,5):
            cell.plot_loc = (0,1)
            cell.use_x_labels = False
            cell.use_y_labels = False
            pass
        elif cell.indices == (3,3):
            cell.plot_loc = (2,2)
            cell.use_x_labels = True
            cell.use_y_labels = False
            pass
        elif cell.indices == (3,4):
            cell.plot_loc = (1,2)
            cell.use_x_labels = False
            cell.use_y_labels = False
            pass
        elif cell.indices == (3,5):
            cell.plot_loc = (0,2)
            cell.use_x_labels = False
            cell.use_y_labels = False
            pass
        elif cell.indices == (4,4):
            cell.plot_loc = (1,3)
            cell.use_x_labels = True
            cell.use_y_labels = False
            pass
        elif cell.indices == (4,5):
            cell.plot_loc = (0,3)
            cell.use_x_labels = False
            cell.use_y_labels = False
            pass
        elif cell.indices == (5,5):
            cell.plot_loc = (0,4)
            cell.use_x_labels = True
            cell.use_y_labels = False
            pass
        else:
            pass
        pass

    C.plots_off = [(4,1),(4,2),(4,3),(4,4),
                  (3,2),(3,3),(3,4),
                  (2,3),(2,4),
                  (1,4)]
    
    C.panel_wspace = -0.37
    C.panel_hspace = 0.3

    C.title_size = 24
    C.title_pad = 4.0
    C.label_size = 20
    C.tick_size = 12
    C.legend_size = 20

    C.y_label_coords = (-0.5, 0.57)

    return None

if __name__ == "__main__":
    # Main

    N = 20
    N2 = 100
    h = (np.pi)/(N-1)
    alpha1 = 0.0
    alpha2 = 0.0
    aharonovbohm1 = 0.0
    aharonovbohm2 = 0.0

    C = DumbbellSubdividedAnyonsHardcore(N, N, N, alpha1, alpha2, aharonovbohm1, aharonovbohm2)

    C.gen_lapl()

    C.lapl_solve(h,1)

    # Load the states
    # Eigenvalues filepath
    states_path = "dumbellsubdivided_hardcore_states/dumbellsubdivided_N" + str(N) + "_alpha1=" + str(alpha1) + "_alpha2=" + str(alpha2) + "_aharonovbohm1=" + str(aharonovbohm1) + "_aharonovbohm2=" + str(aharonovbohm2) + "_"
    # Eigenstates filepath
    eigs_path = "dumbellsubdivided_hardcore_eigenvalues/dumbellsubdivided_N" + str(N) + "_alpha1=" + str(alpha1) + "_alpha2=" + str(alpha2) + "_aharonovbohm1=" + str(aharonovbohm1) + "_aharonovbohm2=" + str(aharonovbohm2)

    # Load the eigenvalues
    #C.load_eigenvalues(eigs_path)
    # Load the eigenstates
    #C.load_states(states_path)

    spec = C.spectrum
    spec.sort()
    print(spec)
    ArrangeDumbbellSubdividedPlots(C)
    #C.plot_states(0, plotting_method="contour", realimag="abs", N_levels=20)
    ri = "abs"
    save_name = "dumbbellsubdivided_N"+str(N)+"_alpha1="+str(alpha1)+"_alpha2="+str(alpha2)+"_aharonovbohm1="+str(aharonovbohm1)+"_aharonovbohm2="+str(aharonovbohm2)+"_"+ri+"_contour.pdf"
    C.plot_states(0, plotting_method="contour", realimag=ri, N_levels=20, save_as=save_name)
    # ri = "phase"
    # save_name = "dumbbellsubdivided_N"+str(N)+"_alpha1="+str(alpha1)+"_alpha2="+str(alpha2)+"_aharonovbohm1="+str(aharonovbohm1)+"_aharonovbohm2="+str(aharonovbohm2)+"_"+ri+"_contour.pdf"
    # C.plot_states(0, plotting_method="contour", realimag=ri, N_levels=20, save_as=save_name)
    #C.plot_states(1)
    #C.plot_states(2)

    for p in [0.25, 0.5, 0.75, 1.0]:

        alpha1 = p
        alpha2 = 0.0
        aharonovbohm1 = 0.0
        aharonovbohm2 = 0.0

        C = DumbbellSubdividedAnyonsHardcore(N, N, N, alpha1, alpha2, aharonovbohm1, aharonovbohm2)

        C.gen_lapl()

        C.lapl_solve(h,1)

        ArrangeDumbbellSubdividedPlots(C)

        ri = "abs"
        save_name = "dumbbellsubdivided_N"+str(N)+"_alpha1="+str(alpha1)+"_alpha2="+str(alpha2)+"_aharonovbohm1="+str(aharonovbohm1)+"_aharonovbohm2="+str(aharonovbohm2)+"_"+ri+"_contour.pdf"
        C.plot_states(0, plotting_method="contour", realimag=ri, N_levels=20, save_as=save_name)

        ri = "phase"
        save_name = "dumbbellsubdivided_N"+str(N)+"_alpha1="+str(alpha1)+"_alpha2="+str(alpha2)+"_aharonovbohm1="+str(aharonovbohm1)+"_aharonovbohm2="+str(aharonovbohm2)+"_"+ri+"_contour.pdf"
        C.plot_states(0, plotting_method="contour", realimag=ri, N_levels=20, save_as=save_name)

        pass

    for p in [0.25, 0.5, 0.75, 1.0]:

        alpha1 = 0.0
        alpha2 = p
        aharonovbohm1 = 0.0
        aharonovbohm2 = 0.0

        C = DumbbellSubdividedAnyonsHardcore(N, N, N, alpha1, alpha2, aharonovbohm1, aharonovbohm2)

        C.gen_lapl()

        C.lapl_solve(h,1)

        ArrangeDumbbellSubdividedPlots(C)

        ri = "abs"
        save_name = "dumbbellsubdivided_N"+str(N)+"_alpha1="+str(alpha1)+"_alpha2="+str(alpha2)+"_aharonovbohm1="+str(aharonovbohm1)+"_aharonovbohm2="+str(aharonovbohm2)+"_"+ri+"_contour.pdf"
        C.plot_states(0, plotting_method="contour", realimag=ri, N_levels=20, save_as=save_name)

        ri = "phase"
        save_name = "dumbbellsubdivided_N"+str(N)+"_alpha1="+str(alpha1)+"_alpha2="+str(alpha2)+"_aharonovbohm1="+str(aharonovbohm1)+"_aharonovbohm2="+str(aharonovbohm2)+"_"+ri+"_contour.pdf"
        C.plot_states(0, plotting_method="contour", realimag=ri, N_levels=20, save_as=save_name)

        pass

    for p in [0.25, 0.5, 0.75, 1.0]:

        alpha1 = 0.0
        alpha2 = 0.0
        aharonovbohm1 = p
        aharonovbohm2 = 0.0

        C = DumbbellSubdividedAnyonsHardcore(N, N, N, alpha1, alpha2, aharonovbohm1, aharonovbohm2)

        C.gen_lapl()

        C.lapl_solve(h,1)

        ArrangeDumbbellSubdividedPlots(C)

        ri = "abs"
        save_name = "dumbbellsubdivided_N"+str(N)+"_alpha1="+str(alpha1)+"_alpha2="+str(alpha2)+"_aharonovbohm1="+str(aharonovbohm1)+"_aharonovbohm2="+str(aharonovbohm2)+"_"+ri+"_contour.pdf"
        C.plot_states(0, plotting_method="contour", realimag=ri, N_levels=20, save_as=save_name)

        ri = "phase"
        save_name = "dumbbellsubdivided_N"+str(N)+"_alpha1="+str(alpha1)+"_alpha2="+str(alpha2)+"_aharonovbohm1="+str(aharonovbohm1)+"_aharonovbohm2="+str(aharonovbohm2)+"_"+ri+"_contour.pdf"
        C.plot_states(0, plotting_method="contour", realimag=ri, N_levels=20, save_as=save_name)

        pass

    for p in [0.25, 0.5, 0.75, 1.0]:

        alpha1 = 0.0
        alpha2 = 0.0
        aharonovbohm1 = 0.0
        aharonovbohm2 = p

        C = DumbbellSubdividedAnyonsHardcore(N, N, N, alpha1, alpha2, aharonovbohm1, aharonovbohm2)

        C.gen_lapl()

        C.lapl_solve(h,1)

        ArrangeDumbbellSubdividedPlots(C)

        ri = "abs"
        save_name = "dumbbellsubdivided_N"+str(N)+"_alpha1="+str(alpha1)+"_alpha2="+str(alpha2)+"_aharonovbohm1="+str(aharonovbohm1)+"_aharonovbohm2="+str(aharonovbohm2)+"_"+ri+"_contour.pdf"
        C.plot_states(0, plotting_method="contour", realimag=ri, N_levels=20, save_as=save_name)

        ri = "phase"
        save_name = "dumbbellsubdivided_N"+str(N)+"_alpha1="+str(alpha1)+"_alpha2="+str(alpha2)+"_aharonovbohm1="+str(aharonovbohm1)+"_aharonovbohm2="+str(aharonovbohm2)+"_"+ri+"_contour.pdf"
        C.plot_states(0, plotting_method="contour", realimag=ri, N_levels=20, save_as=save_name)

        pass