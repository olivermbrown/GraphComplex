import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce
import math
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs
import ygraph as yga
import lasso as lsa
import dumbbell as db

def plot_boundary_wavefunction(complex, gluing, show_plots=True, lasso=False, state_number=None,extrapolate_to_origin=False):
    """
    Plot the wavefunction on the boundary of the given complex.

    Args:
        complex (ConfigurationSpace): The complex to plot the wavefunction on.
        gluing (list): The gluing to plot the wavefunction on.
        show_plots (bool): Whether to show the plots.
        lasso (bool): Whether the complex is a Lasso complex.
        state_number (int): The number of the state to plot.

    Returns:
        None.
    """

    # If the state number is not specified, choose the ground state
    if state_number is None:
        state_number = np.argmin(complex.spectrum)
        pass
    else:
        pass
    
    data_real = complex.plot_states(state_number, return_data=True, show_plots=False, realimag="real")
    data_imag = complex.plot_states(state_number, return_data=True, show_plots=False, realimag="imag")

    data = []
    for i in range(len(data_real)):
        x = data_real[i][0]
        y = data_real[i][1]
        z_real = data_real[i][2]
        z_imag = data_imag[i][2]
        #tuple_real = data_real[i]
        #tuple_imag = data_imag[i]
        #real_part = np.array(tuple_real)
        #imag_part = np.array(tuple_imag)
        z = z_real + 1j*z_imag
        vals_tuple = (x, y, z)
        data.append(vals_tuple)
        #data.append(data_real[i] + 1j*data_imag[i])
        pass


    psis = []
    axis_coords = []
    for edge in gluing:
        
        # Import eigenstate data from cell corresponding to the edge
        cell = edge.domain
        i = complex.cells.index(cell)
        cell_data = data[i]

        xs, ys, zs = cell_data

        # Calculate the length of the edge
        if type(cell) == cls.SquareCell:
            edge_length = (cell.N)-2
            pass
        elif type(cell) == cls.TriangleCell:
            if complex.diagonal_boundary_condition == "neumann" or complex.diagonal_boundary_condition == "robin":
                edge_length = (cell.N)-2
                pass
            elif complex.diagonal_boundary_condition == "dirichlet":
                edge_length = (cell.N)-3
                pass
            else:
                raise Exception("Unknown diagonal boundary condition")
            pass
        else:
            raise Exception("Unknown cell type")

        # Store the axis coordinates
        axis_coords.append(ys[:edge_length])

        # Select the function values on the given edge
        if type(cell) == cls.SquareCell:
            if edge == cell.x0:
                psi = zs[:edge_length]
                pass
            elif edge == cell.x1:
                psi = zs[edge_length**2 - edge_length:edge_length**2]
                pass
            elif edge == cell.y0:
                psi = zs[0:edge_length**2:edge_length]
                pass
            elif edge == cell.y1:
                psi = zs[edge_length-1:edge_length**2:edge_length]
                pass
            else:
                raise Exception("Unknown edge")
                pass
            pass
        elif type(cell) == cls.TriangleCell:
            if edge == cell.x0:
                psi = zs[:edge_length]
                # Add on an extra zero to the beginning of the array
                if complex.diagonal_boundary_condition == "neumann" or complex.diagonal_boundary_condition == "robin":
                    pass
                elif complex.diagonal_boundary_condition == "dirichlet":
                    psi = np.concatenate(([0],psi))
                    pass
                else:
                    raise Exception("Unknown diagonal boundary condition")
                pass
            elif edge == cell.y1:
                ids = []
                c = 0
                for i in range(edge_length):
                    ids.append(edge_length*(i+1) - 1 - c - i)
                    c += i
                    pass
                psi = zs[ids]
                # Add on an extra zero to the end of the array
                if complex.diagonal_boundary_condition == "neumann" or complex.diagonal_boundary_condition == "robin":
                    pass
                elif complex.diagonal_boundary_condition == "dirichlet":
                    psi = np.concatenate((psi,[0]))
                    pass
                else:
                    raise Exception("Unknown diagonal boundary condition")
                pass
            else:
                raise Exception("Unknown edge")
            edge_length += 1
            pass
        else:
            raise Exception("Unknown cell type")

        #psis.append(np.abs(psi))
        psis.append(psi)
        pass

    # Calculate the average value of the wavefunction on the edge
    psi = np.mean(psis, axis=0)

    # Convert wavefunction to absolute value
    psi = np.abs(psi)

    # Enable this line for the lasso graph to reverse the direction of the wavefunction on the edge
    if lasso:
        psi = psi[::-1]
        pass
    else:
        pass

    # Define axis coordinates by choosing the vector from axis with the most elements
    axis_coords = axis_coords[np.argmax([len(axis) for axis in axis_coords])]
    ys = axis_coords

    if extrapolate_to_origin:
        # Extrapolate the wavefunction to the origin
        #poly = np.polyfit(ys[:10], psi[:10], deg=5)
        #psi_0  = np.polyval(poly, 0)
        #print(psi_0)
        #psi = psi - psi_0
        psi = psi - psi[0]
        psi = psi[1:]
        ys = ys - ys[0]
        ys = ys[1:]
        pass
    elif not extrapolate_to_origin:
        pass
    else:
        raise Exception("Unknown extrapolation type")


    # Plot the wavefunction on the edge
    plt.clf()
    z = np.abs(psi)- np.sin(np.pi*ys[:edge_length])*np.max(np.abs(psi))

    # Find largest magnitude of the wavefunction
    psi_min = np.min(psi)
    psi_max = np.max(psi)
    if np.abs(psi_min) > np.abs(psi_max):
        max = psi_min
        pass
    elif np.abs(psi_min) < np.abs(psi_max):
        max = psi_max
        pass
    else:
        max = psi_max
        pass
    approx = np.sin(np.pi*ys[:edge_length])*max
    #plt.plot(ys[:edge_length], np.abs(psi))
    plt.plot(ys[:edge_length], psi)
    #plt.show()
    plt.plot(ys[:edge_length],approx,color='red')

    # Plot the wavefunction on the edge
    if show_plots:
        plt.show()
        pass

    # Plot the log log plot of the wavefunction on the edge
    plt.clf()
    plt.scatter(np.log(ys[:edge_length]),np.log(np.abs(psi)))

    # Calculate the slope of the line from first several points
    degree = 4
    cutoff = 10
    polyfeatures = PolynomialFeatures(degree).fit_transform(ys[:cutoff].reshape(-1,1))
    logfeatures = np.log(np.abs(ys[:cutoff]).reshape(-1,1))
    features = np.concatenate((logfeatures,polyfeatures),axis=1)
    print(features.shape)
    reg = LinearRegression(fit_intercept=False).fit(features, np.log(np.abs(psi[:cutoff])))
    m = reg.coef_
    #print(reg.predict(features))
    #print(np.log(np.abs(psi[:10])))
    r = np.linspace(ys[0], ys[cutoff-1], 100)
    plottingfeatures = PolynomialFeatures(degree).fit_transform(r.reshape(-1,1))
    plottinglogfeatures = np.log(np.abs(r).reshape(-1,1))
    plotting = np.concatenate((plottinglogfeatures,plottingfeatures),axis=1)
    #polynomial = np.dot(plotting, m)
    polynomial = reg.predict(plotting)

    # Calculate the gradient of the polynomial through xs[0]
    derivative = np.polyder(m)
    print(m)

    # Plot show the log log plot with the linear fit through the first three points in red
    plt.plot(np.log(r), polynomial, color='red')
    if show_plots:
        plt.show()
        pass
    else:
        plt.close()
        pass

    return m[0]

def check_boundary_convergence(N_min, N_max):

    Ns = np.arange(N_min, N_max+1,10)
    alpha = 1.0

    ms = []

    for N in Ns:
            
        h = (np.pi)/(N-1)

        C = lsa.LassoAnyons(N, alpha)
        states_filepath = "lasso_states/lasso_N" + str(N) + "_alpha" + str(alpha) + "_"
        eigenvalues_filepath = "lasso_eigenvalues/lasso_N" + str(N) + "_alpha" + str(alpha)
        C.load_states(states_filepath)
        C.load_eigenvalues(eigenvalues_filepath)

        gluing1 = C.gluings[0]
        m = plot_boundary_wavefunction(C, gluing1)
        ms.append(m)

        print("N = ", N, ", m = ", m)

        pass

    # Produce scatter plot of m against N

    plt.clf()
    plt.scatter(Ns, ms)
    plt.show()

    return None

def ygraph_hardcore_interpolation(N):

    alphas = np.linspace(0, 1, 25)
    h = (np.pi)/(N-1)

    CYs = []
    gluing1s = []
    ms = []
    for alpha in alphas:
        C = yga.YgraphAnyonsHardcore(N, alpha)
        states_filepath = "ygraph_hardcore_states/ygraph_N" + str(N) + "_alpha" + str(alpha) + "_"
        eigenvalues_filepath = "ygraph_hardcore_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
        C.load_states(states_filepath)
        C.load_eigenvalues(eigenvalues_filepath)
        CYs.append(C)
        gluing = C.gluings[0]
        print("alpha = ", alpha)
        yga.ArrangeYgraphPlots(C)
        m = plot_boundary_wavefunction(C,gluing,show_plots=False)
        ms.append(m)
        pass

    # Produce scatter plot of m against alpha
    plt.clf()
    plt.plot(alphas, ms)
    plt.scatter(alphas, ms)
    # X axis label
    plt.xlabel(r"$\alpha$")
    # Y axis label
    plt.ylabel(r"$\beta$", rotation=0)
    # Rotate y axis label
    plt.yticks(rotation=0)
    plt.show()

    return None

def ygraph_neumann_interpolation(N):

    #alphas = np.linspace(0, 1, 25)
    alphas = np.linspace(0, 1, 49)
    h = (np.pi)/(N-1)

    CYs = []
    gluing1s = []
    ms = []
    for alpha in alphas:
        C = yga.YgraphAnyonsRobin(N, alpha)
        C.robin_constant = 0
        states_filepath = "ygraph_neumann_states/ygraph_N" + str(N) + "_alpha" + str(alpha) + "_"
        eigenvalues_filepath = "ygraph_neumann_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
        C.load_states(states_filepath)
        C.load_eigenvalues(eigenvalues_filepath)
        CYs.append(C)
        gluing = C.gluings[0]
        print("alpha = ", alpha)
        yga.ArrangeYgraphPlots(C)
        m = plot_boundary_wavefunction(C,gluing,show_plots=True,extrapolate_to_origin=True)
        ms.append(m)
        pass

    # Produce scatter plot of m against alpha
    plt.clf()
    plt.plot(alphas, ms)
    plt.scatter(alphas, ms)
    # X axis label
    plt.xlabel(r"$\alpha$")
    # Y axis label
    plt.ylabel(r"$\beta$", rotation=0)
    # Rotate y axis label
    plt.yticks(rotation=0)
    plt.show()

    return None

def ygraph_robin_interpolation(N):

    alphas = np.linspace(0, 1, 25)
    h = (np.pi)/(N-1)

    CYs = []
    gluing1s = []
    ms = []
    for alpha in alphas:
        C = yga.YgraphAnyonsRobin(N, alpha)
        C.robin_constant = 1
        states_filepath = "ygraph_robin_states/ygraph_N" + str(N) + "_alpha" + str(alpha) + "_"
        eigenvalues_filepath = "ygraph_robin_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
        C.load_states(states_filepath)
        C.load_eigenvalues(eigenvalues_filepath)
        CYs.append(C)
        gluing = C.gluings[0]
        print("alpha = ", alpha)
        yga.ArrangeYgraphPlots(C)
        m = plot_boundary_wavefunction(C,gluing,show_plots=False,extrapolate_to_origin=True)
        ms.append(m)
        pass

    # Produce scatter plot of m against alpha
    plt.clf()
    plt.plot(alphas, ms)
    plt.scatter(alphas, ms)
    # X axis label
    plt.xlabel(r"$\alpha$")
    # Y axis label
    plt.ylabel(r"$\beta$", rotation=0)
    # Rotate y axis label
    plt.yticks(rotation=0)
    plt.show()

    return None

def lasso_hardcore_interpolation(N):

    alphas = np.linspace(0, 1, 25)
    h = (np.pi)/(N-1)

    CYs = []
    ms = []
    for alpha in alphas:
        C = lsa.LassoAnyonsHardcore(N, alpha)
        states_filepath = "lasso_hardcore_states/lasso_N" + str(N) + "_alpha" + str(alpha) + "_"
        eigenvalues_filepath = "lasso_hardcore_eigenvalues/lasso_N" + str(N) + "_alpha" + str(alpha)
        C.load_states(states_filepath)
        C.load_eigenvalues(eigenvalues_filepath)
        CYs.append(C)
        #gluing = C.gluings_with_branch_cut[0][0]
        gluing = C.gluings[0]
        print("alpha = ", alpha)
        m = plot_boundary_wavefunction(C,gluing,show_plots=False,lasso=True)
        ms.append(m)
        pass

    # Produce scatter plot of m against alpha
    plt.clf()
    plt.plot(alphas, ms)
    plt.scatter(alphas, ms)
    plt.show()

    return None

def lasso_neumann_interpolation(N):

    #alphas = np.linspace(0, 1, 25)
    alphas = np.linspace(0, 1, 49)
    h = (np.pi)/(N-1)

    CYs = []
    ms = []
    for alpha in alphas:
        C = lsa.LassoAnyonsContact(N, alpha)
        C.robin_constant = 0
        states_filepath = "lasso_neumann_states/lasso_N" + str(N) + "_alpha" + str(alpha) + "_"
        eigenvalues_filepath = "lasso_neumann_eigenvalues/lasso_N" + str(N) + "_alpha" + str(alpha)
        C.load_states(states_filepath)
        C.load_eigenvalues(eigenvalues_filepath)
        CYs.append(C)
        #gluing = C.gluings_with_branch_cut[0][0]
        gluing = C.gluings[0]
        print("alpha = ", alpha)
        m = plot_boundary_wavefunction(C,gluing,show_plots=True,lasso=True,extrapolate_to_origin=True)
        ms.append(m)
        pass

    # Produce scatter plot of m against alpha
    # Produce scatter plot of m against alpha
    plt.clf()
    plt.plot(alphas, ms)
    plt.scatter(alphas, ms)
    # X axis label
    plt.xlabel(r"$\alpha$")
    # Y axis label
    plt.ylabel(r"$\beta$", rotation=0)
    # Rotate y axis label
    plt.yticks(rotation=0)
    plt.show()

    return None

def lasso_robin_interpolation(N):

    #alphas = np.linspace(0, 1, 25)
    alphas = np.linspace(0, 1, 49)
    h = (np.pi)/(N-1)

    CYs = []
    ms = []
    for alpha in alphas:
        C = lsa.LassoAnyonsContact(N, alpha)
        C.robin_constant = 1
        states_filepath = "lasso_robin_states/lasso_N" + str(N) + "_alpha" + str(alpha) + "_"
        eigenvalues_filepath = "lasso_robin_eigenvalues/lasso_N" + str(N) + "_alpha" + str(alpha)
        C.load_states(states_filepath)
        C.load_eigenvalues(eigenvalues_filepath)
        CYs.append(C)
        #gluing = C.gluings_with_branch_cut[0][0]
        gluing = C.gluings[0]
        print("alpha = ", alpha)
        m = plot_boundary_wavefunction(C,gluing,show_plots=False,lasso=True,extrapolate_to_origin=True)
        ms.append(m)
        pass

    # Produce scatter plot of m against alpha
    # Produce scatter plot of m against alpha
    plt.clf()
    plt.plot(alphas, ms)
    plt.scatter(alphas, ms)
    # X axis label
    plt.xlabel(r"$\alpha$")
    # Y axis label
    plt.ylabel(r"$\beta$", rotation=0)
    # Rotate y axis label
    plt.yticks(rotation=0)
    plt.show()

    return None

def dumbbell_hardcore_interpolation(N):
    
        alphas = np.linspace(0, 1, 25)
        h = (np.pi)/(N-1)
    
        CYs = []
        ms = []
        for alpha1 in alphas:
            alpha2 = 0.0
            C = db.DumbbellAnyonsHardcore(N, alpha1, alpha2)
            states_filepath = "dumbbell_hardcore_states/dumbbell_N" + str(N) + "_alpha1" + str(alpha1) + "_alpha2" + str(alpha2) + "_"
            eigenvalues_filepath = "dumbbell_hardcore_eigenvalues/dumbbell_N" + str(N) + "_alpha1" + str(alpha1) + "_alpha2" + str(alpha2)
            C.load_states(states_filepath)
            C.load_eigenvalues(eigenvalues_filepath)
            #C.gen_lapl()
            #C.lapl_solve(h, N_eigs=2)
            CYs.append(C)
            #gluing = C.gluings_with_branch_cut[0][0]
            gluing = C.gluings[0]
            print("alpha = ", alpha1)
            m = plot_boundary_wavefunction(C,gluing,show_plots=False)
            ms.append(m)
            pass
    
        # Produce scatter plot of m against alpha
        plt.clf()
        plt.plot(alphas, ms)
        plt.scatter(alphas, ms)
        # X axis label
        plt.xlabel(r"$\alpha$")
        # Y axis label
        plt.ylabel(r"$\beta$")
        plt.show()
    
        return None

if __name__ =="__main__":

    #check_boundary_convergence(50, 100)

    N = 100
    h = (np.pi)/(N-1)
    #alpha = 1.0

    #C = db.DumbbellAnyons(N, alpha)

    #C = lsa.LassoAnyons(N, alpha)
    # C.gen_lapl()
    # C.lapl_solve(h, N_eigs=50)

    # states_filepath = "lasso_states/lasso_N" + str(N) + "_alpha" + str(alpha) + "_"
    # eigenvalues_filepath = "lasso_eigenvalues/lasso_N" + str(N) + "_alpha" + str(alpha)
    # C.load_states(states_filepath)
    # C.load_eigenvalues(eigenvalues_filepath)
    # print(C.spectrum)

    # states_filepath = "dumbbell_states/dumbbell_N" + str(N) + "_alpha" + str(alpha) + "_"
    # eigenvalues_filepath = "dumbbell_eigenvalues/dumbbell_N" + str(N) + "_alpha" + str(alpha)
    # C.load_states(states_filepath)
    # C.load_eigenvalues(eigenvalues_filepath)
    # print(C.spectrum)

    #states = C.states
    #gs = states[:,0]

    #D12 = C.cells[3]
    #edge = D12.y0

    #gluing1 = C.gluings_with_branch_cut[0][0]

    #plot_boundary_wavefunction(C, gluing1)

    #ygraph_hardcore_interpolation(N)
    #ygraph_neumann_interpolation(N)
    #ygraph_robin_interpolation(N)
    #lasso_hardcore_interpolation(N)
    #lasso_neumann_interpolation(N)
    #lasso_robin_interpolation(N)
    #dumbbell_interpolation(N)
    dumbbell_hardcore_interpolation(N)

    pass