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

def plot_boundary_wavefunction(complex, gluing, show_plots=True, lasso=False, state_number=None):
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
    data = complex.plot_states(state_number, return_data=True, show_plots=False)

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
            edge_length = (cell.N)-3
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
                psi = np.concatenate(([0],psi))
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
                psi = np.concatenate((psi,[0]))
                pass
            else:
                raise Exception("Unknown edge")
                pass
            edge_length += 1
            pass
        else:
            raise Exception("Unknown cell type")
            pass

        psis.append(np.abs(psi))
        pass

    # Calculate the average value of the wavefunction on the edge
    psi = np.mean(psis, axis=0)

    # Enable this line for the lasso graph to reverse the direction of the wavefunction on the edge
    if lasso:
        psi = psi[::-1]
        pass
    else:
        pass

    # Define axis coordinates by choosing the vector from axis with the most elements
    axis_coords = axis_coords[np.argmax([len(axis) for axis in axis_coords])]
    ys = axis_coords

    # Plot the wavefunction on the edge
    plt.clf()
    z = np.abs(psi)- np.sin(np.pi*ys[:edge_length])*np.max(np.abs(psi))
    approx = np.sin(np.pi*ys[:edge_length])*np.max(np.abs(psi))
    plt.plot(ys[:edge_length], np.abs(psi))
    #plt.show()
    plt.plot(ys[:edge_length], approx,color='red')

    # Plot the wavefunction on the edge
    if show_plots:
        plt.show()
        pass

    # Plot the log log plot of the wavefunction on the edge
    plt.clf()
    plt.scatter(np.log(ys[:edge_length]),np.log(np.abs(psi)))

    # Calculate the slope of the line from first point points
    degree = 2
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

def ygraph_interpolation(N):

    alphas = np.linspace(0, 1, 25)
    h = (np.pi)/(N-1)

    CYs = []
    gluing1s = []
    ms = []
    for alpha in alphas:
        C = yga.YgraphAnyons(N, alpha)
        states_filepath = "ygraph_states/ygraph_N" + str(N) + "_alpha" + str(alpha) + "_"
        eigenvalues_filepath = "ygraph_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
        C.load_states(states_filepath)
        C.load_eigenvalues(eigenvalues_filepath)
        CYs.append(C)
        gluing = C.gluings[0]
        print("alpha = ", alpha)
        m = plot_boundary_wavefunction(C,gluing,show_plots=False)
        ms.append(m)
        pass

    # Produce scatter plot of m against alpha
    plt.clf()
    plt.scatter(alphas, ms)
    plt.show()

    return None

def bosons_ygraph_interpolation(N):

    #alphas = np.linspace(0, 1, 25)
    alphas = [0]
    h = (np.pi)/(N-1)

    CYs = []
    gluing1s = []
    ms = []
    for alpha in alphas:
        C = yga.YgraphBosonsFree(N, alpha)
        states_filepath = "bosons_ygraph_states/ygraph_N" + str(N) + "_alpha" + str(alpha) + "_"
        eigenvalues_filepath = "bosons_ygraph_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
        C.load_states(states_filepath)
        C.load_eigenvalues(eigenvalues_filepath)
        CYs.append(C)
        gluing = C.gluings[0]
        print("alpha = ", alpha)
        m = plot_boundary_wavefunction(C,gluing,show_plots=True,state_number=2)
        ms.append(m)
        pass

    # Produce scatter plot of m against alpha
    plt.clf()
    plt.scatter(alphas, ms)
    plt.show()

    return None

def lasso_interpolation(N):

    alphas = np.linspace(0, 1, 25)
    h = (np.pi)/(N-1)

    CYs = []
    ms = []
    for alpha in alphas:
        C = lsa.LassoAnyons(N, alpha)
        states_filepath = "lasso_states/lasso_N" + str(N) + "_alpha" + str(alpha) + "_"
        eigenvalues_filepath = "lasso_eigenvalues/lasso_N" + str(N) + "_alpha" + str(alpha)
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
    plt.scatter(alphas, ms)
    plt.show()

    return None

def dumbbell_interpolation(N):
    
        alphas = np.linspace(0, 1, 25)
        h = (np.pi)/(N-1)
    
        CYs = []
        ms = []
        for alpha in alphas:
            C = db.DumbbellAnyons(N, alpha)
            states_filepath = "dumbbell_states/dumbbell_N" + str(N) + "_alpha" + str(alpha) + "_"
            eigenvalues_filepath = "dumbbell_eigenvalues/dumbbell_N" + str(N) + "_alpha" + str(alpha)
            C.load_states(states_filepath)
            C.load_eigenvalues(eigenvalues_filepath)
            CYs.append(C)
            #gluing = C.gluings_with_branch_cut[0][0]
            gluing = C.gluings[0]
            print("alpha = ", alpha)
            m = plot_boundary_wavefunction(C,gluing,show_plots=False)
            ms.append(m)
            pass
    
        # Produce scatter plot of m against alpha
        plt.clf()
        plt.scatter(alphas, ms)
        plt.show()
    
        return None

if __name__ =="__main__":

    #check_boundary_convergence(50, 100)

    N = 70
    h = (np.pi)/(N-1)
    alpha = 1.0

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

    #ygraph_interpolation(N)
    #bosons_ygraph_interpolation(N)
    #lasso_interpolation(N)
    dumbbell_interpolation(N)

    pass