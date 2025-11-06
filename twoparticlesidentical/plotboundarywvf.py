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
import squarecell as sc
import twowires as tw
import lassosubdivided as lss

def plot_boundary_wavefunction(complex, gluing, show_plots=True, lasso=False, state_number=None, extrapolate_to_origin=False):
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
    #plt.clf()
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
    if show_plots:
        plt.close('all')
        plt.plot(ys[:edge_length], psi)
        plt.plot(ys[:edge_length],approx,color='red')
        plt.show()
        pass
    else:
        plt.close('all')
        pass

    # Plot the log log plot of the wavefunction on the edge
    plt.clf()
    if show_plots:
        plt.scatter(np.log(ys[:edge_length]),np.log(np.abs(psi)))
        pass
    else:
        plt.close('all')

    # Calculate the slope of the line from first several points
    degree = 1
    start = 3
    cutoff = 35
    polyfeatures = PolynomialFeatures(degree).fit_transform(ys[start:cutoff].reshape(-1,1))
    logfeatures = np.log(np.abs(ys[start:cutoff]).reshape(-1,1))
    features = np.concatenate((logfeatures,polyfeatures),axis=1)
    print(features.shape)
    reg = LinearRegression(fit_intercept=False).fit(features, np.log(np.abs(psi[start:cutoff])))
    m = reg.coef_
    #print(reg.predict(features))
    #print(np.log(np.abs(psi[:10])))
    r = np.linspace(ys[start], ys[cutoff-1], 100)
    plottingfeatures = PolynomialFeatures(degree).fit_transform(r.reshape(-1,1))
    plottinglogfeatures = np.log(np.abs(r).reshape(-1,1))
    plotting = np.concatenate((plottinglogfeatures,plottingfeatures),axis=1)
    #polynomial = np.dot(plotting, m)
    polynomial = reg.predict(plotting)

    # Calculate the gradient of the polynomial through xs[0]
    derivative = np.polyder(m)
    print(m)

    # 1. Compute predictions and residuals
    y_true = np.log(np.abs(psi[start:cutoff]))
    y_pred = reg.predict(features)
    residuals = y_true - y_pred

    # 2. Degrees of freedom
    n_samples, n_features = features.shape
    dof = n_samples - n_features

    # 3. Estimate residual variance
    residual_variance = np.sum(residuals**2) / dof

    # 4. Covariance matrix of coefficients: (XᵀX)^(-1) * σ²
    XTX_inv = np.linalg.inv(features.T @ features)
    cov_matrix = residual_variance * XTX_inv

    # 5. Standard error of the first coefficient (m[0])
    se_m0 = np.sqrt(cov_matrix[0, 0])
    print("Standard error of the first coefficient (m[0]):", se_m0)

    # Plot show the log log plot with the linear fit through the first three points in red
    if show_plots:
        plt.plot(np.log(r), polynomial, color='red')
        plt.show()
        pass
    else:
        plt.close('all')
        pass

    return m[0], se_m0

def plot_hypotenuse_wavefunction(complex, cell, show_plots=True, lasso=False, state_number=None, extrapolate_to_origin=False, start=2, cutoff=25):
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
    
        
    # Import eigenstate data from cell
    i = complex.cells.index(cell)
    cell_data = data[i]

    xs, ys, zs = cell_data

    if type(cell) == cls.TriangleCell:
        raise Exception("This function is only for wavefunctions along the hypotenuse of square cells")
    elif type(cell) == cls.SquareCell:
        # If the cell is a square cell, we can proceed

        N = cell.N - 2
        # Find the indices of the hyptonuse (diagonal) of the square cell
        if N % 2 == 0:
            # If the cell has an even number of nodes, the hypotenuse is the diagonal from (0,0) to (N-1,N-1)
            hypotenuse_indices = np.arange(0, N**2, N+1)

            psi = zs[hypotenuse_indices]

            axis_coords.append(ys[hypotenuse_indices])
            
            pass
        elif N % 2 == 1:
            raise Exception("This function is only for square cells with an even number of nodes")
        else:
            raise Exception("Unknown cell type")
        pass
    else:
        raise Exception("Unknown cell type")

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
    ys = axis_coords*np.sqrt(2)

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
    #plt.clf()
    z = np.abs(psi)- np.sin(np.pi*ys)*np.max(np.abs(psi))

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
    approx = np.sin(np.pi*ys)*max
    #plt.plot(ys[:edge_length], np.abs(psi))
    if show_plots:
        plt.close('all')
        plt.plot(ys, psi)
        sigma = 0.75
        sigmaanalytic = 0.783653
        sigmablowup = 0.9
        psirminussigma = psi*(ys**(-sigma))
        psirminussigmaanalytic = psi*(ys**(-sigmaanalytic))
        psirminussigmablowup = psi*(ys**(-sigmablowup))
        plt.plot(ys, psirminussigma, color='green')
        plt.plot(ys, psirminussigmaanalytic, color='orange')
        plt.plot(ys, psirminussigmablowup, color='purple')
        #plt.plot(ys,approx,color='red')
        #plt.show()
        pass
    else:
        plt.close('all')
        pass

    # Plot the log log plot of the wavefunction on the edge
    plt.clf()
    if show_plots:
        plt.scatter(np.log(ys),np.log(np.abs(psi)), marker='x')
        pass
    else:
        plt.close('all')

    DlogPsi = np.log(np.abs(psi[1])) - np.log(np.abs(psi[0]))
    DlogR = np.log(np.abs(ys[1])) - np.log(np.abs(ys[0]))
    grad = DlogPsi / DlogR
    print("Gradient of the wavefunction at the origin: ", grad)

    # Calculate the slope of the line from first several points
    degree = 1
    #start = 2
    #cutoff = 25
    polyfeatures = PolynomialFeatures(degree).fit_transform(ys[start:cutoff].reshape(-1,1))
    logfeatures = np.log(np.abs(ys[start:cutoff]).reshape(-1,1))
    features = np.concatenate((logfeatures,polyfeatures),axis=1)
    print(features.shape)
    reg = LinearRegression(fit_intercept=False).fit(features, np.log(np.abs(psi[start:cutoff])))
    m = reg.coef_
    #print(reg.predict(features))
    #print(np.log(np.abs(psi[:10])))
    r = np.linspace(ys[start], ys[cutoff-1], 100)
    plottingfeatures = PolynomialFeatures(degree).fit_transform(r.reshape(-1,1))
    plottinglogfeatures = np.log(np.abs(r).reshape(-1,1))
    plotting = np.concatenate((plottinglogfeatures,plottingfeatures),axis=1)
    #polynomial = np.dot(plotting, m)
    polynomial = reg.predict(plotting)

    # Calculate the gradient of the polynomial through xs[0]
    derivative = np.polyder(m)
    print(m)

    # Plot show the log log plot with the linear fit through the first three points in red
    if show_plots:
        plt.plot(np.log(r), polynomial, color='darkorange')
        plt.xlabel(r"$\log(R)$", fontsize=14)
        plt.ylabel(r"$\log(\psi)$", fontsize=14)
        plt.legend(['Data', 'Linear fit'], fontsize=14, loc='lower left')
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.show()
        pass
    else:
        plt.close('all')
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

def ygraph_hardcore_analytic():

    betas = [0.783653, 0.784082, 0.785365, 0.7875, 0.790479, 0.794293, 0.798929, 0.804372, 0.810607, 0.817615, 0.825376, 0.833868, 0.843069, 0.852956, 0.863504, 0.874689, 0.886487, 0.898872, 0.911819, 0.925306, 0.939307, 0.9538, 0.968761, 0.984168, 1]

    return betas

def ygraph_hardcore_interpolation(N):

    alphas = np.linspace(0, 1, 25)
    h = (np.pi)/(N-1)

    CYs = []
    gluing1s = []
    ms = []
    ses = []
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
        m, se = plot_boundary_wavefunction(C,gluing,show_plots=False)
        ms.append(m)
        ses.append(se)
        pass

    # Produce scatter plot of m against alpha
    plt.clf()
    #plt.plot(alphas, ms)
    betas_analytic = ygraph_hardcore_analytic()
    plt.errorbar(alphas, ms, yerr=ses, fmt='.', capsize=5, label='Data with error bars')
    plt.plot(alphas, betas_analytic, label='Analytic solution')
    #plt.scatter(alphas, ms)
    # X axis label
    plt.xlabel(r"$\alpha$")
    # Y axis label
    plt.ylabel(r"$\beta$", rotation=0)
    # Rotate y axis label
    plt.yticks(rotation=0)
    plt.show()

    return None

def ygraph_hardcore_interpolation_hypotenuse(N):

    alphas = np.linspace(0, 1, 25)
    h = (np.pi)/(N-1)

    CYs = []
    ms = []
    us = []
    ls= []
    for alpha in alphas:
        C = yga.YgraphAnyonsHardcore(N, alpha)
        states_filepath = "ygraph_hardcore_states/ygraph_N" + str(N) + "_alpha" + str(alpha) + "_"
        eigenvalues_filepath = "ygraph_hardcore_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
        C.load_states(states_filepath)
        C.load_eigenvalues(eigenvalues_filepath)
        #C.gen_lapl()
        #C.lapl_solve(h, N_eigs=2)
        CYs.append(C)
        cell = C.cells[3]
        print("alpha = ", alpha)
        yga.ArrangeYgraphPlots(C)
        m = plot_hypotenuse_wavefunction(C,cell,show_plots=False, cutoff=20)
        u = plot_hypotenuse_wavefunction(C,cell,show_plots=False, cutoff=30)
        l = plot_hypotenuse_wavefunction(C,cell,show_plots=False, cutoff=10)
        ms.append(m)
        us.append(u)
        ls.append(l)
        pass

    # Set errors to the max of each of us or ls
    us = np.array(us)
    ls = np.array(ls)
    ds = np.maximum(us, ls)
    ers = ds - ms
    ers = np.abs(ers)

    # Produce scatter plot of m against alpha
    plt.clf()
    betas_analytic = ygraph_hardcore_analytic()
    plt.errorbar(alphas, ms, yerr=ers, fmt='.', capsize=5, label='Numerical estimate with error bars')
    plt.plot(alphas, betas_analytic, label='Analytic solution')
    #plt.plot(alphas, ms)
    #plt.scatter(alphas, ms, marker = '+')
    # X axis label
    plt.xlabel(r"$\alpha$", fontsize=14)
    # Y axis label
    plt.ylabel(r"$\sigma$", rotation=0, fontsize=14)
    # Rotate y axis label
    plt.yticks(rotation=0)
    # Fontsize settings
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=14)
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

def lassosubdivided_hardcore_interpolation(N):

    alphas = np.linspace(0, 1, 25)
    h = (np.pi)/(N-1)

    ab = 0.0 # Aharonov-Bohm phase

    CYs = []
    ms = []
    for alpha in alphas:
        C = lss.LassoSubdividedAnyonsHardcore(N, N, alpha, ab)
        states_filepath = "lassosubdivided_hardcore_states/lassosubdivided_N" + str(N) + "_alpha" + str(alpha) + "_aharonovbohm" + str(ab) + "_"
        eigenvalues_filepath = "lassosubdivided_hardcore_eigenvalues/lassosubdivided_N" + str(N) + "_alpha" + str(alpha) + "_aharonovbohm" + str(ab)
        C.load_states(states_filepath)
        C.load_eigenvalues(eigenvalues_filepath)
        CYs.append(C)
        #gluing = C.gluings_with_branch_cut[0][0]
        gluing = C.gluings[0]
        print("alpha = ", alpha)
        m, _ = plot_boundary_wavefunction(C,gluing,show_plots=False,lasso=True)
        ms.append(m)
        pass

    # Produce scatter plot of m against alpha
    plt.clf()
    #plt.plot(alphas, ms)
    #betas_analytic = ygraph_hardcore_analytic()
    #plt.errorbar(alphas, ms, yerr=ses, fmt='.', capsize=5, label='Data with error bars')
    plt.plot(alphas, ms, label='Analytic solution')
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

def singlesquarecell_hardcore_interpolation(N):

    alphas = np.linspace(0, 1, 25)
    h = (np.pi)/(N-1)

    Cs = []
    gluing1s = []
    ms = []
    for alpha in alphas:
        C = sc.Square(N, alpha)
        #states_filepath = "ygraph_hardcore_states/ygraph_N" + str(N) + "_alpha" + str(alpha) + "_"
        #eigenvalues_filepath = "ygraph_hardcore_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
        #C.load_states(states_filepath)
        #C.load_eigenvalues(eigenvalues_filepath)
        C.gen_lapl()
        C.lapl_solve(h, N_eigs=2)
        Cs.append(C)
        print("alpha = ", alpha)
        C.plot_states(0, plotting_method="surface", realimag="abs", N_levels=20)
        yga.ArrangeYgraphPlots(C)
        m = plot_hypotenuse_wavefunction(C, C.cells[0], show_plots=True)
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

def twowires_hardcore_interpolation(N):

    alphas = np.linspace(0, 1, 25)
    h = (np.pi)/(N-1)

    ab = 0.0 # Aharonov-Bohm phase

    CYs = []
    ms = []
    for alpha in alphas:
        C = tw.TwoWiresWithPhases(N, alpha, ab)
        states_filepath = "twowires_hardcore_states/twowires_N" + str(N) + "_alpha1" + str(ab) + "_alpha2" + str(alpha) + "_"
        eigenvalues_filepath = "twowires_hardcore_eigenvalues/twowires_N" + str(N) + "_alpha1" + str(ab) + "_alpha2" + str(alpha)
        C.load_states(states_filepath)
        C.load_eigenvalues(eigenvalues_filepath)
        CYs.append(C)
        gluing = C.gluings_with_branch_cut[0][0]
        #gluing = C.gluings[0]
        print("alpha = ", alpha)
        m, _ = plot_boundary_wavefunction(C,gluing,show_plots=False,lasso=False)
        ms.append(m)
        pass

    # Produce scatter plot of m against alpha
    plt.clf()
    #plt.plot(alphas, ms)
    #betas_analytic = ygraph_hardcore_analytic()
    #plt.errorbar(alphas, ms, yerr=ses, fmt='.', capsize=5, label='Data with error bars')
    plt.plot(alphas, ms, label='Analytic solution')
    plt.scatter(alphas, ms)
    # X axis label
    plt.xlabel(r"$\alpha$")
    # Y axis label
    plt.ylabel(r"$\beta$", rotation=0)
    # Rotate y axis label
    plt.yticks(rotation=0)
    plt.show()

    return None

def model_config_space(N):

    D = cls.SquareCell(N)
    D.indices = (1,1)

    C = configs.ConfigurationSpace([D])

    for e in D.edges:
        e.eliminated = True

    return C

def model_function(N, alpha):
    """
    Plot the function f(r,theta) = r^sigma ( e^(i pi alpha) + e^(-i pi alpha) ) on a square grid
    """

    # Parameters
    sigma = 0.5     # Power of r

    # Create a square Cartesian grid
    x = np.linspace(0, 1, N)
    y = np.linspace(0, 1, N)
    x = x[1:-1]
    y = y[1:-1]
    X, Y = np.meshgrid(x, y)

    # Convert to polar coordinates
    R = np.sqrt(X**2 + Y**2)
    Theta = np.arctan2(Y, X)  # Not used in this case, but useful if needed

    
    # Compute the function
    f = (R ** sigma) * (1 + R) * (2 * np.cos(np.pi * alpha * Theta))

    # Flatten to a column vector of shape (N^2, 1)
    flat = f.reshape((N-2)**2, 1)  # or a.flatten().reshape(-1, 1)

    return flat

def model_function_interpolation(N):

    alphas = np.linspace(0, 1, 25)
    h = (np.pi)/(N-1)

    Cs = []
    ms = []
    for alpha in alphas:
        C = model_config_space(N)
        cell = C.cells[0]
        cell.eigenstates = model_function(N, alpha)
        C.spectrum = [1]
        Cs.append(C)
        print("alpha = ", alpha)
        m = plot_hypotenuse_wavefunction(C,cell,show_plots=True)
        ms.append(m)
        pass

    # Produce scatter plot of m against alpha
    plt.clf()
    #plt.plot(alphas, ms)
    #plt.errorbar(alphas, ms, fmt='.', capsize=5, label='Data with error bars')
    plt.plot(alphas, ms, label='Model function')
    #plt.scatter(alphas, ms)
    # X axis label
    plt.xlabel(r"$\alpha$")
    # Y axis label
    plt.ylabel(r"$\beta$", rotation=0)
    # Rotate y axis label
    plt.yticks(rotation=0)
    plt.show()

    return None

if __name__ =="__main__":

    #check_boundary_convergence(50, 100)

    N = 150
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
    #dumbbell_hardcore_interpolation(N)
    ygraph_hardcore_interpolation_hypotenuse(N)
    #singlesquarecell_hardcore_interpolation(N)
    #lassosubdivided_hardcore_interpolation(N)
    #twowires_hardcore_interpolation(N)
    #model_function_interpolation(N)

    pass