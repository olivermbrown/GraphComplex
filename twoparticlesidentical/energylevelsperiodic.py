# This code is used to check if the energy levels of the Hamiltonian are periodic in alpha.

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

def plot_energy_levels(Cs, alphas, show_plots=True):
    """
    Create a scatter plot showing the energy levels on the graph for each value of alpha.

    Args:
        Cs (list): A list of ConfigurationSpace objects.
        show_plots (bool): Whether to display the plots.
    
    Returns:
        None.
    """

    fig, ax = plt.subplots()

    ys = []
    for C in Cs:
        spec = C.spectrum
        spec = np.sort(spec)
        ys.append(spec)
        pass

    minimum_eigenvalue = 0
    maximum_eigenvalue = 2

    ys = np.array(ys)[:,minimum_eigenvalue:maximum_eigenvalue]

    fig.set_size_inches(6, 10)

    for i in range(ys.shape[1]):
        #ax.plot(alphas, ys[:,i])
        ax.plot(alphas, ys[:,i], label="Energy level " + str(i))

    ax.set_xlabel(r"$\alpha$")
    ax.set_ylabel("Eigenenergy")
    # Locate legend outside of plot
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    if show_plots:
        plt.show()
        pass



    return None

def check_for_stationary_point(Cs, alphas, x0):

    if x0 == 0.0:
        pass
    elif x0 == 1.0:
        # Reverse the order of the alphas
        alphas = alphas[::-1]
        # Reverse the order of the Cs
        Cs = Cs[::-1]
        pass
    else:
        raise ValueError("x0 must be either 0.0 or 1.0.")

    # Create a list of the energy levels for each value of alpha
    ys = []
    for C in Cs:
        spec = C.spectrum
        spec = np.sort(spec)
        ys.append(spec)
        pass

    minimum_eigenvalue = 0
    maximum_eigenvalue = 30

    ys = np.array(ys)[:,minimum_eigenvalue:maximum_eigenvalue]

    derivs = []
    for alpha in alphas[1:]:
        idx = np.where(alphas == alpha)[0][0]
        diff = ys[idx] - ys[idx-1]
        delta_alpha = alpha - alphas[0]
        deriv = diff/delta_alpha
        derivs.append(deriv)
        pass

    derivs = np.array(derivs)

    # Plot the derivatives
    fig, ax = plt.subplots()
    for i in range(derivs.shape[1]):
        ax.plot(alphas[1:], derivs[:,i], label="Energy level " + str(i))
        pass
    ax.set_xlabel(r"$\alpha$")
    ax.set_ylabel("Derivative of eigenenergy")

    # Locate legend outside of plot
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.show()


    

    return None

if __name__ =="__main__":


    N = 100
    h = (np.pi)/(N-1)
    alphas = np.linspace(0, 0.1, 101)

    Cs = []

    for alpha in alphas:
        C = yga.YgraphAnyonsHardcore(N, alpha)
        states_filepath = "ygraph_hardcore_states/ygraph_N" + str(N) + "_alpha" + str(alpha) + "_"
        eigenvalues_filepath = "ygraph_hardcore_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
        #C.load_states(states_filepath)
        C.load_eigenvalues(eigenvalues_filepath)
        Cs.append(C)
        pass

    # for alpha in alphas:
    #     C = yga.YgraphAnyonsRobin(N, alpha)
    #     C.robin_constant = 0
    #     states_filepath = "ygraph_neumann_states/ygraph_N" + str(N) + "_alpha" + str(alpha) + "_"
    #     eigenvalues_filepath = "ygraph_neumann_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
    #     #C.load_states(states_filepath)
    #     C.load_eigenvalues(eigenvalues_filepath)
    #     Cs.append(C)
    #     pass

    # alphas = np.linspace(0.2, 0.3, 101)
    # for alpha in alphas:
    #     C = yga.YgraphAnyonsHardcore(N, alpha)
    #     states_filepath = "ygraph_zoom_hardcore_states/ygraph_N" + str(N) + "_alpha" + str(alpha) + "_"
    #     eigenvalues_filepath = "ygraph_zoom_hardcore_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
    #     #C.load_states(states_filepath)
    #     C.load_eigenvalues(eigenvalues_filepath)
    #     Cs.append(C)
    #     pass

    # for alpha in alphas:
    #     C = lsa.LassoAnyonsHardcore(N, alpha)
    #     states_filepath = "lasso_hardcore_states/lasso_N" + str(N) + "_alpha" + str(alpha) + "_"
    #     eigenvalues_filepath = "lasso_hardcore_eigenvalues/lasso_N" + str(N) + "_alpha" + str(alpha)
    #     #C.load_states(states_filepath)
    #     C.load_eigenvalues(eigenvalues_filepath)
    #     Cs.append(C)
    #     pass

    # alpha1 = 0.5
    # for alpha2 in alphas:
    #     C = db.DumbbellAnyonsHardcore(N, alpha1, alpha2)
    #     states_filepath = "dumbbell_hardcore_states/dumbbell_N" + str(N) + "_alpha1" + str(alpha1) + "_alpha2" + str(alpha2) + "_"
    #     eigenvalues_filepath = "dumbbell_hardcore_eigenvalues/dumbbell_N" + str(N) + "_alpha1" + str(alpha1) + "_alpha2" + str(alpha2)
    #     C.load_eigenvalues(eigenvalues_filepath)
    #     Cs.append(C)
    #     pass

    #plot_energy_levels(Cs, alphas)
    check_for_stationary_point(Cs, alphas, 0.0)

    pass