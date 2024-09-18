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
    maximum_eigenvalue = 50

    ys = np.array(ys)[:,minimum_eigenvalue:maximum_eigenvalue]

    fig.set_size_inches(6, 10)

    for i in range(ys.shape[1]):
        ax.plot(alphas, ys[:,i])
        #ax.scatter(alphas, ys[:,i])#, label="Energy level " + str(i+1))

    ax.set_xlabel("alpha")
    ax.set_ylabel("Eigenenergy")
    # Locate legend outside of plot
    #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    if show_plots:
        plt.show()
        pass



    return None

if __name__ =="__main__":


    N = 70
    h = (np.pi)/(N-1)
    alphas = np.linspace(0, 1, 25)

    Cs = []

    # for alpha in alphas:
    #     C = yga.YgraphAnyons(N, alpha)
    #     states_filepath = "ygraph_states/ygraph_N" + str(N) + "_alpha" + str(alpha) + "_"
    #     eigenvalues_filepath = "ygraph_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
    #     C.load_states(states_filepath)
    #     C.load_eigenvalues(eigenvalues_filepath)
    #     Cs.append(C)
    #     pass

    # for alpha in alphas:
    #     C = lsa.LassoAnyons(N, alpha)
    #     states_filepath = "lasso_states/lasso_N" + str(N) + "_alpha" + str(alpha) + "_"
    #     eigenvalues_filepath = "lasso_eigenvalues/lasso_N" + str(N) + "_alpha" + str(alpha)
    #     C.load_states(states_filepath)
    #     C.load_eigenvalues(eigenvalues_filepath)
    #     Cs.append(C)
    #     pass

    for alpha in alphas:
        C = db.DumbbellAnyons(N, alpha)
        states_filepath = "dumbbell_states/dumbbell_N" + str(N) + "_alpha" + str(alpha) + "_"
        eigenvalues_filepath = "dumbbell_eigenvalues/dumbbell_N" + str(N) + "_alpha" + str(alpha)
        C.load_states(states_filepath)
        C.load_eigenvalues(eigenvalues_filepath)
        Cs.append(C)
        pass

    plot_energy_levels(Cs, alphas)

    pass