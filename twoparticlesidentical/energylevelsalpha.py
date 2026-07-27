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
import lassosubdivided as lss
import dumbbell as db
import dumbbellsubdivided as dbs

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
    maximum_eigenvalue = 20

    ys = np.array(ys)[:,minimum_eigenvalue:maximum_eigenvalue]

    fig.set_size_inches(10, 8)

    for i in range(ys.shape[1]):
        #ax.plot(alphas, ys[:,i])
        ax.plot(alphas, ys[:,i])

    ax.set_xlabel(r"$\alpha$", fontsize=12)
    ax.set_ylabel("Energy / "+r'$\pi^2$', fontsize=12)
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.set_label_coords(-0.08, 0.5)
    # Locate legend outside of plot
    #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    if show_plots:
        plt.show()
        pass

    # Get user input for filename to save the plot
    filename = input("Enter filename to save the plot (or press Enter to skip saving): ")
    if filename:
        fig.savefig(filename+ ".pdf", dpi=300)
        print(f"Plot saved as {filename}")
        pass
    else:
        print("Plot not saved.")
        pass



    return None

def plot_energy_levels_3D(Cs, alphas, show_plots=True):
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
    maximum_eigenvalue = 1

    ys = np.array(ys)[:,minimum_eigenvalue:maximum_eigenvalue]

    fig.set_size_inches(6, 10)
    ax = plt.axes(projection='3d')

    for i in range(ys.shape[1]):

        # Create 3D plot

        # Create a meshgrid of alphas
        X, Y = np.meshgrid(alphas, alphas)
        Z = ys[:,i].T

        # Reshape Z to be 2D
        Z = np.reshape(Z, (len(alphas), len(alphas)))

        # Plot the surface
        ax.plot_surface(X, Y, Z, label="Energy level " + str(i+minimum_eigenvalue), alpha=0.5)

        pass

    ax.set_xlabel(r"$\alpha1$", fontsize=12)
    ax.set_ylabel(r'$\alpha2$', fontsize=12)
    ax.set_zlabel("Ground State Energy", fontsize=12)
    # Locate legend outside of plot
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    #plt.tight_layout()
    if show_plots:
        plt.show()
        pass



    return None

if __name__ =="__main__":


    N = 100
    h = (np.pi)/(N-1)
    alphas = np.linspace(0, 1, 25)

    Cs = []

    # for alpha in alphas:
    #     C = yga.YgraphAnyonsHardcore(N, alpha)
    #     states_filepath = "ygraph_hardcore_states/ygraph_N" + str(N) + "_alpha" + str(alpha) + "_"
    #     eigenvalues_filepath = "ygraph_hardcore_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
    #     #C.load_states(states_filepath)
    #     C.load_eigenvalues(eigenvalues_filepath)
    #     Cs.append(C)
    #     pass

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
    #     C = lss.LassoSubdividedAnyonsHardcore(N, N, alpha, 0.0)
    #     states_filepath = "lassosubdivided_hardcore_states/lassosubdivided_N" + str(N) + "_alpha" + str(alpha) + "_aharonovbohm0.0_"
    #     eigenvalues_filepath = "lassosubdivided_hardcore_eigenvalues/lassosubdivided_N" + str(N) + "_alpha" + str(alpha) + "_aharonovbohm0.0"
    #     #C.load_states(states_filepath)
    #     C.load_eigenvalues(eigenvalues_filepath)
    #     Cs.append(C)
    #     pass

    for alpha in alphas:
        C = dbs.DumbbellSubdividedAnyonsHardcore(N, N, N, alpha, alpha, 0.0, 0.0)
        states_filepath = "dumbbellsubdivided_hardcore_states/dumbbellsubdivided_N" + str(N) + "_alpha1=" + str(alpha) + "_alpha2=" + str(alpha) + "_aharonovbohm1=0.0_" + "_aharonovbohm2=0.0_"
        eigenvalues_filepath = "dumbbellsubdivided_hardcore_eigenvalues/dumbbellsubdivided_N" + str(N) + "_alpha1=" + str(alpha) + "_alpha2=" + str(alpha) + "_aharonovbohm1=0.0" + "_aharonovbohm2=0.0"
        #C.load_states(states_filepath)
        C.load_eigenvalues(eigenvalues_filepath)
        Cs.append(C)
        pass

    # Create list of pairs of alphas
    # alphas_pairs = [(alpha1, alpha2) for alpha1 in alphas for alpha2 in alphas]
    # for (alpha1, alpha2) in alphas_pairs:
    #     C = db.DumbbellAnyonsHardcore(N, alpha1, alpha2)
    #     eigenvalues_filepath = "dumbbell_hardcore_eigenvalues/dumbbell_N" + str(N) + "_alpha1" + str(alpha1) + "_alpha2" + str(alpha2)
    #     C.load_eigenvalues(eigenvalues_filepath, override_directory_path=True)
    #     Cs.append(C)
    #     pass

    plot_energy_levels(Cs, alphas)
    #plot_energy_levels_3D(Cs, alphas)

    #fig = plt.gcf()
    # Save plot
    #fig.savefig("energy_levels_ygraph_anyons_hardcore.pdf", dpi=300)

    pass