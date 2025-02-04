import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce
import math
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import os

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs
import ygraph as yga
import lasso as lsa
import dumbbell as db

def load_ygraph_wavefunctions(N, alphas):

    Cs = []

    for alpha in alphas:
        C = yga.YgraphAnyonsHardcore(N, alpha)
        states_filepath = "ygraph_zoom_hardcore_states/ygraph_N" + str(N) + "_alpha" + str(alpha) + "_"
        eigenvalues_filepath = "ygraph_zoom_hardcore_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
        #states_filepath = "ygraph_hardcore_states/ygraph_N" + str(N) + "_alpha" + str(alpha) + "_"
        #eigenvalues_filepath = "ygraph_hardcore_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
        C.load_states(states_filepath,override_directory_path=False)
        C.load_eigenvalues(eigenvalues_filepath,override_directory_path=False)
        Cs.append(C)
        pass

    return Cs

def calculate_geometric_tensors(n, Cs, alphas):

    states = []
    for i, C in enumerate(Cs):

        # Find the correct eigenstate
        eigvals_sorted = np.sort(C.spectrum)
        #if i == len(Cs)-1:
        #    idx = np.where(C.spectrum == eigvals_sorted[n+1])[0][0]
        #else:
        idx = np.where(C.spectrum == eigvals_sorted[n])[0][0]

        state = []
        for cell in C.cells:
            state.append(cell.eigenstates[:,idx])
            pass
        state = np.concatenate(state)
        states.append(state)
        pass

    Gs = []
    for i in range(len(Cs)-1):
        delta_alpha = alphas[i+1] - alphas[i]

        # Calculate the geometric tensor
        state1 = states[i]
        state2 = states[i+1]

        # Calculate the overlap
        overlap = np.dot(state1.conj(), state2).flatten()[0]

        # Calculate the geometric tensor
        G = (np.abs(overlap)**2 - 1)/(delta_alpha**2)
        Gs.append(G)
        print(i, "G calculated")
        pass

    return Gs

if __name__ =="__main__":

    N = 100
    h = (np.pi)/(N-1)
    alpha = 1.0

    alphas = np.linspace(0.9,1.0,25)
    #alphas = np.linspace(0, 1, 25)


    Cs = load_ygraph_wavefunctions(N, alphas)

    Gs0 = calculate_geometric_tensors(0, Cs, alphas)
    # Gs1 = calculate_geometric_tensors(1, Cs, alphas)
    # Gs16 = calculate_geometric_tensors(16, Cs, alphas)
    # Gs17 = calculate_geometric_tensors(17, Cs, alphas)

    # Plot the geometric tensors
    fig, ax = plt.subplots()
    fig.set_size_inches(6, 10)
    
    ax.plot(alphas[:-1], Gs0, label="n=0")
    # ax.plot(alphas[:-1], Gs1, label="n=1")
    # ax.plot(alphas[:-1], Gs16, label="n=16")
    # ax.plot(alphas[:-1], Gs17, label="n=17")
    ax.legend()

    ax.set_xlabel("alpha")
    ax.set_ylabel("Geometric tensor")
    plt.tight_layout()

    plt.show()


    pass