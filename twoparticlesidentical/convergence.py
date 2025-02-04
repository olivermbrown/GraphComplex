import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce
import math
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import os
from itertools import combinations, product, combinations_with_replacement

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs
import ygraph as yga
import lasso as lsa
import dumbbell as db

def check_convergence_YgraphHardcore(alpha, n_eigs):

    Ns = np.arange(10, 160, 10)
    Cs = []
    spectra = []
    for N in Ns:
        C = yga.YgraphAnyonsHardcore(N, alpha)
        eigenvalues_filepath = "ygraph_hardcore_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
        C.load_eigenvalues(eigenvalues_filepath,override_directory_path=False)
        Cs.append(C)
        spectra.append(np.sort(C.spectrum)[:n_eigs])
        pass

    # Plot the spectra as a function of N
    fig, ax = plt.subplots()
    fig.set_size_inches(6, 6)
    for i in range(n_eigs):
        ys = [s[i] for s in spectra]
        ax.plot(Ns, ys)
        pass

    # Calculate the physical spectrum
    single_particle_spectrum = [0.5, 1.0, 1.0, 1.5, 2.0, 2.0, 2.5, 3.0, 3.0, 3.5, 4.0, 4.0, 4.5, 5.0, 5.0, 5.5, 6.0, 6.0, 6.5, 7.0, 7.0, 7.5, 8.0, 8.0, 8.5, 9.0, 9.0]
    two_particle_spectrum = []
    comb_list = list(combinations(single_particle_spectrum,2))
    for pair in comb_list:
        two_particle_spectrum.append(pair[0]**2 + pair[1]**2)
        pass
    two_particle_spectrum = np.sort(two_particle_spectrum)[:n_eigs]
    # Plot physical spectrum with dotted lines
    for i in range(n_eigs):
        ax.axhline(y=two_particle_spectrum[i], color='gray', linestyle='--')
        pass

    ax.set_xlabel("N")
    ax.set_ylabel("Eigenenergy")
    plt.title("Convergence of the spectrum for non-interacting fermions (Dirichlet boundary conditions)", wrap=True)
    plt.tight_layout()
    plt.show()


    # Calculate error when N = 100
    C = Cs[np.where(Ns == 100)[0][0]]
    spec = np.sort(C.spectrum)[:n_eigs]
    error = np.abs(spec - two_particle_spectrum)
    print("Error for N = 100: ", error)

    return None

def check_convergence_YgraphNeumann(alpha, n_eigs):

    Ns = np.arange(10, 160, 10)

    Cs = []
    spectra = []
    for N in Ns:
        C = yga.YgraphAnyonsHardcore(N, alpha)
        eigenvalues_filepath = "ygraph_neumann_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
        C.load_eigenvalues(eigenvalues_filepath,override_directory_path=False)
        Cs.append(C)
        spectra.append(np.sort(C.spectrum)[:n_eigs])
        pass

    # Plot the spectra as a function of N
    fig, ax = plt.subplots()
    fig.set_size_inches(6, 6)
    for i in range(n_eigs):
        ys = [s[i] for s in spectra]
        ax.plot(Ns, ys)
        pass

    # Calculate the physical spectrum
    single_particle_spectrum = [0.5, 1.0, 1.0, 1.5, 2.0, 2.0, 2.5, 3.0, 3.0, 3.5, 4.0, 4.0, 4.5, 5.0, 5.0, 5.5, 6.0, 6.0, 6.5, 7.0, 7.0, 7.5, 8.0, 8.0, 8.5, 9.0, 9.0]
    two_particle_spectrum = []
    comb_list = list(combinations_with_replacement(single_particle_spectrum,2))
    for pair in comb_list:
        two_particle_spectrum.append(pair[0]**2 + pair[1]**2)
        pass
    two_particle_spectrum = np.sort(two_particle_spectrum)[:n_eigs]
    # Plot physical spectrum with dotted lines
    for i in range(n_eigs):
        ax.axhline(y=two_particle_spectrum[i], color='gray', linestyle='--')
        pass

    ax.set_xlabel("N")
    ax.set_ylabel("Eigenenergy")
    plt.title("Convergence of the spectrum for non-interacting bosons (Neumann boundary conditions)", wrap=True)
    plt.tight_layout()
    plt.show()

    # Calculate error when N = 100
    C = Cs[np.where(Ns == 100)[0][0]]
    spec = np.sort(C.spectrum)[:n_eigs]
    error = np.abs(spec - two_particle_spectrum[:n_eigs])
    print("Error for N = 100: ", error)

    return None

if __name__ =="__main__":

    alpha = 1.0

    check_convergence_YgraphHardcore(alpha, 20)


    pass