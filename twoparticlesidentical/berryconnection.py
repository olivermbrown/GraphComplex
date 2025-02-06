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

    n_eigs = 25

    Cs = []

    for alpha in alphas:
        C = yga.YgraphAnyonsHardcore(N, alpha)
        states_filepath = "ygraph_zoom_hardcore_states/ygraph_N" + str(N) + "_alpha" + str(alpha) + "_"
        eigenvalues_filepath = "ygraph_zoom_hardcore_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
        #states_filepath = "ygraph_hardcore_states/ygraph_N" + str(N) + "_alpha" + str(alpha) + "_"
        #eigenvalues_filepath = "ygraph_hardcore_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
        C.load_states(states_filepath,override_directory_path=True)
        C.load_eigenvalues(eigenvalues_filepath,override_directory_path=True)
        Cs.append(C)
        pass

    return Cs

def calculate_berry_connections(n, Cs, alphas):

    rhos = []
    for i, C in enumerate(Cs):

        # Find the correct eigenstate
        eigvals_sorted = np.sort(C.spectrum)
        idx = np.where(C.spectrum == eigvals_sorted[n])[0][0]

        state = []
        for cell in C.cells:
            state.append(cell.eigenstates[:,idx])
            pass
        state = np.concatenate(state)

        projector = np.outer(state, state.conj())
        
        rhos.append(projector)
        pass

    As = []
    for i in range(len(Cs)-1):
        delta_alpha = alphas[i+1] - alphas[i]
        prod = rhos[i+1] @ rhos[i]
        A = (np.trace(prod)-1)/(2*delta_alpha)
        As.append(A)
        print(i, "A calculated")
        pass
    
    return As

def check_convergence(n, Cs, alphas):

    rhos = []
    for i, C in enumerate(Cs):

        # Find the correct eigenstate
        eigvals_sorted = np.sort(C.spectrum)
        idx = np.where(C.spectrum == eigvals_sorted[n])[0][0]

        state = []
        for cell in C.cells:
            state.append(cell.eigenstates[:,idx])
            pass
        state = np.concatenate(state)

        projector = np.outer(state, state.conj())
        
        rhos.append(projector)
        pass

    As = []
    for i in range(len(Cs)-1):
        delta_alpha = alphas[i+1] - alphas[0]
        prod = rhos[i+1] @ rhos[0]
        A = (np.trace(prod)-1)/(2*delta_alpha)
        As.append(A)
        print(i, "A calculated")
        pass
    
    return As

def check_convergence_new(n, Cs, alphas):

    states = []
    for i, C in enumerate(Cs):

        # Find the correct eigenstate
        eigvals_sorted = np.sort(C.spectrum)
        idx = np.where(C.spectrum == eigvals_sorted[n])[0][0]

        state = []
        for cell in C.cells:
            state.append(cell.eigenstates[:,idx])
            pass
        state = np.concatenate(state)

        states.append(state)
        pass

    As = []
    for i in range(len(Cs)-1):
        delta_alpha = alphas[i+1] - alphas[0]
        prod = np.dot(states[i+1], states[0])
        A = (prod-1)/(2*delta_alpha)
        As.append(A)
        print(i, "A calculated")
        pass
    
    return As

def save_convergence(N, As, N_eig):
    
        folder = "convergence/"
        filename = "convergence_N"+str(N)+"_N_eig"+str(N_eig)
        file_path = folder+filename+"_"
    
        # Create folder if it doesn't exist
        if not os.path.exists(folder):
            try:
                os.makedirs(folder)
                pass
            except FileExistsError:
                pass
            pass
    
        # Save the values to a txt file
        np.savetxt(file_path, As)
    
        return None

def load_convergence(N, N_eig):

    folder = "convergence/"
    filename = "convergence_N"+str(N)+"_N_eig"+str(N_eig)
    file_path = folder+filename+"_"

    As = np.loadtxt(file_path)

    return As

def save_berry_connections(N, As, N_eig):

    folder = "berry_connections/"
    filename = "berry_connections_N"+str(N)+"_N_eig"+str(N_eig)
    file_path = folder+filename+"_"

    # Create folder if it doesn't exist
    if not os.path.exists(folder):
        try:
            os.makedirs(folder)
            pass
        except FileExistsError:
            pass
        pass

    # Save the values to a txt file
    np.savetxt(file_path, As)

    return None

def load_berry_connections(N, N_eig):

    folder = "berry_connections/"
    filename = "berry_connections_N"+str(N)+"_N_eig"+str(N_eig)
    file_path = folder+filename+"_"

    As = np.loadtxt(file_path)

    return As

if __name__ =="__main__":

    N = 50
    h = (np.pi)/(N-1)
    alpha = 1.0

    alphas = np.linspace(0.2,0.3,101)


    Cs = load_ygraph_wavefunctions(N, alphas)

    # As_15 = calculate_berry_connections(15,Cs, alphas)
    # As_16 = calculate_berry_connections(16,Cs, alphas)

    As_0 = check_convergence(0, Cs, alphas)
    # As_1 = check_convergence(1, Cs, alphas)

    save_convergence(N, As_0, N_eig=0)

    #As_conv = load_convergence(N, N_eig=15)

    # save_berry_connections(N, As_15, N_eig=15)
    # save_berry_connections(N, As_16, N_eig=16)

    #As_0 = load_berry_connections(N, N_eig=0)
    #As_1 = load_berry_connections(N, N_eig=1)

    #xs = alphas[:-1] - alphas[0]

    #plt.plot(xs, As_conv, label="n=15")
    #plt.ylabel("Berry connection")
    #plt.xlabel(r"$\Delta \alpha$")
    #plt.plot(alphas[:-1], As_0, label="n=0")
    #plt.plot(alphas[:-1], As_1, label="n=1")
    #plt.legend()
    #plt.show()


    pass