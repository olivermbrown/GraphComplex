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
        #states_filepath = "ygraph_zoom_hardcore_states/ygraph_N" + str(N) + "_alpha" + str(alpha) + "_"
        #eigenvalues_filepath = "ygraph_zoom_hardcore_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
        states_filepath = "ygraph_hardcore_states/ygraph_N" + str(N) + "_alpha" + str(alpha) + "_"
        eigenvalues_filepath = "ygraph_hardcore_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
        C.load_states(states_filepath,override_directory_path=False)
        C.load_eigenvalues(eigenvalues_filepath,override_directory_path=False)
        Cs.append(C)
        pass

    return Cs

def berry_curvature_matrix(C, alpha):

    N = C.cells[0].N

    C.gen_lapl()

    L = C.sL

    # Find Hamiltonian without real part of e^(1j*alpha)

    C_no_real = lsa.LassoAnyonsHardcore(N, 0.5)

    C_no_real.gen_lapl()

    L_no_real = C_no_real.sL

    H = L - np.real(L_no_real)

    # Calculate the Berry curvature matrix

    im_H = np.imag(H)

    mask = np.abs(im_H) > 1e-6

    mask_new = np.abs(im_H) < 1e-6

    H_masked = mask*H

    M = - np.tan(np.pi*alpha) * np.real(H_masked) - 1j * np.imag(H_masked) * 1/(np.tan(np.pi*alpha))

    #C.L = M
    #C.print_eqs()

    print(np.sum(np.abs(M)))

    ss = []
    hs = np.linspace(0.001, 0.01, 10)
    for h in hs:
        #h = 0.001
        Ch = lsa.LassoAnyonsHardcore(N, alpha+h)

        Ch.gen_lapl()

        Hh = Ch.sL

        H_deriv = (Hh - L)/(np.pi*h)

        M_approx = H_deriv*mask

        s = np.sum(np.abs(M_approx - M))
        print(s)
        ss.append(s)
        pass

    plt.plot(hs, ss)
    # X axis label
    plt.xlabel('h')
    # Y axis label
    plt.ylabel('Sum ||M_approx - M||')
    plt.show()



    return M

if __name__ =="__main__":

    N = 30
    h = (np.pi)/(N-1)
    alpha = 0.2

    CY = lsa.LassoAnyonsHardcore(N, alpha)

    M = berry_curvature_matrix(CY, alpha)

    


    #Cs = load_ygraph_wavefunctions(N, alphas)

    #alphas = np.linspace(0,1,25)
    # for alpha in alphas:
    #     C = yga.YgraphAnyonsHardcore(N, alpha)
    #     berry_curvature_matrix(C, alpha)
    #     pass

    # As_15 = calculate_berry_connections(15,Cs, alphas)
    # As_16 = calculate_berry_connections(16,Cs, alphas)

    #As_0 = check_convergence(0, Cs, alphas)
    # As_1 = check_convergence(1, Cs, alphas)

    #save_convergence(N, As_0, N_eig=0)

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