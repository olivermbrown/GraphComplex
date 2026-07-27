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
import lassosubdivided as lss
import dumbbell as db
import hgraph as hga

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
    fig.set_size_inches(10,8)
    for i in range(n_eigs):
        ys = [s[i] for s in spectra]
        # Interpolate the data to get a smooth curve
        # Use scipy Spline to interpolate the data
        # Interpolate
        f = sp.interpolate.interp1d(Ns, ys, kind='cubic')  # or kind='linear'

        x_new = np.linspace(min(Ns), max(Ns), 200)
        y_new = f(x_new)

        # Plot
        color = f"C{i}"
        ax.plot(Ns, ys, 'o', color=color, markersize=4, label="Numerical data" if i == 0 else "_nolegend_")
        ax.plot(x_new, y_new, '-', color=color, label="Interpolated" if i == 0 else "_nolegend_")
        #ax.plot(Ns, ys)
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
        ax.axhline(y=two_particle_spectrum[i], color='gray', linestyle='--', label="Exact spectrum" if i == 0 else "_nolegend_")
        pass

    ax.set_xlabel("M", fontsize=12)
    ax.set_ylabel("Energy / "+r'$\pi^2$', fontsize=12)
    # Rotate y-axis label to horizontal
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.set_label_coords(-0.08, 0.5)
    plt.legend(loc='lower right', bbox_to_anchor=(1, 0.5))
    #plt.tight_layout()
    plt.savefig("convergence_ygraph_noninteracting_fermions.pdf", dpi=300)
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
    fig.set_size_inches(8, 8)
    for i in range(n_eigs):
        ys = [s[i] for s in spectra]
        f = sp.interpolate.interp1d(Ns, ys, kind='cubic')
        x_new = np.linspace(min(Ns), max(Ns), 200)
        y_new = f(x_new)
        color = f"C{i}"
        ax.plot(x_new, y_new, '-', color=color)
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

    ax.set_xlabel("N", fontsize=12)
    ax.set_ylabel("Energy / "+r'$\pi$', fontsize=12)
    # Rotate y-axis label to horizontal
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.set_label_coords(-0.1, 0.5)
    plt.tight_layout()
    plt.show()

    # Calculate error when N = 100
    C = Cs[np.where(Ns == 100)[0][0]]
    spec = np.sort(C.spectrum)[:n_eigs]
    error = np.abs(spec - two_particle_spectrum[:n_eigs])
    print("Error for N = 100: ", error)

    return None

def check_convergence_HgraphHardcore(alpha, n_eigs):

    Ns = np.arange(10, 160, 10)
    Cs = []
    spectra = []
    for N in Ns:
        C = hga.HgraphAnyonsHardcore(N, alpha, alpha)
        eigenvalues_filepath = "hgraph_hardcore_eigenvalues/hgraph_N" + str(N) + "_alpha1=" + str(alpha) + "_alpha2=" + str(alpha)
        C.load_eigenvalues(eigenvalues_filepath,override_directory_path=False)
        Cs.append(C)
        spectra.append(np.sort(C.spectrum)[:n_eigs])
        pass

    # Plot the spectra as a function of N
    fig, ax = plt.subplots()
    fig.set_size_inches(10,8)
    for i in range(n_eigs):
        ys = [s[i] for s in spectra]
        # Interpolate the data to get a smooth curve
        # Use scipy Spline to interpolate the data
        # Interpolate
        f = sp.interpolate.interp1d(Ns, ys, kind='cubic')  # or kind='linear'

        x_new = np.linspace(min(Ns), max(Ns), 200)
        y_new = f(x_new)

        # Plot
        color = f"C{i}"
        ax.plot(Ns, ys, 'o', color=color, markersize=4, label="Numerical data" if i == 0 else "_nolegend_")
        ax.plot(x_new, y_new, '-', color=color, label="Interpolated" if i == 0 else "_nolegend_")
        #ax.plot(Ns, ys)
        pass

    # Calculate the physical spectrum
    single_particle_spectrum = [0.15353, 0.36987, 1.0, 1.0, 1.0, 1.93718, 2.58622, 4.0, 4.0, 4.0, 5.72083]
    two_particle_spectrum = []
    comb_list = list(combinations(single_particle_spectrum,2))
    for pair in comb_list:
        two_particle_spectrum.append(pair[0] + pair[1])
        pass
    two_particle_spectrum = np.sort(two_particle_spectrum)[:n_eigs]
    # Plot physical spectrum with dotted lines
    for i in range(n_eigs):
        ax.axhline(y=two_particle_spectrum[i], color='gray', linestyle='--', label="Exact spectrum" if i == 0 else "_nolegend_")
        pass

    ax.set_xlabel("M", fontsize=12)
    ax.set_ylabel("Energy / "+r'$\pi^2$', fontsize=12)
    # Rotate y-axis label to horizontal
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.set_label_coords(-0.08, 0.5)
    plt.legend(loc='lower right', bbox_to_anchor=(1, 0.5))
    #plt.tight_layout()
    plt.savefig("convergence_hgraph_noninteracting_fermions.pdf", dpi=300)
    plt.show()


    # Calculate error when N = 100
    C = Cs[np.where(Ns == 100)[0][0]]
    spec = np.sort(C.spectrum)[:n_eigs]
    error = np.abs(spec - two_particle_spectrum)
    print("Error for N = 100: ", error)

    return None

def check_convergence_LassosubdividedHardcore(alpha, n_eigs):

    Ns = np.arange(10, 160, 10)
    Cs = []
    spectra = []
    aharonovbohm = 0.0
    for N in Ns:
        C = lss.LassoSubdividedAnyonsHardcore(N, N, alpha, aharonovbohm)
        eigenvalues_filepath = "lassosubdivided_hardcore_eigenvalues/lassosubdivided_N" + str(N) + "_alpha" + str(alpha) + "_aharonovbohm" + str(aharonovbohm)
        C.load_eigenvalues(eigenvalues_filepath,override_directory_path=False)
        Cs.append(C)
        spectra.append(np.sort(C.spectrum)[:n_eigs])
        pass

    # Plot the spectra as a function of N
    fig, ax = plt.subplots()
    fig.set_size_inches(10,8)
    for i in range(n_eigs):
        ys = [s[i] for s in spectra]
        # Interpolate the data to get a smooth curve
        # Use scipy Spline to interpolate the data
        # Interpolate
        f = sp.interpolate.interp1d(Ns, ys, kind='cubic')  # or kind='linear'

        x_new = np.linspace(min(Ns), max(Ns), 200)
        y_new = f(x_new)

        # Plot
        color = f"C{i}"
        ax.plot(Ns, ys, 'o', color=color, markersize=4, label="Numerical data" if i == 0 else "_nolegend_")
        ax.plot(x_new, y_new, '-', color=color, label="Interpolated" if i == 0 else "_nolegend_")
        #ax.plot(Ns, ys)
        pass

    # Calculate the physical spectrum
    single_particle_spectrum = [0.03838, 0.64656, 1, 1.43021, 3.25473, 4, 4.82204, 7.86290, 9]
    two_particle_spectrum = []
    comb_list = list(combinations(single_particle_spectrum,2))
    for pair in comb_list:
        two_particle_spectrum.append(pair[0] + pair[1])
        pass
    two_particle_spectrum = np.sort(two_particle_spectrum)[:n_eigs]
    # Plot physical spectrum with dotted lines
    for i in range(n_eigs):
        ax.axhline(y=two_particle_spectrum[i], color='gray', linestyle='--', label="Exact spectrum" if i == 0 else "_nolegend_")
        pass

    ax.set_xlabel("M", fontsize=12)
    ax.set_ylabel("Energy / "+r'$\pi^2$', fontsize=12)
    # Rotate y-axis label to horizontal
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.set_label_coords(-0.08, 0.5)
    plt.legend(loc='lower right', bbox_to_anchor=(1, 0.5))
    #plt.tight_layout()
    plt.savefig("convergence_lassosubdivided_noninteracting_fermions.pdf", dpi=300)
    plt.show()


    # Calculate error when N = 100
    C = Cs[np.where(Ns == 100)[0][0]]
    spec = np.sort(C.spectrum)[:n_eigs]
    error = np.abs(spec - two_particle_spectrum)
    print("Error for N = 100: ", error)

    return None

if __name__ =="__main__":

    alpha = 1.0

    check_convergence_LassosubdividedHardcore(alpha, 20)


    pass