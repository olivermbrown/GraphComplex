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

# Change the default font family globally
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

import cells as cls
import configs
import ygraph as yga
import lasso as lsa
import lassosubdivided as lss
import dumbbell as db
import hgraph as hga
import cross as cga

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
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
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
        ax.plot(Ns, ys, 'o', color=color, markersize=8, label="Numerical data" if i == 0 else "_nolegend_")
        ax.plot(x_new, y_new, '-', color=color, label="Interpolated" if i == 0 else "_nolegend_", linewidth=2)
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
        ax.axhline(y=two_particle_spectrum[i], color='gray', linestyle='--', label="Exact spectrum" if i == 0 else "_nolegend_", linewidth=2)
        pass

    # Set tick sizes
    ax.tick_params(axis='both', which='both', labelsize=16)
    ax.set_xlabel("M", fontsize=20, labelpad=12)
    ax.set_ylabel("Energy / "+r'$\pi^2$', fontsize=20)
    # Rotate y-axis label to horizontal
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.set_label_coords(-0.2, 0.5)
    ax.yaxis.label.set_horizontalalignment('right')
    ax.yaxis.label.set_verticalalignment('center')
    plt.legend(loc='center right', fontsize=20)
    #plt.tight_layout()
    plt.savefig("convergence_Ygraph_noninteracting_fermions.pdf", dpi=300, bbox_inches='tight', pad_inches=0.02)
    plt.show()


    # Calculate error when N = 100
    C = Cs[np.where(Ns == 100)[0][0]]
    spec = np.sort(C.spectrum)[:n_eigs]
    error = np.abs(spec - two_particle_spectrum)
    print("Error for N = 100: ", error)

    return None

def log_convergence_YgraphHardcore_single_eigval(alpha, n_eigval):

    Ns = np.arange(10, 160, 10)
    Cs = []
    eigvals = []
    for N in Ns:
        C = yga.YgraphAnyonsHardcore(N, alpha)
        eigenvalues_filepath = "ygraph_hardcore_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
        C.load_eigenvalues(eigenvalues_filepath,override_directory_path=False)
        Cs.append(C)
        eigvals.append(np.sort(C.spectrum)[n_eigval])
        pass

        # Calculate the physical spectrum
    single_particle_spectrum = [0.5, 1.0, 1.0, 1.5, 2.0, 2.0, 2.5, 3.0, 3.0, 3.5, 4.0, 4.0, 4.5, 5.0, 5.0, 5.5, 6.0, 6.0, 6.5, 7.0, 7.0, 7.5, 8.0, 8.0, 8.5, 9.0, 9.0]
    two_particle_spectrum = []
    comb_list = list(combinations(single_particle_spectrum,2))
    for pair in comb_list:
        two_particle_spectrum.append(pair[0]**2 + pair[1]**2)
        pass
    exact_eigval = np.sort(two_particle_spectrum)[n_eigval]

    # Calculate the error for each N
    errors = [np.abs(eigval - exact_eigval) for eigval in eigvals]

    hs = [1.0/N for N in Ns]

    # Plot the nth eigenvalue as a function of N
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)

    # Choose the range over which to fit
    fit_end = 10     # use first 10 data points

    logh = np.log(hs[:fit_end])
    loge = np.log(errors[:fit_end])

    # Linear fit in log-log space
    coeffs = np.polyfit(logh, loge, 1)
    m, c = coeffs

    # Plot the numerical erros
    ax.loglog(hs, errors, 'D', label='Numerical error', color='darkblue', markersize=8)

    # Plot fitted line ONLY over the fitted region
    h_fit = hs[:fit_end]
    error_fit = np.exp(c) * h_fit**m

    ax.loglog(h_fit, error_fit, '--',
            label=f'Fit (gradient = {m:.2f})', color='orange', linewidth=2)

    # # Shade the implementation-limited region
    # plt.axvspan(hs[9], hs[-1], alpha=0.15, color='grey',
    #         label='Implementation-limited regime')

    ax.legend(loc='center right', fontsize=20)

    # Set tick sizes
    ax.tick_params(axis='both', which='both', labelsize=16)

    ax.set_xlabel(r"$h$", fontsize=20, labelpad=12)
    ax.set_ylabel(r"$|E_0^h - E_0| \, / \, \pi^2$",
              fontsize=20)
    # Rotate y-axis label to horizontal
    ax.yaxis.label.set_rotation(0)
    # Adjust the position of the y-axis label
    ax.yaxis.label.set_horizontalalignment('right')
    ax.yaxis.label.set_verticalalignment('center')
    ax.yaxis.set_label_coords(-0.2, 0.5)
    ax.legend(loc='center right', fontsize=20)
    filename = f"log_convergence_Ygraph_noninteracting_fermions_eigval{n_eigval}.pdf"
    plt.savefig(filename, dpi=300, bbox_inches='tight', pad_inches=0.02)
    plt.show()

    # Print all errors
    print("Errors for N = 10, 20, ..., 150: ", errors)


    return None

def log_convergence_YgraphHardcore(alpha, n_eigs, fit_end=10):

    Ns = np.arange(10, 160, 10)
    hs = np.array([1.0/N for N in Ns])

    # ------------------------------------------------------------------
    # Load numerical spectra
    # ------------------------------------------------------------------
    spectra = []

    for N in Ns:
        C = yga.YgraphAnyonsHardcore(N, alpha)
        eigenvalues_filepath = (
            f"ygraph_hardcore_eigenvalues/ygraph_N{N}_alpha{alpha}"
        )
        C.load_eigenvalues(eigenvalues_filepath,
                           override_directory_path=False)

        spectra.append(np.sort(C.spectrum))

    # ------------------------------------------------------------------
    # Exact two-particle spectrum
    # ------------------------------------------------------------------
    single_particle_spectrum = [
        0.5, 1.0, 1.0, 1.5, 2.0, 2.0, 2.5, 3.0, 3.0,
        3.5, 4.0, 4.0, 4.5, 5.0, 5.0, 5.5, 6.0, 6.0,
        6.5, 7.0, 7.0, 7.5, 8.0, 8.0, 8.5, 9.0, 9.0
    ]

    two_particle_spectrum = np.sort([
        a**2 + b**2
        for a, b in combinations(single_particle_spectrum, 2)
    ])

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)

    markers = ['D', 'o', 's', '^', 'v', '<', '>', 'p', '*', 'X']

    for n in range(n_eigs):

        exact = two_particle_spectrum[n]

        errors = np.array([
            abs(spec[n] - exact)
            for spec in spectra
        ])

        # Fit over first fit_end points
        coeffs = np.polyfit(
            np.log(hs[:fit_end]),
            np.log(errors[:fit_end]),
            1
        )

        m, c = coeffs

        # Numerical data
        ax.loglog(
            hs,
            errors,
            marker=markers[n % len(markers)],
            linestyle='None',
            markersize=8,
            label=fr"$\lambda_{{{n}}}$"
        )

        # Fitted line
        h_fit = hs[:fit_end]
        error_fit = np.exp(c) * h_fit**m

        ax.loglog(
            h_fit,
            error_fit,
            '--',
            linewidth=1.5,
            label=fr"Fit $\lambda_{{{n}}}$ ({m:.2f})"
        )

    ax.set_xlabel(r"$h$", fontsize=20, labelpad=12)

    ax.set_ylabel(
        r"$|E_0^h - E_0| \, / \, \pi^2$",
        fontsize=20,
        labelpad=12
    )

    ax.yaxis.label.set_rotation(0)
    ax.yaxis.label.set_horizontalalignment('right')
    ax.yaxis.label.set_verticalalignment('center')
    ax.yaxis.set_label_coords(-0.2, 0.5)

    ax.legend(fontsize=20, ncol=2)

    plt.savefig(
        "log_convergence_Ygraph_noninteracting_fermions.pdf",
        dpi=300,
        bbox_inches='tight',
        pad_inches=0.02
    )

    plt.show()

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
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
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

    ax.set_xlabel("N", fontsize=20)
    ax.set_ylabel("Energy / "+r'$\pi$', fontsize=20)
    # Rotate y-axis label to horizontal
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.set_label_coords(-0.2, 0.5)
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
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    
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
        ax.plot(Ns, ys, 'o', color=color, markersize=8, label="Numerical data" if i == 0 else "_nolegend_")
        ax.plot(x_new, y_new, '-', color=color, label="Interpolated" if i == 0 else "_nolegend_", linewidth=2)
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
        ax.axhline(y=two_particle_spectrum[i], color='gray', linestyle='--', label="Exact spectrum" if i == 0 else "_nolegend_", linewidth=2)
        pass

    # Set tick sizes
    ax.tick_params(axis='both', which='both', labelsize=16)

    ax.set_xlabel("M", fontsize=20, labelpad=12)
    ax.set_ylabel("Energy / "+r'$\pi^2$', fontsize=20)
    # Rotate y-axis label to horizontal
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.set_label_coords(-0.2, 0.5)
    ax.yaxis.label.set_horizontalalignment('right')
    ax.yaxis.label.set_verticalalignment('center')
    plt.legend(loc='center right', fontsize=20)
    #plt.tight_layout()
    plt.savefig("convergence_Hgraph_noninteracting_fermions.pdf", dpi=300, bbox_inches='tight', pad_inches=0.02)
    plt.show()


    # Calculate error when N = 100
    C = Cs[np.where(Ns == 100)[0][0]]
    spec = np.sort(C.spectrum)[:n_eigs]
    error = np.abs(spec - two_particle_spectrum)
    print("Error for N = 100: ", error)

    return None

def log_convergence_HgraphHardcore_single_eigval(alpha, n_eigval):
    Ns = np.arange(10, 160, 10)
    Cs = []
    eigvals = []
    for N in Ns:
        C = hga.HgraphAnyonsHardcore(N, alpha, alpha)
        eigenvalues_filepath = "hgraph_hardcore_eigenvalues/hgraph_N" + str(N) + "_alpha1=" + str(alpha) + "_alpha2=" + str(alpha)
        C.load_eigenvalues(eigenvalues_filepath,override_directory_path=False)
        Cs.append(C)
        eigvals.append(np.sort(C.spectrum)[n_eigval])
        pass

    # Calculate the physical spectrum
    single_particle_spectrum = [0.15353, 0.36987, 1.0, 1.0, 1.0, 1.93718, 2.58622, 4.0, 4.0, 4.0, 5.72083]
    two_particle_spectrum = []
    comb_list = list(combinations(single_particle_spectrum,2))
    for pair in comb_list:
        two_particle_spectrum.append(pair[0] + pair[1])
        pass
    exact_eigval = np.sort(two_particle_spectrum)[n_eigval]

    # Calculate the error for each N
    errors = [np.abs(eigval - exact_eigval) for eigval in eigvals]

    hs = [1.0/N for N in Ns]

    # Plot the nth eigenvalue as a function of N
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)

    # Choose the range over which to fit
    fit_end = 15     # use first 10 data points

    logh = np.log(hs[:fit_end])
    loge = np.log(errors[:fit_end])

    # Linear fit in log-log space
    coeffs = np.polyfit(logh, loge, 1)
    m, c = coeffs

    # Plot the numerical erros
    ax.loglog(hs, errors, 'D', label='Numerical error', color='darkblue', markersize=8)

    # Plot fitted line ONLY over the fitted region
    h_fit = hs[:fit_end]
    error_fit = np.exp(c) * h_fit**m

    ax.loglog(h_fit, error_fit, '--', label=f'Fitted line (m={m:.2f})', color='orange', linewidth=2)

    # Set tick sizes
    ax.tick_params(axis='both', which='both', labelsize=16)

    ax.set_xlabel(r"$h$", fontsize=20, labelpad=12)
    ax.set_ylabel(r"$|E_0^h - E_0| \, / \, \pi^2$",
              fontsize=20,
              labelpad=12)
    # Rotate y-axis label to horizontal
    ax.yaxis.label.set_rotation(0)
    # Adjust the position of the y-axis label
    ax.yaxis.label.set_horizontalalignment('right')
    ax.yaxis.label.set_verticalalignment('center')
    ax.yaxis.set_label_coords(-0.2, 0.5)
    ax.legend(loc='center right', fontsize=20)
    #plt.tight_layout()
    plt.savefig(f"log_convergence_Hgraph_noninteracting_fermions_eigval{n_eigval}.pdf", dpi=300, bbox_inches='tight', pad_inches=0.02)
    plt.show()

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
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    
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
        ax.plot(Ns, ys, 'o', color=color, markersize=8, label="Numerical data" if i == 0 else "_nolegend_")
        ax.plot(x_new, y_new, '-', color=color, label="Interpolated" if i == 0 else "_nolegend_", linewidth=2)
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
        ax.axhline(y=two_particle_spectrum[i], color='gray', linestyle='--', label="Exact spectrum" if i == 0 else "_nolegend_", linewidth=2)
        pass

    # Set tick sizes
    ax.tick_params(axis='both', which='both', labelsize=16)

    ax.set_xlabel("M", fontsize=20, labelpad=12)
    ax.set_ylabel("Energy / "+r'$\pi^2$', fontsize=20)
    # Rotate y-axis label to horizontal
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.set_label_coords(-0.2, 0.5)
    ax.yaxis.label.set_horizontalalignment('right')
    ax.yaxis.label.set_verticalalignment('center')
    plt.legend(loc='center right', fontsize=20)
    #plt.tight_layout()
    plt.savefig("convergence_Lassosubdivided_noninteracting_fermions.pdf", dpi=300, bbox_inches='tight', pad_inches=0.02)
    plt.show()


    # Calculate error when N = 100
    C = Cs[np.where(Ns == 100)[0][0]]
    spec = np.sort(C.spectrum)[:n_eigs]
    error = np.abs(spec - two_particle_spectrum)
    print("Error for N = 100: ", error)

    return None

def log_convergence_LassosubdividedHardcore_single_eigval(alpha, n_eigval):

    Ns = np.arange(10, 160, 10)
    Cs = []
    eigvals = []
    aharonovbohm = 0.0
    for N in Ns:
        C = lss.LassoSubdividedAnyonsHardcore(N, N, alpha, aharonovbohm)
        eigenvalues_filepath = "lassosubdivided_hardcore_eigenvalues/lassosubdivided_N" + str(N) + "_alpha" + str(alpha) + "_aharonovbohm" + str(aharonovbohm)
        C.load_eigenvalues(eigenvalues_filepath,override_directory_path=False)
        Cs.append(C)
        eigvals.append(np.sort(C.spectrum)[n_eigval])
        pass

    # Calculate the physical spectrum
    single_particle_spectrum = [0.03838, 0.64656, 1, 1.43021, 3.25473, 4, 4.82204, 7.86290, 9]
    two_particle_spectrum = []
    comb_list = list(combinations(single_particle_spectrum,2))
    for pair in comb_list:
        two_particle_spectrum.append(pair[0] + pair[1])
        pass
    exact_eigval = np.sort(two_particle_spectrum)[n_eigval]

    # Calculate the error for each N
    errors = [np.abs(eigval - exact_eigval) for eigval in eigvals]

    hs = [1.0/N for N in Ns]

    # Plot the nth eigenvalue as a function of N
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)

    # Choose the range over which to fit
    fit_end = 15     # use first 10 data points

    logh = np.log(hs[:fit_end])
    loge = np.log(errors[:fit_end])

    # Linear fit in log-log space
    coeffs = np.polyfit(logh, loge, 1)
    m, c = coeffs

    # Plot the numerical erros
    ax.loglog(hs, errors, 'D', label='Numerical error', color='darkblue', markersize=8)

    # Plot fitted line ONLY over the fitted region
    h_fit = hs[:fit_end]
    error_fit = np.exp(c) * h_fit**m

    ax.loglog(h_fit, error_fit, '--', label=f'Fitted line (m={m:.2f})', color='orange', linewidth=2)

    # Set tick sizes
    ax.tick_params(axis='both', which='both', labelsize=16)
    
    ax.set_xlabel(r"$h$", fontsize=20, labelpad=12)
    ax.set_ylabel(r"$|E_0^h - E_0| \, / \, \pi^2$",
              fontsize=20,
              labelpad=12)
    # Rotate y-axis label to horizontal
    ax.yaxis.label.set_rotation(0)
    # Adjust the position of the y-axis label
    ax.yaxis.label.set_horizontalalignment('right')
    ax.yaxis.label.set_verticalalignment('center')
    ax.yaxis.set_label_coords(-0.2, 0.5)
    ax.legend(loc='center right', fontsize=20)
    #plt.tight_layout()
    plt.savefig(f"log_convergence_Lassosubdivided_noninteracting_fermions_eigval{n_eigval}.pdf", dpi=300, bbox_inches='tight', pad_inches=0.02)
    plt.show()

    return None

def check_convergence_DumbbellsubdividedHardcore(alpha, n_eigs):

    Ns = np.arange(10, 160, 10)
    Cs = []
    spectra = []
    for N in Ns:
        C = db.DumbbellAnyonsHardcore(N, alpha, alpha)
        eigenvalues_filepath = "dumbbellsubdivided_hardcore_eigenvalues/dumbbellsubdivided_N" + str(N) + "_alpha1=" + str(alpha) + "_alpha2=" + str(alpha) + "_aharonovbohm1=0.0" + "_aharonovbohm2=0.0"
        C.load_eigenvalues(eigenvalues_filepath,override_directory_path=False)
        Cs.append(C)
        spectra.append(np.sort(C.spectrum)[:n_eigs])
        pass

    # Plot the spectra as a function of N
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    
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
        ax.plot(Ns, ys, 'o', color=color, markersize=8, label="Numerical data" if i == 0 else "_nolegend_")
        ax.plot(x_new, y_new, '-', color=color, label="Interpolated" if i == 0 else "_nolegend_", linewidth=2)
        #ax.plot(Ns, ys)
        pass

    # Calculate the physical spectrum
    single_particle_spectrum = [0.0, 0.07167, 0.53623, 1.0, 1.0, 1.0, 1.60712, 3.00079, 4.0, 4.0, 4.0, 5.14256, 7.46535, 9.0, 9.0, 9.0]
    two_particle_spectrum = []
    comb_list = list(combinations(single_particle_spectrum,2))
    for pair in comb_list:
        two_particle_spectrum.append(pair[0] + pair[1])
        pass
    two_particle_spectrum = np.sort(two_particle_spectrum)[:n_eigs]
    # Plot physical spectrum with dotted lines
    for i in range(n_eigs):
        ax.axhline(y=two_particle_spectrum[i], color='gray', linestyle='--', label="Exact spectrum" if i == 0 else "_nolegend_", linewidth=2)
        pass

    # Set tick sizes
    ax.tick_params(axis='both', which='both', labelsize=16)

    ax.set_xlabel("M", fontsize=20, labelpad=12)
    ax.set_ylabel("Energy / "+r'$\pi^2$', fontsize=20)
    # Rotate y-axis label to horizontal
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.set_label_coords(-0.2, 0.5)
    ax.yaxis.label.set_horizontalalignment('right')
    ax.yaxis.label.set_verticalalignment('center')
    plt.legend(loc='center right', fontsize=20)
    #plt.tight_layout()
    plt.savefig("convergence_Dumbbellsubdivided_noninteracting_fermions.pdf", dpi=300, bbox_inches='tight', pad_inches=0.02)
    plt.show()

    # Calculate error when N = 100
    C = Cs[np.where(Ns == 100)[0][0]]
    spec = np.sort(C.spectrum)[:n_eigs]
    error = np.abs(spec - two_particle_spectrum)
    print("Error for N = 100: ", error)

    return None

def log_convergence_DumbbellsubdividedHardcore_single_eigval(alpha, n_eigval):

    Ns = np.arange(10, 160, 10)
    Cs = []
    eigvals = []
    for N in Ns:
        C = db.DumbbellAnyonsHardcore(N, alpha, alpha)
        eigenvalues_filepath = "dumbbellsubdivided_hardcore_eigenvalues/dumbbellsubdivided_N" + str(N) + "_alpha1=" + str(alpha) + "_alpha2=" + str(alpha) + "_aharonovbohm1=0.0" + "_aharonovbohm2=0.0"
        C.load_eigenvalues(eigenvalues_filepath,override_directory_path=False)
        Cs.append(C)
        eigvals.append(np.sort(C.spectrum)[n_eigval])
        pass

    # Calculate the physical spectrum
    single_particle_spectrum = [0.0, 0.07167, 0.53623, 1.0, 1.0, 1.0, 1.60712, 3.00079, 4.0, 4.0, 4.0, 5.14256, 7.46535, 9.0, 9.0, 9.0]
    two_particle_spectrum = []
    comb_list = list(combinations(single_particle_spectrum,2))
    for pair in comb_list:
        two_particle_spectrum.append(pair[0] + pair[1])
        pass
    exact_eigval = np.sort(two_particle_spectrum)[n_eigval]

    # Calculate the error for each N
    errors = [np.abs(eigval - exact_eigval) for eigval in eigvals]

    hs = [1.0/N for N in Ns]

    # Plot the nth eigenvalue as a function of N
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)

    # Choose the range over which to fit
    fit_end = 15     # use first 10 data points

    logh = np.log(hs[:fit_end])
    loge = np.log(errors[:fit_end])

    # Linear fit in log-log space
    coeffs = np.polyfit(logh, loge, 1)
    m, c = coeffs

    # Plot the numerical erros
    ax.loglog(hs, errors, 'D', label='Numerical error', color='darkblue', markersize=8)

    # Plot fitted line ONLY over the fitted region
    h_fit = hs[:fit_end]
    error_fit = np.exp(c) * h_fit**m

    ax.loglog(h_fit, error_fit, '--', label=f'Fitted line (m={m:.2f})', color='orange', linewidth=2)

    # Set tick sizes
    ax.tick_params(axis='both', which='both', labelsize=16)

    ax.set_xlabel(r"$h$", fontsize=20, labelpad=12)
    ax.set_ylabel(r"$|E_0^h - E_0| \, / \, \pi^2$",
              fontsize=20,
                labelpad=12)
    # Rotate y-axis label to horizontal
    ax.yaxis.label.set_rotation(0)
    # Adjust the position of the y-axis label
    ax.yaxis.label.set_horizontalalignment('right')
    ax.yaxis.label.set_verticalalignment('center')
    ax.yaxis.set_label_coords(-0.2, 0.5)
    ax.legend(loc='center right', fontsize=20)
    #plt.tight_layout()
    plt.savefig(f"log_convergence_Dumbbellsubdivided_noninteracting_fermions_eigval{n_eigval}.pdf", dpi=300, bbox_inches='tight', pad_inches=0.02)
    plt.show()

    return None

def check_convergence_CrossgraphHardcore(alpha, n_eigs):

    Ns = np.arange(10, 160, 10)
    Cs = []
    spectra = []
    for N in Ns:
        C = cga.CrossgraphAnyonsHardcore(N, alpha, alpha, alpha)
        eigenvalues_filepath = "cross_all_alphas_varying_hardcore_eigenvalues/cross_N" + str(N) + "_alpha" + str(alpha)
        C.load_eigenvalues(eigenvalues_filepath,override_directory_path=False)
        Cs.append(C)
        spectra.append(np.sort(C.spectrum)[:n_eigs])
        pass

    # Plot the spectra as a function of N
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    
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
        ax.plot(Ns, ys, 'o', color=color, markersize=8, label="Numerical data" if i == 0 else "_nolegend_")
        ax.plot(x_new, y_new, '-', color=color, label="Interpolated" if i == 0 else "_nolegend_", linewidth=2)
        #ax.plot(Ns, ys)
        pass

    # Calculate the physical spectrum
    single_particle_spectrum = [0.25, 1.0, 1.0, 1.0, 2.25, 4.0, 4.0, 4.0, 6.25, 9.0, 9.0, 9.0]
    two_particle_spectrum = []
    comb_list = list(combinations(single_particle_spectrum,2))
    for pair in comb_list:
        two_particle_spectrum.append(pair[0] + pair[1])
        pass
    two_particle_spectrum = np.sort(two_particle_spectrum)[:n_eigs]
    # Plot physical spectrum with dotted lines
    for i in range(n_eigs):
        ax.axhline(y=two_particle_spectrum[i], color='gray', linestyle='--', label="Exact spectrum" if i == 0 else "_nolegend_", linewidth=2)
        pass

    # Set tick sizes
    ax.tick_params(axis='both', which='both', labelsize=16)
    
    ax.set_xlabel("M", fontsize=20, labelpad=12)
    ax.set_ylabel("Energy / "+r'$\pi^2$', fontsize=20)
    # Rotate y-axis label to horizontal
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.set_label_coords(-0.2, 0.5)
    ax.yaxis.label.set_horizontalalignment('right')
    ax.yaxis.label.set_verticalalignment('center')
    plt.legend(loc='center right', fontsize=20)
    #plt.tight_layout()
    plt.savefig("convergence_Crossgraph_noninteracting_fermions.pdf", dpi=300, bbox_inches='tight', pad_inches=0.02)
    plt.show()

    # Calculate error when N = 100
    C = Cs[np.where(Ns == 100)[0][0]]
    spec = np.sort(C.spectrum)[:n_eigs]
    error = np.abs(spec - two_particle_spectrum)
    print("Error for N = 100: ", error)

    return None

def log_convergence_CrossgraphHardcore_single_eigval(alpha, n_eigval):

    Ns = np.arange(10, 160, 10)
    Cs = []
    eigvals = []
    for N in Ns:
        C = cga.CrossgraphAnyonsHardcore(N, alpha, alpha, alpha)
        eigenvalues_filepath = "cross_all_alphas_varying_hardcore_eigenvalues/cross_N" + str(N) + "_alpha" + str(alpha)
        C.load_eigenvalues(eigenvalues_filepath,override_directory_path=False)
        Cs.append(C)
        eigvals.append(np.sort(C.spectrum)[n_eigval])
        pass

    # Calculate the physical spectrum
    single_particle_spectrum = [0.25, 1.0, 1.0, 1.0, 2.25, 4.0, 4.0, 4.0, 6.25, 9.0, 9.0, 9.0]
    two_particle_spectrum = []
    comb_list = list(combinations(single_particle_spectrum,2))
    for pair in comb_list:
        two_particle_spectrum.append(pair[0] + pair[1])
        pass
    exact_eigval = np.sort(two_particle_spectrum)[n_eigval]

    # Calculate the error for each N
    errors = [np.abs(eigval - exact_eigval) for eigval in eigvals]

    hs = [1.0/N for N in Ns]

    # Plot the nth eigenvalue as a function of N
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)

    # Choose the range over which to fit
    fit_end = 15     # use first 10 data points

    logh = np.log(hs[:fit_end])
    loge = np.log(errors[:fit_end])

    # Linear fit in log-log space
    coeffs = np.polyfit(logh, loge, 1)
    m, c = coeffs

    # Plot the numerical erros
    ax.loglog(hs, errors, 'D', label='Numerical error', color='darkblue', markersize=8)

    # Plot fitted line ONLY over the fitted region
    h_fit = hs[:fit_end]
    error_fit = np.exp(c) * h_fit**m

    ax.loglog(h_fit, error_fit, '--', label=f'Fitted line (m={m:.2f})', color='orange', linewidth=2)

    # Set tick sizes
    ax.tick_params(axis='both', which='both', labelsize=16)

    ax.set_xlabel(r"$h$", fontsize=20, labelpad=12)
    ax.set_ylabel(r"$|E_0^h - E_0| \, / \, \pi^2$",
              fontsize=20,
              labelpad=12)
    # Rotate y-axis label to horizontal
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.set_label_coords(-0.2, 0.5)
    ax.yaxis.label.set_horizontalalignment('right')
    ax.yaxis.label.set_verticalalignment('center')
    plt.legend(loc='center right', fontsize=20)
    plt.savefig(f"log_convergence_Crossgraph_noninteracting_fermions_eigval{n_eigval}.pdf", dpi=300, bbox_inches='tight', pad_inches=0.02)
    plt.show()

    return None

if __name__ =="__main__":

    alpha = 1.0

    log_convergence_YgraphHardcore_single_eigval(alpha, 0)

    pass