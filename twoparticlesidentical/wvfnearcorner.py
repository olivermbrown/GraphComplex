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

# def compute_arc_integrals(xs, ys, zs, dx, rs, theta_min=0, theta_max=np.pi/2):
#     # Precompute R and Theta for all points
#     R = np.sqrt(xs**2 + ys**2)
#     Theta = np.mod(np.arctan2(ys, xs), 2*np.pi)

#     # Apply static angle mask once
#     angle_mask = (Theta >= theta_min) & (Theta <= theta_max)

#     # Masked R, zs
#     R_masked = R[angle_mask]
#     zs_masked = zs[angle_mask]

#     # Sort masked R and zs by R once
#     sort_idx = np.argsort(R_masked)
#     R_sorted = R_masked[sort_idx]
#     zs_sorted = zs_masked[sort_idx]
#     psi2sorted = np.abs(zs_sorted)**2

#     # Use searchsorted to find cutoff for each r
#     cutoff_indices = np.searchsorted(R_sorted, rs)

#     # Cumulative sum of ψ values
#     zs_cumsum = np.cumsum(psi2sorted)

#     # Extract cumulative sums at each radius cutoff
#     psi_sums = zs_cumsum[cutoff_indices - 1]  # Note: rs[0] must be > 0

#     # Final integral values
#     psi_integrals = np.abs(psi_sums) * dx**2

#     return psi_integrals

# def compute_arc_integrals(xs, ys, zs, dx, rs, theta_min=0, theta_max=np.pi/2):
#     # Polar coordinates
#     R = np.sqrt(xs**2 + ys**2)
#     Theta = np.mod(np.arctan2(ys, xs), 2 * np.pi)

#     # Angle mask
#     angle_mask = (Theta >= theta_min) & (Theta <= theta_max)
#     R_masked = R[angle_mask]
#     zs_masked = zs[angle_mask]

#     # Sort by radius
#     sort_idx = np.argsort(R_masked)
#     R_sorted = R_masked[sort_idx]
#     zs_sorted = zs_masked[sort_idx]
#     psi2_sorted = np.abs(zs_sorted)**2

#     # Cumulative sum
#     cumsum = np.cumsum(psi2_sorted)

#     # Index of cutoff point for each radius
#     cutoff_indices = np.searchsorted(R_sorted, rs)

#     # Initialize results
#     psi_integrals = np.zeros_like(rs, dtype=np.float64)

#     for i, cutoff in enumerate(cutoff_indices):
#         if cutoff == 0:
#             psi_integrals[i] = 0.0
#         else:
#             # Sum up to cutoff-1 fully
#             full_sum = cumsum[cutoff - 1] if cutoff > 1 else psi2_sorted[0]
#             # Add half the last ring if it exists
#             edge_correction = 0.5 * psi2_sorted[cutoff] if cutoff < len(psi2_sorted) else 0.0
#             psi_integrals[i] = (full_sum + edge_correction) * dx**2

#     return psi_integrals

# def cell_weight(x, y, dx, r, samples_per_side=5):
#     """Approximate fractional area of square centered at (x, y) that lies inside circle of radius r."""
#     half = dx / 2
#     xs = np.linspace(x - half, x + half, samples_per_side)
#     ys = np.linspace(y - half, y + half, samples_per_side)
#     Xs, Ys = np.meshgrid(xs, ys, indexing='ij')
#     Rs = np.sqrt(Xs**2 + Ys**2)
#     return np.mean(Rs <= r)

# def compute_arc_integrals_weighted(xs, ys, zs, dx, rs, theta_min=0, theta_max=np.pi/2, samples_per_side=5):
#     # Polar coordinates
#     R = np.sqrt(xs**2 + ys**2)
#     Theta = np.mod(np.arctan2(ys, xs), 2*np.pi)

#     # Angle mask
#     angle_mask = (Theta >= theta_min) & (Theta <= theta_max)
#     xs_masked = xs[angle_mask]
#     ys_masked = ys[angle_mask]
#     zs_masked = zs[angle_mask]

#     psi2_masked = np.abs(zs_masked)**2
#     n_pts = len(psi2_masked)

#     # Initialize output
#     psi_integrals = np.zeros_like(rs, dtype=np.float64)

#     # For each radius, weight each point based on partial cell area
#     for i, r in enumerate(rs):
#         weights = np.array([
#             cell_weight(xs_masked[j], ys_masked[j], dx, r, samples_per_side)
#             for j in range(n_pts)
#         ])
#         weighted_sum = np.sum(psi2_masked * weights)
#         psi_integrals[i] = weighted_sum * dx**2

#     return psi_integrals

def compute_arc_integrals_fast(xs, ys, zs, dx, rs, theta_min=0, theta_max=np.pi/2):
    # Compute polar coordinates
    R = np.sqrt(xs**2 + ys**2)
    Theta = np.mod(np.arctan2(ys, xs), 2 * np.pi)

    # Angle mask: static
    angle_mask = (Theta >= theta_min) & (Theta <= theta_max)
    R_masked = R[angle_mask]
    psi2_masked = np.abs(zs[angle_mask])**2

    safe_r = np.maximum(rs, 1e-10)  # avoids divide-by-zero
    partial_weights = (rs**2 - (rs - dx/2)**2) / (dx * safe_r)

    # Prepare output
    psi_integrals = np.zeros_like(rs, dtype=np.float64)

    for i, r in enumerate(rs):
        # Full interior: R < r - dx/2
        full_mask = R_masked < (r - dx / 2)
        partial_mask = ((R_masked >= (r - dx / 2)) & (R_masked < (r + dx / 2)))

        full_sum = np.sum(psi2_masked[full_mask])
        partial_sum = np.sum(psi2_masked[partial_mask])

        # Weight partial cells by partial weight
        psi_integrals[i] = dx**2 * (full_sum + partial_weights[i] * partial_sum)

        pass


    return psi_integrals


def plot_corner_behaviour(complex, cells, show_plots=True, lasso=False, state_number=None, extrapolate_to_origin=False, start=3, cutoff=25):
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

    max_node = cutoff
    min_node = 0
    psi_integrals = np.zeros(shape=(len(cells), max_node - min_node + 1), dtype="complex128")
    
    for cell in cells:
        # Import eigenstate data from cell
        i = complex.cells.index(cell)
        cell_data = data[i]

        xs, ys, zs = cell_data

        if lasso:
            if cell == complex.cells[0]:
                xs = 1 - xs
                ys = 1 - ys
                pass
            elif cell == complex.cells[3]:
                xs = 1 - xs
                pass
            else:
                pass
            pass
        else:
            pass

        N = cell.N - 2
        dx = 1/(cell.N - 1)

        # Calculate range of radii
        rs = np.arange(min_node, max_node+1, 1) * dx

        
        if type(cell) == cls.SquareCell or type(cell) == cls.TriangleCell:
            # If the cell is a square cell, we can proceed

            
            # Find the indices of the nodes near to the corner of the square cell
            if N % 2 == 0:
                # If the cell has an even number of nodes, the hypotenuse is the diagonal from (0,0) to (N-1,N-1)

                psi_integrals[cells.index(cell)] = compute_arc_integrals_fast(xs, ys, zs, dx, rs, theta_min=0, theta_max=np.pi/2)
                
                pass
            elif N % 2 == 1:
                raise Exception("This function is only for square cells with an even number of nodes")
            else:
                raise Exception("Unknown cell type")
            pass
        else:
            raise Exception("Unknown cell type")
        pass
    else:
        pass

    psi_integrals = np.sum(psi_integrals, axis=0)

    # Remove first few integrals
    # rem = psi_integrals[start-1]
    # psi_integrals = psi_integrals[start-1:]
    # rs = rs[start-1:]
    # psi_integrals -= rem

    # Plot the log log plot of the integrated wavefunction
    plt.clf()
    if show_plots:
        plt.scatter(np.log(rs),np.log(np.abs(psi_integrals)), label='Integrated wavefunction', marker='x', color='blue')
        #plt.scatter(np.log(rs), -1.5*rs, label='Radii', color='orange')
        plt.xlabel(r"$\log(R)$", fontsize=14)
        plt.ylabel(r"$\log(I(R))$", fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        pass
    else:
        plt.close('all')

    #DlogPsi = np.log(np.abs(psi[1])) - np.log(np.abs(psi[0]))
    #DlogR = np.log(np.abs(ys[1])) - np.log(np.abs(ys[0]))
    #grad = DlogPsi / DlogR
    #print("Gradient of the wavefunction at the origin: ", grad)

    # Calculate the slope of the line from first several points

    # Fit a polynomial to the data
    degree = 0
    polyfeatures = PolynomialFeatures(degree).fit_transform(rs[start:].reshape(-1,1))
    logfeatures = np.log(np.abs(rs[start:]).reshape(-1,1))
    features = np.concatenate((logfeatures,polyfeatures),axis=1)
    #print(features.shape)
    reg = LinearRegression(fit_intercept=False).fit(features, np.log(np.abs(psi_integrals[start:])))
    m = reg.coef_
    #print(reg.predict(features))
    #print(np.log(np.abs(psi[:10])))
    r = np.linspace(rs[start], rs[-1], 100)
    plottingfeatures = PolynomialFeatures(degree).fit_transform(r.reshape(-1,1))
    plottinglogfeatures = np.log(np.abs(r).reshape(-1,1))
    plotting = np.concatenate((plottinglogfeatures,plottingfeatures),axis=1)
    #polynomial = np.dot(plotting, m)
    polynomial = reg.predict(plotting)

    # Calculate the gradient of the polynomial through xs[0]
    derivative = np.polyder(m)
    print(m)
    print("beta = ", m[0]/2 - 1)
    # Print blank line
    print()

    # Plot show the log log plot with the linear fit through the first three points in red
    if show_plots:
        plt.plot(np.log(r), polynomial, color='darkorange', label='Linear fit')

        plt.legend(fontsize=14, loc='upper left')

        plt.show()
        pass
    else:
        plt.close('all')
        pass

    return m[0]/2 - 1

def ygraph_hardcore_analytic():

    betas = [0.783653, 0.784082, 0.785365, 0.7875, 0.790479, 0.794293, 0.798929, 0.804372, 0.810607, 0.817615, 0.825376, 0.833868, 0.843069, 0.852956, 0.863504, 0.874689, 0.886487, 0.898872, 0.911819, 0.925306, 0.939307, 0.9538, 0.968761, 0.984168, 1]

    return betas

def lassosubdivided_hardcore_analytic():

    betas = [0.783653, 0.784201, 0.78584, 0.78855, 0.792301, 0.797053, 0.802755, 0.80935, 0.816776, 0.824969, 0.833863, 0.843396, 0.853512, 0.864158, 0.875295, 0.886889, 0.898922, 0.911387, 0.924293, 0.937658, 0.951504, 0.965804, 0.980287, 0.993531, 1.]

    return betas

def ygraph_hardcore_interpolation_corner(N):

    alphas = np.linspace(0, 1, 25)
    h = (np.pi)/(N-1)

    CYs = []
    ms = []
    eus = []
    els = []
    for alpha in alphas:
        C = yga.YgraphAnyonsHardcore(N, alpha)
        states_filepath = "ygraph_hardcore_states/ygraph_N" + str(N) + "_alpha" + str(alpha) + "_"
        eigenvalues_filepath = "ygraph_hardcore_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
        C.load_states(states_filepath)
        C.load_eigenvalues(eigenvalues_filepath)
        #C.gen_lapl()
        #C.lapl_solve(h, N_eigs=2)
        CYs.append(C)
        cells = C.cells
        print("alpha = ", alpha)
        yga.ArrangeYgraphPlots(C)
        m = plot_corner_behaviour(C, cells, show_plots=True, cutoff=20)
        eu = plot_corner_behaviour(C, cells, show_plots=False, cutoff=25)
        el = plot_corner_behaviour(C, cells, show_plots=False, cutoff=15)
        ms.append(m)
        eus.append(eu)
        els.append(el)
        pass

    # Calculate error bars
    eus = np.array(eus)
    els = np.array(els)
    ers = np.maximum(eus, els)

    # Produce scatter plot of m against alpha
    plt.clf()
    betas_analytic = ygraph_hardcore_analytic()
    plt.errorbar(alphas, ms, yerr=ers, fmt='x', capsize=5, label='Data with error bars')
    plt.plot(alphas, betas_analytic, label='Analytic estimate', color='orange')
    #plt.scatter(alphas, ms, label='Numerical', color='blue')
    # X axis label
    plt.xlabel(r"$\alpha$", fontsize=14)
    # Y axis label
    plt.ylabel(r"$\beta$", rotation=0, fontsize=14)
    # Rotate y axis label
    plt.yticks(rotation=0, fontsize=14)
    plt.xticks(fontsize=14)
    # Add legend
    plt.legend(fontsize=14)
    plt.show()

    return None

def lassosubdivided_hardcore_interpolation_corner(N):

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
        cells = list(np.array(C.cells)[[0,1,3]])
        print("alpha = ", alpha)
        m = plot_corner_behaviour(C,cells,show_plots=False,lasso=True)
        ms.append(m)
        pass

    # Produce scatter plot of m against alpha
    plt.clf()
    #plt.plot(alphas, ms)
    #betas_analytic = ygraph_hardcore_analytic()
    #plt.errorbar(alphas, ms, yerr=ses, fmt='.', capsize=5, label='Data with error bars')
    betas_analytic = lassosubdivided_hardcore_analytic()
    plt.plot(alphas, betas_analytic, label='Analytic', color='orange')
    plt.scatter(alphas, ms, label='Numerical', color='blue')
    # X axis label
    plt.xlabel(r"$\alpha$")
    # Y axis label
    plt.ylabel(r"$\beta$", rotation=0)
    # Rotate y axis label
    plt.yticks(rotation=0)
    plt.show()

    return None



if __name__ =="__main__":

    N = 150
    h = (np.pi)/(N-1)
    

    ygraph_hardcore_interpolation_corner(N)
    #lassosubdivided_hardcore_interpolation_corner(N)

    pass