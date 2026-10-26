import os
import sys

# Threading controls
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mpi4py.futures import MPIPoolExecutor

import ygraph as yga
import hgraph as hga
import lassosubdivided as lss
import dumbbellsubdivided as dbs
import cross as cga

# Font formatting
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'


def plot_energy_levels(spectra, params, filename, xlabel=r"$\alpha$"):
    fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)

    ys = [np.sort(spec) for spec in spectra]
    ys = np.array(ys)[:, 0:20]

    for i in range(ys.shape[1]):
        ax.plot(params, ys[:, i])

    ax.tick_params(axis='both', which='both', labelsize=14)
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(r"Energy / $\pi^2$", fontsize=18)
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.set_label_coords(-0.14, 0.5)

    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight', pad_inches=0.02)
        print(f"[SUCCESS] Plot saved as {filename}")

    plt.close(fig)

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



# ==========================================
# Exact Worker Functions
# ==========================================

def solve_Ygraph_alpha(args):
    N, alpha = args
    CY = yga.YgraphAnyonsHardcore(N, alpha)
    CY.gen_lapl()
    CY.lapl_solve((np.pi) / (N - 1), N_eigs=20)
    return CY.spectrum

def solve_Hgraph_single_alpha_varying(args):
    N, alpha = args
    CH = hga.HgraphAnyonsHardcore(N, alpha, 0.0)
    CH.gen_lapl()
    CH.lapl_solve((np.pi) / (N - 1), N_eigs=20)
    return CH.spectrum

def solve_Hgraph_all_alphas_varying(args):
    N, alpha = args
    CH = hga.HgraphAnyonsHardcore(N, alpha, alpha)
    CH.gen_lapl()
    CH.lapl_solve((np.pi) / (N - 1), N_eigs=20)
    return CH.spectrum

def solve_Lassosubdivided_alpha_varying(args):
    N, alpha = args
    CL = lss.LassoSubdividedAnyonsHardcore(N, N, alpha, 0.0)
    CL.gen_lapl()
    CL.lapl_solve((np.pi) / (N - 1), N_eigs=20)
    return CL.spectrum

def solve_Lassosubdivided_beta_varying(args):
    N, beta = args
    CL = lss.LassoSubdividedAnyonsHardcore(N, N, 0.0, beta)
    CL.gen_lapl()
    CL.lapl_solve((np.pi) / (N - 1), N_eigs=20)
    return CL.spectrum

def solve_Dumbbellsubdivided_single_alpha_varying(args):
    N, alpha = args
    CD = dbs.DumbbellSubdividedAnyonsHardcore(N, N, N, alpha, 0.0, 0.0, 0.0)
    CD.gen_lapl()
    CD.lapl_solve((np.pi) / (N - 1), N_eigs=20)
    return CD.spectrum

def solve_Dumbbellsubdivided_all_alphas_varying(args):
    N, alpha = args
    CD = dbs.DumbbellSubdividedAnyonsHardcore(N, N, N, alpha, alpha, 0.0, 0.0)
    CD.gen_lapl()
    CD.lapl_solve((np.pi) / (N - 1), N_eigs=20)
    return CD.spectrum

def solve_Dumbbellsubdivided_single_beta_varying(args):
    N, beta = args
    CD = dbs.DumbbellSubdividedAnyonsHardcore(N, N, N, 0.0, 0.0, beta, 0.0)
    CD.gen_lapl()
    CD.lapl_solve((np.pi) / (N - 1), N_eigs=20)
    return CD.spectrum

def solve_Dumbbellsubdivided_all_betas_varying(args):
    N, beta = args
    CD = dbs.DumbbellSubdividedAnyonsHardcore(N, N, N, 0.0, 0.0, beta, beta)
    CD.gen_lapl()
    CD.lapl_solve((np.pi) / (N - 1), N_eigs=20)
    return CD.spectrum

def solve_Cross_single_alpha_varying(args):
    N, alpha = args
    CC = cga.CrossgraphAnyonsHardcore(N, alpha, 0.0, 0.0)
    CC.gen_lapl()
    # Request 30 eigenvalues instead of 20
    CC.lapl_solve((np.pi) / (N - 1), N_eigs=30)

    # Sort low-to-high so each index consistently tracks the same energy level
    sorted_spec = np.sort(CC.spectrum)
    return sorted_spec[:20]  # Slice the clean lowest 20

def solve_Cross_all_alphas_varying(args):
    N, alpha = args
    CC = cga.CrossgraphAnyonsHardcore(N, alpha, alpha, alpha)
    CC.gen_lapl()
    # Request 30 eigenvalues instead of 20
    CC.lapl_solve((np.pi) / (N - 1), N_eigs=30)

    # Sort low-to-high so each index consistently tracks the same energy level
    sorted_spec = np.sort(CC.spectrum)
    return sorted_spec[:20]  # Slice the clean lowest 20


# ==========================================
# Main Execution Dispatcher
# ==========================================

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python run_single_graph.py <task_id_1_to_11>")
        sys.exit(1)

    task_id = int(sys.argv[1])
    N = 100
    alphas = np.linspace(-1, 1, 49)
    betas = np.linspace(-1, 1, 49)

    # Task mapping dictionary
    TASKS = {
        1: (solve_Ygraph_alpha, alphas, "Energylevels_Ygraph", r"$\alpha$"),
        2: (solve_Hgraph_single_alpha_varying, alphas, "Energylevels_Hgraph_single_alpha", r"$\alpha_1$"),
        3: (solve_Hgraph_all_alphas_varying, alphas, "Energylevels_Hgraph_all_alphas", r"$\alpha$"),
        4: (solve_Lassosubdivided_alpha_varying, alphas, "Energylevels_Lassosubdivided_alpha", r"$\alpha$"),
        5: (solve_Lassosubdivided_beta_varying, betas, "Energylevels_Lassosubdivided_beta", r"$\beta$"),
        6: (solve_Dumbbellsubdivided_single_alpha_varying, alphas, "Energylevels_Dumbbellsubdivided_single_alpha", r"$\alpha_1$"),
        7: (solve_Dumbbellsubdivided_all_alphas_varying, alphas, "Energylevels_Dumbbellsubdivided_all_alphas", r"$\alpha$"),
        8: (solve_Dumbbellsubdivided_single_beta_varying, betas, "Energylevels_Dumbbellsubdivided_single_beta", r"$\beta_1$"),
        9: (solve_Dumbbellsubdivided_all_betas_varying, betas, "Energylevels_Dumbbellsubdivided_all_betas", r"$\beta$"),
        10: (solve_Cross_single_alpha_varying, alphas, "Energylevels_Cross_single_alpha", r"$\alpha_1$"),
        11: (solve_Cross_all_alphas_varying, alphas, "Energylevels_Cross_all_alphas", r"$\alpha$"),
    }

    if task_id not in TASKS:
        raise ValueError(f"Invalid task_id: {task_id}. Must be between 1 and 11.")

    func, params, name, xlabel = TASKS[task_id]

    print(f"=== STARTING TASK {task_id}/11: {name} (N={N}) ===")

    with MPIPoolExecutor() as executor:
        spectra = list(executor.map(func, [(N, p) for p in params]))

    # Save raw array data
    np.save(f"{name}_data.npy", spectra)
    print(f"[SUCCESS] Spectra data saved as {name}_data.npy")

    # Generate PDF plot
    plot_energy_levels(spectra, params, filename=f"{name}.pdf", xlabel=xlabel)

    print(f"=== TASK {task_id} COMPLETED SUCCESSFULLY ===")