import os

# Cap underlying library threads to 1
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Headless backend
import matplotlib.pyplot as plt

from mpi4py.futures import MPIPoolExecutor

# Custom graph modules
import ygraph as yga
import hgraph as hga
import lassosubdivided as lss
import dumbbellsubdivided as dbs
import cross as cga

# Plot styling
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'


def plot_energy_levels(spectra, alphas, filename, xlabel=r"$\alpha$"):
    fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)

    ys = [np.sort(spec) for spec in spectra]
    ys = np.array(ys)[:, 0:20]

    for i in range(ys.shape[1]):
        ax.plot(alphas, ys[:, i])

    ax.tick_params(axis='both', which='both', labelsize=14)
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(r"Energy / $\pi^2$", fontsize=18)
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.set_label_coords(-0.14, 0.5)

    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight', pad_inches=0.02)
        print(f"[TEST SUCCESS] Plot saved as {filename}")

    plt.close(fig)


# Worker Solver Functions
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

def solve_Cross_all_alphas_varying(args):
    N, alpha = args
    CC = cga.CrossgraphAnyonsHardcore(N, alpha, alpha, alpha)
    CC.gen_lapl()
    CC.lapl_solve((np.pi) / (N - 1), N_eigs=20)
    return CC.spectrum


if __name__ == "__main__":
    print("=== STARTING QUICK HPC TEST (N=20) ===")
    
    # Test settings
    N = 20
    alphas = np.linspace(-1, 1, 20)

    # Test running across distributed nodes via MPI
    with MPIPoolExecutor() as executor:
        print("Computing Y-Graph spectrum...")
        CYs = list(executor.map(solve_Ygraph_alpha, [(N, a) for a in alphas]))

        print("Computing H-Graph spectrum...")
        CHs_single = list(executor.map(solve_Hgraph_single_alpha_varying, [(N, a) for a in alphas]))

        print("Computing Cross-Graph spectrum...")
        CXs_all = list(executor.map(solve_Cross_all_alphas_varying, [(N, a) for a in alphas]))

    # Test plot generation
    plot_energy_levels(CYs, alphas, filename="TEST_Ygraph.pdf", xlabel=r"$\alpha$")
    plot_energy_levels(CHs_single, alphas, filename="TEST_Hgraph.pdf", xlabel=r"$\alpha_1$")
    plot_energy_levels(CXs_all, alphas, filename="TEST_Crossgraph.pdf", xlabel=r"$\alpha$")

    print("=== TEST FINISHED SUCCESSFULLY ===")