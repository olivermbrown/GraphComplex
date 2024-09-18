import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce
import os
import sys
import multiprocessing as mp

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs
import ygraph as yg
import lasso as lsa
import dumbbell as db

def save_ygraph_eigenstates(N, alpha, n_eigs):
    """
    Save the eigenstates of the ConfigurationSpace object for the Y graph to a file.

    Args:
        N (int): The number of points in the discretization.
        alpha (float): The parameter in the boundary conditions.
        n_eigs (int): The number of eigenstates to save.

    Returns:
        None.
    """

    #N, alpha, n_eigs = args_tuple
    
    h = (np.pi)/(N-1)

    folder = "ygraph_states/"
    filename = "ygraph_N"+str(N)+"_alpha"+str(alpha)
    file_path = folder+filename+"_"

    # Create folder if it doesn't exist
    if not os.path.exists(folder):
        try:
            os.makedirs(folder)
            pass
        except FileExistsError:
            pass
        pass

    eigs_folder = "ygraph_eigenvalues/"
    if not os.path.exists(eigs_folder):
        try:
            os.makedirs(eigs_folder)
            pass
        except FileExistsError:
            pass
        pass
    eigs_path = eigs_folder+filename

    print("alpha = ", alpha)

    # Create the ConfigurationSpace object.
    CY = yg.YgraphAnyons(N, alpha)

    # Generate the Laplacian matrix.
    CY.gen_lapl()

    # Solve the eigenvalue problem.
    CY.lapl_solve(h, N_eigs=n_eigs)

    # Save the values to a file
    CY.save_eigenvalues(eigs_path)
    print("Eigenvalues saved to ", eigs_path)

    # Save the eigenstates to the file.
    CY.save_states(file_path)
    print("Eigenstates saved to ", file_path)

    return None

def save_lasso_eigenstates(N, alpha, n_eigs):
    """
    Save the eigenstates of the ConfigurationSpace object for the lasso graph to a file.

    Args:
        N (int): The number of points in the discretization.
        alpha (float): The parameter in the boundary conditions.
        n_eigs (int): The number of eigenstates to save.

    Returns:
        None.
    """

    #N, alpha, n_eigs = args_tuple
    
    h = (np.pi)/(N-1)

    folder = "lasso_states/"
    filename = "lasso_N"+str(N)+"_alpha"+str(alpha)
    file_path = folder+filename+"_"

    # Create folder if it doesn't exist
    if not os.path.exists(folder):#
        try:
            os.makedirs(folder)
            pass
        except FileExistsError:
            pass
        pass

    eigs_folder = "lasso_eigenvalues/"
    if not os.path.exists(eigs_folder):
        try:
            os.makedirs(eigs_folder)
            pass
        except FileExistsError:
            pass
        pass
    eigs_path = eigs_folder+filename

    print("alpha = ", alpha)

    # Create the ConfigurationSpace object.
    CY = lsa.LassoAnyons(N, alpha)

    # Generate the Laplacian matrix.
    CY.gen_lapl()

    # Solve the eigenvalue problem.
    CY.lapl_solve(h, N_eigs=n_eigs)

    # Save the values to a file
    CY.save_eigenvalues(eigs_path)
    print("Eigenvalues saved to ", eigs_path)

    # Save the eigenstates to the file.
    CY.save_states(file_path)
    print("Eigenstates saved to ", file_path)

    return None

def save_dumbbell_eigenstates(N, alpha, n_eigs):
    """
    Save the eigenstates of the ConfigurationSpace object for the lasso graph to a file.

    Args:
        N (int): The number of points in the discretization.
        alpha (float): The parameter in the boundary conditions.
        n_eigs (int): The number of eigenstates to save.

    Returns:
        None.
    """

    #N, alpha, n_eigs = args_tuple
    
    h = (np.pi)/(N-1)

    folder = "dumbbell_states/"
    filename = "dumbbell_N"+str(N)+"_alpha"+str(alpha)
    file_path = folder+filename+"_"

    # Create folder if it doesn't exist
    if not os.path.exists(folder):
        os.makedirs(folder)
        pass

    eigs_folder = "dumbbell_eigenvalues/"
    if not os.path.exists(eigs_folder):
        os.makedirs(eigs_folder)
        pass
    eigs_path = eigs_folder+filename

    # Create the ConfigurationSpace object.
    C = db.DumbbellAnyons(N, alpha)

    # Generate the Laplacian matrix.
    C.gen_lapl()

    # Solve the eigenvalue problem.
    C.lapl_solve(h, N_eigs=n_eigs)

    # Save the values to a file
    C.save_eigenvalues(eigs_path)

    # Save the eigenstates to the file.
    C.save_states(file_path)

    return None

if __name__ == "__main__":
    # Main
    alphas = np.linspace(0, 1, 25)
    Ns = np.linspace(110, 110, 5)
    args = [(int(N), float(alpha), 2) for N in Ns for alpha in alphas]

    # i defined as system argument
    i = int(sys.argv[1])

    # Take ith entry of args
    arg = args[i]

    # Set up the pool of workers
    # with mp.Pool() as pool:
    #     # Parallel execution using starmap which directly supports multiple arguments
    #     pool.starmap(save_lasso_eigenstates, args)

    save_ygraph_eigenstates(*arg)

    pass