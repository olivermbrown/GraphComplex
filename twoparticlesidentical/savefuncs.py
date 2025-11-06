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
import twowires as tw
import lassosubdivided as lss
import cross as crs
import star as star

def save_ygraph_eigenstates_hardcore(N, alpha, n_eigs):
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

    folder = "ygraph_hardcore_states/"
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

    eigs_folder = "ygraph_hardcore_eigenvalues/"
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
    CY = yg.YgraphAnyonsHardcore(N, alpha)

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

def save_ygraph_zoom_eigenstates_hardcore(N, alpha, n_eigs):
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

    folder = "ygraph_zoom_hardcore_states/"
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

    eigs_folder = "ygraph_zoom_hardcore_eigenvalues/"
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
    CY = yg.YgraphAnyonsHardcore(N, alpha)

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

def save_lasso_eigenstates_hardcore(N, alpha, n_eigs):
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

    folder = "lasso_hardcore_states/"
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

    eigs_folder = "lasso_hardcore_eigenvalues/"
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
    CY = lsa.LassoAnyonsHardcore(N, alpha)

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

def save_dumbbell_eigenstates_hardcore(N, alpha1, alpha2, n_eigs):
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

    folder = "dumbbell_hardcore_states/"
    filename = "dumbbell_N"+str(N)+"_alpha1"+str(alpha1)+"_alpha2"+str(alpha2)
    file_path = folder+filename+"_"

    # Create folder if it doesn't exist
    if not os.path.exists(folder):
        os.makedirs(folder)
        pass

    eigs_folder = "dumbbell_hardcore_eigenvalues/"
    if not os.path.exists(eigs_folder):
        os.makedirs(eigs_folder)
        pass
    eigs_path = eigs_folder+filename

    # Create the ConfigurationSpace object.
    C = db.DumbbellAnyonsHardcore(N, alpha1, alpha2)

    # Generate the Laplacian matrix.
    C.gen_lapl()

    # Solve the eigenvalue problem.
    C.lapl_solve(h, N_eigs=n_eigs)

    # Save the values to a file
    C.save_eigenvalues(eigs_path)

    # Save the eigenstates to the file.
    C.save_states(file_path)
    print("Eigenstates saved to ", file_path)

    return None

def save_ygraph_eigenstates_neumann(N, alpha, n_eigs):
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

    folder = "ygraph_neumann_states/"
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

    eigs_folder = "ygraph_neumann_eigenvalues/"
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
    CY = yg.YgraphAnyonsRobin(N, alpha)

    CY.robin_constant = 0

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

def save_ygraph_eigenstates_robin(N, alpha, n_eigs):
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

    folder = "ygraph_robin_states/"
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

    eigs_folder = "ygraph_robin_eigenvalues/"
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
    CY = yg.YgraphAnyonsRobin(N, alpha)

    CY.robin_constant = 1

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

def save_lasso_eigenstates_neumann(N, alpha, n_eigs):
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

    folder = "lasso_neumann_states/"
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

    eigs_folder = "lasso_neumann_eigenvalues/"
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
    CY = lsa.LassoAnyonsContact(N, alpha)

    CY.robin_constant = 0

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

def save_lasso_eigenstates_robin(N, alpha, n_eigs):
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

    folder = "lasso_robin_states/"
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

    eigs_folder = "lasso_robin_eigenvalues/"
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
    CY = lsa.LassoAnyonsContact(N, alpha)

    CY.robin_constant = 1

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

def save_dumbbell_eigenstates_neumann(N, alpha1, alpha2, n_eigs):
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

    folder = "dumbbell_neumann_states/"
    filename = "dumbbell_N"+str(N)+"_alpha1"+str(alpha1)+"_alpha2"+str(alpha2)
    file_path = folder+filename+"_"

    # Create folder if it doesn't exist
    if not os.path.exists(folder):
        os.makedirs(folder)
        pass

    eigs_folder = "dumbbell_neumann_eigenvalues/"
    if not os.path.exists(eigs_folder):
        os.makedirs(eigs_folder)
        pass
    eigs_path = eigs_folder+filename

    # Create the ConfigurationSpace object.
    C = db.DumbbellAnyonsNeumann(N, alpha1, alpha2)

    C.robin_constant = 0

    # Generate the Laplacian matrix.
    C.gen_lapl()

    # Solve the eigenvalue problem.
    C.lapl_solve(h, N_eigs=n_eigs)

    # Save the values to a file
    C.save_eigenvalues(eigs_path)

    # Save the eigenstates to the file.
    C.save_states(file_path)
    print("Eigenstates saved to ", file_path)

    return None

def save_twowires_eigenstates_hardcore(N, alpha1, alpha2, n_eigs):
    """
    Save the eigenstates of the ConfigurationSpace object for the two wires graph to a file.

    Args:
        N (int): The number of points in the discretization.
        alpha1 (float): The parameter in the boundary conditions for wire 1.
        alpha2 (float): The parameter in the boundary conditions for wire 2.
        n_eigs (int): The number of eigenstates to save.

    Returns:
        None.
    """

    #N, alpha1, alpha2, n_eigs = args_tuple
    
    h = (np.pi)/(N-1)

    folder = "twowires_hardcore_states/"
    filename = "twowires_N"+str(N)+"_alpha1"+str(alpha1)+"_alpha2"+str(alpha2)
    file_path = folder+filename+"_"

    # Create folder if it doesn't exist
    if not os.path.exists(folder):
        os.makedirs(folder)
        pass

    eigs_folder = "twowires_hardcore_eigenvalues/"
    if not os.path.exists(eigs_folder):
        os.makedirs(eigs_folder)
        pass
    eigs_path = eigs_folder+filename

    # Create the ConfigurationSpace object.
    C = tw.TwoWiresWithPhases(N, alpha1, alpha2)

    # Generate the Laplacian matrix.
    C.gen_lapl()

    # Solve the eigenvalue problem.
    C.lapl_solve(h, N_eigs=n_eigs)

    # Save the values to a file
    C.save_eigenvalues(eigs_path)

    # Save the eigenstates to the file.
    C.save_states(file_path)
    print("Eigenstates saved to ", file_path)

    return None

def save_lassosubdivided_eigenstates_hardcore(N, alpha, aharonovbohm, n_eigs):
    """
    Save the eigenstates of the ConfigurationSpace object for the lasso subdivided graph to a file.

    Args:
        N (int): The number of points in the discretization for the first wire.
        N2 (int): The number of points in the discretization for the second wire.
        alpha (float): The parameter in the boundary conditions.
        aharonovbohm (float): The Aharonov-Bohm phase.
        n_eigs (int): The number of eigenstates to save.

    Returns:
        None.
    """

    #N, N2, alpha, aharonovbohm, n_eigs = args_tuple
    
    h = (np.pi)/(N-1)

    folder = "lassosubdivided_hardcore_states/"
    filename = "lassosubdivided_N"+str(N)+"_alpha"+str(alpha)+"_aharonovbohm"+str(aharonovbohm)
    file_path = folder+filename+"_"

    # Create folder if it doesn't exist
    if not os.path.exists(folder):
        os.makedirs(folder)
        pass

    eigs_folder = "lassosubdivided_hardcore_eigenvalues/"
    if not os.path.exists(eigs_folder):
        os.makedirs(eigs_folder)
        pass
    eigs_path = eigs_folder+filename

    # Create the ConfigurationSpace object.
    C = lss.LassoSubdividedAnyonsHardcore(N, N, alpha, aharonovbohm)

    # Generate the Laplacian matrix.
    C.gen_lapl()

    # Solve the eigenvalue problem.
    C.lapl_solve(h, N_eigs=n_eigs)

    # Save the values to a file
    C.save_eigenvalues(eigs_path)

    # Save the eigenstates to the file.
    C.save_states(file_path)
    print("Eigenstates saved to ", file_path)

    return None

def save_cross_eigenstates_hardcore(N, alpha, n_eigs):
    """
    Save the eigenstates of the ConfigurationSpace object for the cross graph to a file.

    Args:
        N (int): The number of points in the discretization.
        alpha (float): The parameter in the boundary conditions.
        n_eigs (int): The number of eigenstates to save.

    Returns:
        None.
    """

    #N, alpha, n_eigs = args_tuple
    
    h = (np.pi)/(N-1)

    folder = "cross_hardcore_states/"
    filename = "cross_N"+str(N)+"_alpha"+str(alpha)
    file_path = folder+filename+"_"

    # Create folder if it doesn't exist
    if not os.path.exists(folder):
        try:
            os.makedirs(folder)
            pass
        except FileExistsError:
            pass
        pass

    eigs_folder = "cross_hardcore_eigenvalues/"
    if not os.path.exists(eigs_folder):
        try:
            os.makedirs(eigs_folder)
            pass
        except FileExistsError:
            pass
        pass
    eigs_path = eigs_folder+filename

    print("alpha = ", alpha)

    alpha2 = 0.0
    alpha3 = 0.0

    # Create the ConfigurationSpace object.
    CY = crs.CrossgraphAnyonsHardcore(N, alpha, alpha2, alpha3)

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

def save_star_eigenstates_hardcore(N, alpha, n_eigs):
    """
    Save the eigenstates of the ConfigurationSpace object for the star graph to a file.

    Args:
        N (int): The number of points in the discretization.
        alpha (float): The parameter in the boundary conditions.
        n_eigs (int): The number of eigenstates to save.

    Returns:
        None.
    """

    #N, alpha, n_eigs = args_tuple
    
    h = (np.pi)/(N-1)

    folder = "star_hardcore_states/"
    filename = "star_N"+str(N)+"_alpha"+str(alpha)
    file_path = folder+filename+"_"

    # Create folder if it doesn't exist
    if not os.path.exists(folder):
        try:
            os.makedirs(folder)
            pass
        except FileExistsError:
            pass
        pass

    eigs_folder = "star_hardcore_eigenvalues/"
    if not os.path.exists(eigs_folder):
        try:
            os.makedirs(eigs_folder)
            pass
        except FileExistsError:
            pass
        pass
    eigs_path = eigs_folder+filename

    print("alpha = ", alpha)

    alpha2 = 0.0
    alpha3 = 0.0

    # Create the ConfigurationSpace object.
    CY = star.StargraphBosonsHardcore(N)

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

if __name__ == "__main__":
    # Main
    alphas = [0.0]#np.linspace(0, 1, 25)
    N = 150
    args = [(int(N), float(alpha), 50) for alpha in alphas]

    # i defined as system argument
    i = int(sys.argv[1])

    # Take ith entry of args
    arg = args[i]

    # Set up the pool of workers
    # with mp.Pool() as pool:
    #     # Parallel execution using starmap which directly supports multiple arguments
    #     pool.starmap(save_lasso_eigenstates, args)

    save_star_eigenstates_hardcore(*arg)

    pass