import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce
import math

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs
import ygraph as yga
import lasso as lsa
import dumbbell as db



if __name__ =="__main__":

    N = 150
    h = (np.pi)/(N-1)
    alpha = 0.0

    n = N-2

    curve = [0*n + 9, 1*n + 9 - 1 - 1, 2*n + 9 - 1 - 2 - 2, 3*n + 9 - 1 - 2 - 3 - 1 - 3, 3*n + 9 - 1 - 2 - 3 - 3, 4*n + 9 - 1 - 2 - 3 - 4 - 1 - 4, 5*n + 9 - 1 - 2 - 3 - 4 - 5 - 2 - 5, 6*n + 9 - 1 - 2 - 3 - 4 - 5 - 6 - 3 - 6, 6*n + 9 - 1 - 2 - 3 - 4 - 5 - 6 - 2 - 6]
    print(curve)

    # Angles
    thetas = []
    thetas.append(0)
    thetas.append(np.arctan(1/10))
    thetas.append(np.arctan(2/10))
    thetas.append(np.arctan(3/9.5))
    thetas.append(np.arctan(4/9))
    thetas.append(np.arctan(5/8))
    thetas.append(np.arctan(6/7))
    thetas.append(np.pi/4)
    print(thetas)

    CY = yga.YgraphAnyonsHardcore(N, alpha)
    states_filepath = "ygraph_hardcore_states/ygraph_N" + str(N) + "_alpha" + str(alpha) + "_"
    eigenvalues_filepath = "ygraph_hardcore_eigenvalues/ygraph_N" + str(N) + "_alpha" + str(alpha)
    CY.load_states(states_filepath)
    CY.load_eigenvalues(eigenvalues_filepath)

    CY.plot_states(0, plotting_method="surface", realimag="real", N_levels=20)

    D11 = CY.cells[0]

    state = D11.eigenstates[:,0]

    coords = D11.non_elim_coords
    coords = np.array(coords)
    print(coords)

    
    angular = state[curve]
    angular[3] = (angular[3]+ angular[4])/2
    angular[7] = (angular[7]+ angular[8])/2
    angular = np.delete(angular, [4,7], 0)
    angular = np.append(angular, 0)


    plt.plot(thetas, np.real(angular), 'o', label='Angular behaviour')
    # Label x axis with greek theta brackets radians
    plt.xlabel(r'$\theta$ (radians)')
    # Label y axis with greek psi_11
    plt.ylabel(r'$\psi_{11}$')
    # Give plot title - angular behaviour on D11
    plt.title('Angular behaviour on D11')
    plt.show()

    pass