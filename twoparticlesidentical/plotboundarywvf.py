import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce
import math

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs
import ygraphanyons as yga
import lassoanyons as lsa

def plot_boundary_wavefunction(complex, edge):

    cell = edge.domain
    i = complex.cells.index(cell)

    data = complex.plot_states(0, return_data=True, show_plots=False)
    cell_data = data[i]

    xs, ys, zs = cell_data
    if type(cell) == cls.SquareCell:
        edge_length = (cell.N)-2
        pass
    elif type(cell) == cls.TriangleCell:
        edge_length = (cell.N)-3
        pass
    else:
        raise Exception("Unknown cell type")
        pass

    # Select the function values on the given edge
    if type(cell) == cls.SquareCell:
        if edge == cell.x0:
            psi = zs[:edge_length]
            pass
        elif edge == cell.x1:
            psi = zs[edge_length**2 - edge_length:edge_length**2]
            pass
        elif edge == cell.y0:
            psi = zs[0:edge_length**2:edge_length]
            pass
        elif edge == cell.y1:
            psi = zs[edge_length-1:edge_length**2:edge_length]
            pass
        else:
            raise Exception("Unknown edge")
            pass
        pass
    elif type(cell) == cls.TriangleCell:
        if edge == cell.x0:
            psi = zs[:edge_length]
            pass
        elif edge == cell.y1:
            ids = []
            c = 0
            for i in range(edge_length):
                ids.append(edge_length*(i+1) - 1 - c - i)
                c += i
                pass
            psi = zs[ids]
            pass
        else:
            raise Exception("Unknown edge")
            pass
        pass
    else:
        raise Exception("Unknown cell type")
        pass

    # Plot the wavefunction on the edge
    plt.clf()
    plt.plot(ys[:edge_length], np.abs(psi))
    plt.show()

    # Plot the log log plot of the wavefunction on the edge
    plt.clf()
    plt.scatter(np.log(ys[:edge_length]),np.log(np.abs(psi)))

    # Calculate the slope of the line from first three points
    m = np.polyfit(np.log(ys[:3]), np.log(np.abs(psi[:3])), 1)
    r = np.linspace(np.log(ys[0]), np.log(ys[2]), 100)
    polynomial = np.polyval(m, r)

    # Calculate the gradient of the polynomial through xs[0]
    derivative = np.polyder(m)
    print(np.polyval(derivative, np.log(ys[0])))

    # Plot show the log log plot with the linear fit through the first three points in red
    plt.plot(r, polynomial, color='red')
    plt.show()

    return None

if __name__ =="__main__":

    N = 100
    h = (np.pi)/(N-1)
    alpha = 0

    C = lsa.LassoAnyons(N, alpha)
    C.lapl_solve(h,2,20)

    states = C.states
    gs = states[:,0]

    D12 = C.cells[1]
    edge = D12.x1

    plot_boundary_wavefunction(C, edge)

    pass