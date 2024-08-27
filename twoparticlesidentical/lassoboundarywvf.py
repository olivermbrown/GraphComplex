import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce
import math

import matplotlib.pyplot as plt
from matplotlib import cm

import cells as cls
import configs
import lassoanyons as lsa

def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

if __name__=="__main__":

    N = 60
    h = (np.pi)/(N-1)
    alpha = 1

    CY = lsa.LassoAnyons(N, alpha)
    CY.lapl_solve(h,2,20)

    states = CY.states
    gs = states[:,0]

    D21 = CY.cells[1]
    edge = D21.y0

    data = CY.plot_states(0, return_data=True, show_plots=True)
    D21_data = data[1]
    xs, ys, zs = D21_data
    edge_length = N-2
    plt.clf()
    plt.plot(xs[:edge_length], np.abs(zs[:edge_length]))
    plt.show()
    plt.clf()
    plt.scatter(np.log(xs[:edge_length]),np.log(np.abs(zs[:edge_length])))
    # Calculate the slope of the line from first three points
    nearest = find_nearest(np.log(xs), -2)
    cutoff = list(np.log(xs)).index(nearest)
    m = np.polyfit(np.log(xs[:cutoff]), np.log(np.abs(zs[:cutoff])), 5)
    r = np.linspace(np.log(xs[0]), np.log(xs[cutoff-1]), 100)
    polynomial = np.polyval(m, r)
    print(m)
    # Calculate the gradient of the polynomial through xs[0]
    derivative = np.polyder(m)
    print(np.polyval(derivative, np.log(xs[0])))
    plt.plot(r, polynomial, color='red')
    plt.show()
    pass