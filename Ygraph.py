"""
Code to solve the Laplacian on the distinguishable configuration space of two particles on the Y graph
"""

import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import squareFDM
import twocomplex

def twoparticlesYgraph(N):
        # Evaluate the spectrum of the Laplacian on the distinguishable configuration space of two particles on the Y graph
        # Scaling factor
        h = (np.pi)/(N-1)
        
        global N_eigs
    
        D11 = squareFDM.Domain(N)
        D22 = squareFDM.Domain(N)
        D33 = squareFDM.Domain(N)
        D12 = squareFDM.Domain(N)
        D13 = squareFDM.Domain(N)
        D23 = squareFDM.Domain(N)
        D21 = squareFDM.Domain(N)
        D31 = squareFDM.Domain(N)
        D32 = squareFDM.Domain(N)
        D11.split_domain("dirichlet")
        D22.split_domain("dirichlet")
        D33.split_domain("dirichlet")
        
        cells = [D11,D22,D33,D12,D13,D23,D21,D31,D32]
        
        Network = twocomplex.Complex(cells)
        
        gluingA = [(D11, D11.x0),(D13, D13.x0),(D12, D12.x0)]
        gluingB = [(D22, D22.y0),(D12, D12.y0),(D32, D32.y0)]
        gluingC = [(D33, D33.x0),(D32, D32.x0),(D31, D31.x0)]
        gluingD = [(D11, D11.y0),(D31, D31.y0),(D21, D21.y0)]
        gluingE = [(D22, D22.x0),(D21, D21.x0),(D23, D23.x0)]
        gluingF = [(D33, D33.y0),(D23, D23.y0),(D13, D13.y0)]
        
        Network.diagonal_bc("dirichlet")
        
        Network.glue(gluingA)
        Network.glue(gluingB)
        Network.glue(gluingC)
        Network.glue(gluingD)
        Network.glue(gluingE)
        Network.glue(gluingF)
        
        Network.exterior_bc("dirichlet")
        
        Network.simplify_lapl()
        spectrum = Network.lapl_spectrum(h,2,N_eigs)
        
        return spectrum

if __name__ == "__main__":
        # Main
        
        #N = 80 # The length of a wire
        
        Nmax = 130
        
        global N_eigs
        N_eigs = 30
        
        eval = []
        
        Ns = [x for x in range (10,Nmax,10)]
        
        #for N in Ns:
        #        spectrum = twoparticlesYgraph(N,N_eigs)
        #        eval.append(2*spectrum)
        #        pass
        
        eval = list(map(twoparticlesYgraph,Ns))
        
        print(eval)
        
        eval = 2*np.array(eval)

        for i in range (0,N_eigs):
                plt.plot(Ns,eval[:,i])
                pass
        plt.show()
        pass
"""
    #print(spectrum)
    #print("Scaling factor:")
    #print(2*spectrum)
    #Network.print_eqs()

    spectrum, states = Network.lapl_solve(h,2)
    print(spectrum)
    n = 4
    for cell in cells:
            plot = twocomplex.plot_state(cell, states[:,n])
            plot.show()
            pass
"""