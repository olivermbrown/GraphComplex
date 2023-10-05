"""
Glue together pieces to form a large triangle
"""

import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import squareFDM
import twocomplex

if __name__ == "__main__":
    # Main
    
    N = 50 # The length of a wire
    
    # Scaling factor
    h = (np.pi)/(N-1)
    
    D1 = squareFDM.Domain(N)
    D2 = squareFDM.Domain(N)
    D2.split_domain()
    
    cells = [D1,D2]
    
    Network = twocomplex.Complex(cells)
    
    gluingA = [(D1, D1.x0),(D2, D2.x1)]
    gluingB = [(D2, D2.y0),(D1, D1.y1)]
    
    
    Network.glue(gluingA)
    Network.glue(gluingB)
    
    Network.exterior_bc("dirichlet")
    Network.diagonal_bc("dirichlet")
    
    N = Network.simplify_lapl()
    new = Network.sL
    print(new)
    spectrum = Network.lapl_spectrum(h,2)
    print(spectrum)
    
    pass