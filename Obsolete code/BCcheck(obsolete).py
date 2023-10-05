"""
Code to check that the boundary conditions on a 2-Complex are being implemented correctly
"""

import numpy as np
import scipy as sp
import sympy as sym

import squareFDM
import triangleFDM
from twocomplex import Complex

def gen_vector(Complex):
    # Generate a symbol vector of function values
    
    L = Complex.L
    
    dim = L.shape[0]
    
    els = []
    
    for i in range (0,dim):
        symbol = "f" + str(i)
        els.append(symbol)
        pass
    
    v = sym.Matrix(els)
    
    return v

if __name__ == "__main__":
    # Main
    
    N = 10 # The length of a wire
    
    # Scaling factor
    h = (np.pi)/(N-1)
    
    # Build the Complex from two square domains
    D11 = squareFDM.Domain(N)
    D12 = squareFDM.Domain(N)
    D21 = squareFDM.Domain(N)
    D22 = squareFDM.Domain(N)
    D11.split_domain("dirichlet")
    D22.split_domain("dirichlet")
    
    cells = [D11,D12,D21,D22]
    #cells = [D11]
    
    Network = Complex(cells)
    
    gluingA = [(D11, D11.x1),(D12, D12.x0)]
    gluingB = [(D12, D12.y1),(D22, D22.y0)]
    gluingC = [(D22, D22.x0),(D21, D21.x1)]
    gluingD = [(D21, D21.y0),(D11, D11.y1)]
    
    Network.diagonal_bc("dirichlet")
    
    Network.glue(gluingA)
    Network.glue(gluingB)
    Network.glue(gluingC)
    Network.glue(gluingD)
    
    Network.exterior_bc("dirichlet")
    
    #M = Network.L.toarray()
    
    Network.print_eqs()
    
    """
    v = gen_vector(Network)
    
    LHS = M * v
    
    c = 0
    for row in LHS:
        print("Equation " + str(c) + ":")
        print(row)
        c += 1
    """
    
    pass