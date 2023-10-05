"""
Code to solve the Laplacian on a CW 1-Complex
"""

import numpy as np
import scipy as sp
import networkx as nx

import matplotlib.pyplot as plt
from matplotlib import cm
import lineFDM

class Complex:
    # A class for building a 1-Complex
    
    def __init__(self,cells):
        self.cells = cells
        self.lapl_gen = False
        self.removed_nodes = []
        self.eliminated_vars = []
        return None
    
    def construct_indices(self):
        # Assign indices for the complex
        
        return None
    
    def lapl(self):
        # Generate the Laplacian matrix on the unglued complex
        
        cells = self.cells
        for l in cells:
            l.lapl()
            pass
        
        Lapls = tuple([cell.L for cell in cells])
        
        L = sp.sparse.block_diag((Lapls),format="csr")
        
        self.L = L

        self.lapl_gen = True
        
        return None
    
    def glue(self,list):
        # Glue together a list of endpoints
        
        if self.lapl_gen == True:
            pass
        else:
            self.lapl()
            pass
        
        L = self.L
        cells = self.cells
        el = self.eliminated_vars
        
        I = []
        X = []
        
        for l in list:
            (line, end) = l
            I.append(cells.index(line))
            X.append(end)
            pass
        
        for cell in cells:
            for i in range (0,len(I)):
                if cells.index(cell) < I[i]:
                    X[i] += cell.N
                    pass
                else:
                    pass
                pass
            pass
        
        dim = L.shape[0]
        
        L = L.astype('float32')
        
        v = np.zeros(dim)
        
        n = len(X)
        
        for l in list:
            (line, end) = l
            i = list.index(l)
            x = X[i]
            if end == line.start:
                v[x+1] += 1/n
                pass
            elif end == line.end:
                v[x-1] += 1/n
                pass
            else:
                pass
            pass
        
        x0 = X[0]
        for x1 in X[1:]:
            for row in range (0, dim):
                if L[row, x1] != 0:
                    e = L[row, x1]
                    L[row] += e*v
                    L[row, x1] = 0
                    pass
                else:
                    pass
        
        for x in X:
            el.append(x)
            L[x] = 0
            L[:,x] = 0
            pass
        
        el.sort()
        
        self.L = L
        
        self.eliminated_vars = el
        
        return None
    
    def simplify_lapl(self):
        # Remove all zero rows and columns from the Laplacian matrix
        
        el = self.eliminated_vars
        L = self.L
        
        c = 0
        for e in el:
            i = e - c
            L = sp.sparse.vstack([L[:i,:],L[i+1:,:]],format='csr')
            L = sp.sparse.hstack([L[:,:i],L[:,i+1:]],format='csr')
            c += 1
        
        self.sL = L
        
        return None
    
    def lapl_spectrum(self,h,dps):
        # Calculate the spectrum of a matrix M
        # Scaling factor h
        # Decimal places to round to is dps
        
        matrix = self.sL
    
        #matrix = matrix.astype('float32')
        eigs = sp.sparse.linalg.eigs(matrix,which='SM',k=20,return_eigenvectors=False)
    
        spec = (h)**(-1) * np.sqrt(np.abs(eigs))
        #spec = (h)**(-2) * np.abs(eigs)
        spec.sort()
        return np.round(spec,dps)

if __name__ == "__main__":
    # Main
    
    N = 500 # The length of the line-chain
    
    # Scaling factor
    h = (np.pi)/(N-1)
    
    l1 = lineFDM.Line(N)
    l2 = lineFDM.Line(N)
    l3 = lineFDM.Line(N)
    
    cells = [l1,l2,l3]
    
    Network = Complex(cells)
    
    Network.lapl()
    
    print(Network.L)
    M = Network.L.toarray()
    print(M)
    
    gluings = [(l1,l1.start),(l2,l2.start),(l3,l3.start)]

    Network.glue(gluings)
    #print(Network.L)
    M = Network.L.toarray()
    #print(M)
    #print(M[5])
    #print(M[9])
    Network.simplify_lapl()
    SD = Network.sL.toarray()
    print(SD)
    
    spectrum = Network.lapl_spectrum(h,2)
    
    print(spectrum)
    
    pass

"""
def continuity(self,list):
        
        # Apply the continuity part of the gluing map to the Laplacian matrix

        L = self.L
        cells = self.cells
        removed_nodes = self.removed_nodes
        
        I = []
        X = []
        
        for l in list:
            (line, end) = l
            I.append(cells.index(line))
            X.append(end)
            pass
        
        for cell in cells:
            for i in range (0,len(I)):
                if cells.index(cell) < I[i]:
                    X[i] += cell.N
                    pass
                else:
                    pass
                pass
            pass
        
        for rn in removed_nodes:
            for x in X:
                if x > rn:
                    x -= 1
                    pass
                else:
                    pass
        
        X.sort()
        
        dim = L.shape[0]

        # Modify Laplacian to incorporate continuity
        
        x0 = X[0]
        for x1 in X[1:]:
            for row in range (0, dim):
                if L[row, x1] != 0:
                    e = L[row, x1]
                    L[row, x0] = e
                    L[row, x1] = 0
                    pass
                else:
                    pass
        
        for x1 in X[1:]:
            # Remove rows and columns corresponding to eliminated variables
            i = x1
            i -= len(removed_nodes)
            L = sp.sparse.vstack([L[:i,:],L[i+1:,:]],format='csr')
            L = sp.sparse.hstack([L[:,:i],L[:,i+1:]],format='csr')
            removed_nodes.append(x1)
        
        self.L = L
        
        self.removed_nodes = removed_nodes
        
        return None
    
    def current(self,list):
        # Apply the current conservation part of the gluing map to the Laplacian matrix
        
        L = self.L
        cells = self.cells
        removed_nodes = self.removed_nodes
        
        I = []
        X = []
        lines = []
        
        for l in list:
            (line, end) = l
            I.append(cells.index(line))
            X.append(end)
            lines.append(line)
            pass
        
        for cell in cells:
            for i in range (0,len(I)):
                if cells.index(cell) < I[i]:
                    X[i] += cell.N
                    pass
                else:
                    pass
                pass
            pass
        
        # Create vector for current boundary condition
        
        dim = L.shape[0]
        
        L = L.astype('float32')
        
        v = np.zeros(dim)
        
        pos = X
        
        for rn in removed_nodes:
            for x in X:
                if x >= rn:
                    pos[X.index(x)] -= 1
                    pass
                else:
                    pass
                pass
            pass
        
        pos.sort()
        
        for l in list:
            (line, end) = l
            i = list.index(l)
            s = pos[i]
            if end == line.start:
                v[s+1] += 1/3
                pass
            elif end == line.end:
                v[s-1] += 1/3
                pass
            else:
                pass
            pass
        
        x0 = pos[0]
        
        L[x0] = v
        
        L = sp.sparse.vstack([L[:x0,:],L[x0+1:,:]],format='csr')
        L = sp.sparse.hstack([L[:,:x0],L[:,x0+1:]],format='csr')
        
        self.L = L
        
        return None
"""