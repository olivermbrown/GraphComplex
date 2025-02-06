"""
Code to solve the Laplacian on a CW 1-Complex
"""

import numpy as np
import scipy as sp
import sympy as sym
import networkx as nx

import matplotlib.pyplot as plt
from matplotlib import cm
from itertools import groupby
import lineFDM

class Complex:
    # A class for building a 1-Complex
    
    def __init__(self,cells):
        self.cells = cells
        #self.create_domain_dict()
        self.lapl_gen = False
        self.removed_nodes = []
        self.eliminated_vars = []
        self.h = 1
        self.find_scaling()
        self.edges_with_bc_appl = []
        self.gluings = []
        self.solved = False
        return None

    def set_scaling(self,N):
        # Set a scaling factor for the Laplacian

        self.h = (np.pi)/(N-1)

        return None
    
    def all_equal(self, iterable):
        # Check if all elements in an iterable are equal

        g = groupby(iterable)

        return next(g, True) and not next(g, False)
    
    def find_scaling(self):
        # If all cells are of the same length,
        # calculate the scaling factor for the Laplacian

        N = self.cells[0].N

        N_same = self.all_equal([cell.N for cell in self.cells])

        if N_same == True:
            self.set_scaling(N)
            pass
        else:
            pass

        return None
    
    def construct_indices(self):
        # Assign indices for the complex

        # TODO - remove?
        
        return None

    def exterior_bc(self, condition):
        # Set boundary conditions on the exterior of the complex

        self.exterior_bc = condition

        return None
    
    def glue(self, gluing):
        # Set gluing to a particular edge in the complex

        self.gluings.append(gluing)

        return None
    
    def gen_lapl(self):
        # Generate the Laplacian matrix on the glued complex

        # Generate the Laplacian matrix on the unglued complex
        self.lapl()

        # Apply the gluing map to the Laplacian matrix
        for gluing in self.gluings:
            #g = []
            #g.append(gluing.line, gluing.node)
            self.apply_gluing(gluing)
            pass

        # Apply exterior boundary conditions to the Laplacian matrix
        if self.exterior_bc == "dirichlet":
            for edge in self.free_edges:
                self.apply_dirichlet(edge)
                pass
            pass
        elif self.exterior_bc == "neumann":
            pass
        elif self.exterior_bc == "robin":
            pass
        else:
            pass

        # Remove all zero rows and columns from the Laplacian matrix
        self.simplify_lapl()

        return None
    
    def lapl(self):
        # Generate the Laplacian matrix on the unglued complex
        
        cells = self.cells
        for l in cells:
            l.lapl()
            pass
        
        Lapls = tuple([cell.L for cell in cells])
        
        L = sp.sparse.block_diag((Lapls),format="lil")
        
        self.L = L

        self.lapl_gen = True

        self.free_edges = []

        for cell in cells:
            for edge in cell.endpoints:
                if edge not in self.edges_with_bc_appl:
                    self.free_edges.append(edge)
                    pass
                else:
                    pass
                pass
            pass
        
        return None
    
    def apply_dirichlet(self, endpoint):
        # Apply Dirichlet boundary conditions to a particular edge in the complex

        line = endpoint.line
        end = endpoint.end
        
        N = line.N
        L = self.L
        cells = self.cells
        el = self.eliminated_vars
        
        nodes_list = list(line.G.nodes())
        
        i = cells.index(line)
        for cell in cells:
            if cells.index(cell) < i:
                t = end
                t += cell.G.number_of_nodes()
                end = t
                pass
            else:
                pass
            pass
        
        
        L[end] = 0
        L[:,end] = 0
        
        self.L = L
        
        el.append(end)
        el.sort()
        self.eliminated_vars = el

        self.edges_with_bc_appl.append(endpoint)
        
        return None
    
    def apply_gluing(self,list):
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
            line = l.line
            end = l.end
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
            line = l.line
            end = l.end
            i = list.index(l)
            x = X[i]
            if end == line.start.end:
                v[x+1] += 1/n
                pass
            elif end == line.end.end:
                v[x-1] += 1/n
                pass
            else:
                pass
            pass

        # Convert L to csr matrix
        L = L.tocsr()
        
        for x1 in X:
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

        for l in list:
            self.edges_with_bc_appl.append(l)
            self.free_edges.remove(l)
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
    
    def lapl_spectrum(self,dps):
        # Calculate the spectrum of a matrix M
        # Decimal places to round to is dps
        
        matrix = self.sL
        h = self.h
    
        #matrix = matrix.astype('float32')
        eigs = sp.sparse.linalg.eigs(matrix,which='SM',k=20,return_eigenvectors=False)
        
        spec = (h)**(-1) * np.sqrt(np.abs(eigs))
        spec.sort()

        self.spectrum = spec

        return np.round(spec,dps)
    
    def lapl_solve(self,dps,n_eigs):
        # Calculate the spectrum of a matrix M
        # Scaling factor h
        # Decimal places to round to is dps
        
        matrix = self.sL
        h = self.h

        matrix = matrix.astype('float32')
        matrix.tocsr()
        eigvals, eigvecs = sp.sparse.linalg.eigs(matrix,which='SM',k=n_eigs,return_eigenvectors=True)
        
        spec = (h)**(-1) * np.sqrt(np.abs(eigvals))

        self.spectrum = spec
        self.states = eigvecs
        self.solved = True

        return np.round(spec,dps), eigvecs
    
    def print_eqs(self):
        # Print out the system of equations described by the Laplacian matrix
        # to check that the boundary conditions have been correctly implemented
        
        # Generate a symbol vector of function values
    
        L = self.L
        
        dim = L.shape[0]
        
        els = []
        
        for i in range (0,dim):
            symbol = "f" + str(i)
            els.append(symbol)
            pass
        
        v = sym.Matrix(els)
        
        M = self.L.toarray()
        
        LHS = M * v
    
        c = 0
        for row in LHS:
            print("Equation " + str(c) + ":")
            print(row)
            c += 1
            pass
        
        return None
    
    def plot_states(self, n):
        # Plot the states of the system

        if self.solved == True:
            pass
        else:
            self.lapl_solve(2)
            pass
        
        cells = self.cells
        states = self.states
        
        for line in cells:
            i = cells.index(line)

            N = line.N
            coords = np.linspace(0,1,N-1,endpoint=False)

            xs = coords[1:]
            length = len(xs)
            state = states[i*length:(i+1)*length,n]
            state = np.real(state)

            plt.plot(xs, state)
            plt.xlabel("x")
            plt.ylabel(r'$\psi$'+str(i)+"(x)")
            plt.show()
            pass
        
        return None

if __name__ == "__main__":
    # Main
    
    pass
