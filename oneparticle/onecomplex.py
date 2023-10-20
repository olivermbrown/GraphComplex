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
        return None
    
    """
    def create_domain_dict(self):
        # Create a dictionary referencing cells and endpoints
        
        domain_dict = {}    # Empty dictionary of cells and endpoints

        for cell in self.cells:
            for e in cell.endpoints:
                domain_dict[e] = cell
                pass
            pass

        self.domain_dict = domain_dict

        return None
    """

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
        
        return None
    
    """
    def set_bc(self, end, condition):
        # Set boundary conditions on a particular edge in the complex

        if condition == "dirichlet":
            pass
        elif condition == "neumann":
            pass
        elif condition == "robin":
            pass

        return None
    """

    def exterior_bc(self, condition):
        # Set boundary conditions on the exterior of the complex

        cells = self.cells

        self.free_edges = []

        self.exterior_bc = condition

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
    
    def glue(self, gluing):
        # Set gluing to a particular edge in the complex

        #dict = self.domain_dict

        #g = []

        #for e in gluing:
        #    g.append((dict[e],e))
        #    pass

        self.gluings.append(gluing)

        return None
    
    def gen_lapl(self):
        # Generate the Laplacian matrix on the glued complex

        # Import domain dictionary
        #dict = self.domain_dict

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
                line = edge.line
                end = edge.node
                self.apply_dirichlet(line, end)
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
        
        return None
    
    def apply_dirichlet(self, endpoint):
        # Apply Dirichlet boundary conditions to a particular edge in the complex
        
        N = line.N
        L = self.L
        cells = self.cells
        el = self.eliminated_vars

        line = endpoint.line
        end = endpoint.node
        
        nodes_list = list(line.G.nodes())
        #edge_coords = [nodes_list.index(x) for x in end]
        
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
        """
        for node in edge:
            domain.removed_nodes.append(node)
            pass
        """
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
            end = l.node
            #(line, end) = l
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
            end = l.node
            #(line, end) = l
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
            #(line, end) = l
            self.edges_with_bc_appl.append(l)
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
        return np.round(spec,dps)
    
    def lapl_solve(self,dps):
        # Calculate the spectrum of a matrix M
        # Scaling factor h
        # Decimal places to round to is dps
        
        self.simplify_lapl()
        matrix = self.sL
        h = self.h

        matrix = matrix.astype('float32')
        matrix.tocsr()
        eigvals, eigvecs = sp.sparse.linalg.eigs(matrix,which='SM',k=20,return_eigenvectors=True)
        
        spec = (h)**(-1) * np.sqrt(np.abs(eigvals))
        #spec.sort()
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

def plot_state(l, state):
    # Plot a state that solves the Schrodinger equation on wire l
    
    N = l.N
    #pos = nx.spring_layout(D.G)
    
    """
    coords = []
    for node in l.G:
        if node not in l.removed_nodes:
            x = node/N
            #y = node[1]/N
            coords.append([x])
            pass
        else:
            pass
        pass
    """
    
    # Define coordinates
    #coords = np.array(coords)
    coords = np.linspace(0,1,N-1,endpoint=False)
    
    # Plot the surface
    xs = coords[1:]
    #ys = coords[:,1]
    
    # Plot 2d line
    fig, ax = plt.subplots()
    #fig, ax = plt.subplots(subplot_kw={"projection": "2d"})
    state = np.real(state)
    ax.plot(xs, state)
    #ax.plot_trisurf(xs, state, vmin=state.min() * 2, cmap=cm.autumn)
    
    
    return plt

if __name__ == "__main__":
    # Main
    
    N = 50 # The length of the line-chain
    
    # Scaling factor
    h = (np.pi)/(N-1)
    
    l1 = lineFDM.Line(N)
    l2 = lineFDM.Line(N)
    l3 = lineFDM.Line(N)
    
    cells = [l1,l2,l3]
    #cells = [l1,l2]
    
    Network = Complex(cells)
    
    Network.lapl()
    
    print(Network.L)
    M = Network.L.toarray()
    print(M)
    
    gluings = [(l1,l1.start),(l2,l2.start),(l3,l3.start)]
    #gluings = [(l1,l1.start),(l2,l2.start)]
    
    Network.apply_dirichlet(l1,l1.end)
    Network.apply_dirichlet(l2,l2.end)
    Network.apply_dirichlet(l3,l3.end)

    Network.glue(gluings)
    #print(Network.L)
    M = Network.L.toarray()
    #print(M)
    #print(M[5])
    #print(M[9])
    Network.simplify_lapl()
    SD = Network.sL.toarray()
    #print(SD)
    
    Network.print_eqs()
    
    spectrum, states = Network.lapl_solve(h,2)

    print(spectrum)
    
    # Choose state to plot
    n = 3
    vector = states[:,n]
    states = np.array_split(vector, 3)
    
    for line in cells:
        i = cells.index(line)
        plot = plot_state(line, states[i])
        plot.show()
        pass
    
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