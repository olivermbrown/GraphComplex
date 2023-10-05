"""
Code to solve the Laplacian on a CW 2-Complex formed from 2 particles on a graph
"""

import numpy as np
import scipy as sp
import sympy as sym
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import squareFDM

class Complex:
    # A class for building a 2-Complex
    
    def __init__(self,cells):
        self.cells = cells
        self.lapl_gen = False
        self.eliminated_vars = []
        self.find_corners()
        self.glueings = []
        self.edges_with_bc_appl = []
        return None
    
    def find_corners(self):
        # Find the locations of every corner in the complex as indices
        
        cells = self.cells
        
        corners = []
        
        c = 0
        
        for cell in cells:
            i = cells.index(cell)
            nodes_list = list(cell.G.nodes())
            cell_corners = np.array([nodes_list.index(cs) for cs in cell.corners])
            cell_corners += c
            corners.append(cell_corners.tolist())
            c += cell.G.number_of_nodes()
            pass
        
        self.corners = reduce(lambda x,y: x+y, corners)
        
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
        
        self.lapl_delete_corners()
        
        return None
    
    def lapl_delete_corners(self):
        # Set to zero any rows and columns in the Laplacian matrix corresponding to corners in cells
        
        L = self.L
        cells = self.cells
        corners = self.corners
        el = self.eliminated_vars
        
        for c in corners:
            L[c] = 0
            L[:,c] = 0
            pass
        
        for cell in cells:
            for corner in cell.corners:
                cell.removed_nodes.append(corner)
                pass
            pass
        
        el += corners
        el.sort()
        self.eliminated_vars = el
        
        L = self.L
        
        return None
    
    def get_edge_coords(self,domain,edge):
        # Get the coordinates in the complex of the locations of nodes on a particular edge
        
        cells = self.cells
        
        nodes_list = list(domain.G.nodes())
        edge_coords = [nodes_list.index(x) for x in edge]
        
        for cell in cells:
            i = cells.index(domain)
            if cells.index(cell) < i:
                t = np.array(edge_coords)
                t += cell.G.number_of_nodes()
                edge_coords = t.tolist()
                pass
            else:
                pass
            pass
        
        return edge_coords
    
    def glue(self,gl):
        # Glue together a list of endpoints
        
        self.check_if_lapl_gen()
        
        L = self.L
        cells = self.cells
        el = self.eliminated_vars
        
        locs = []
        
        blocks = []
        edge_coords = []
        
        for l in gl:
            (domain, edge) = l
            nodes_list = list(domain.G.nodes())
            locs.append([nodes_list.index(x) for x in nodes_list])
            blocks.append(cells.index(domain))
            edge_coords.append([nodes_list.index(x) for x in edge])
            pass
        
        for cell in cells:
            for block in blocks:
                i = blocks.index(block)
                if cells.index(cell) < block:
                    s = np.array(locs[i])
                    t = np.array(edge_coords[i])
                    s += cell.G.number_of_nodes()
                    t += cell.G.number_of_nodes()
                    locs[i] = s.tolist()
                    edge_coords[i] = t.tolist()
                    pass
                else:
                    pass
                pass
            pass
        
        dim = L.shape[0]
        
        L = L.astype('float32')
        
        # Construct vector to implement the gluing condition
        
        v = np.zeros((np.array(edge_coords).shape[1],dim))
        
        n = len(gl)
        
        for l in gl:
            (domain, edge) = l
            i = gl.index(l)
            for x in edge_coords[i]:
                j = edge_coords[i].index(x)
                #
                if edge == domain.x0 or edge == domain.x0inv:
                    s = x + domain.N
                    if s not in el:
                        v[j,s] += 1/n
                        pass
                    else:
                        pass
                    pass
                elif edge == domain.y0 or edge == domain.y0inv:
                    s = x + 1
                    if s not in el:
                        v[j,s] += 1/n
                        pass
                    else:
                        pass
                    pass
                elif edge == domain.x1 or edge == domain.x1inv:
                    s = x - domain.N
                    if s not in el:
                        v[j,s] += 1/n
                        pass
                    else:
                        pass
                    pass
                elif edge == domain.y1 or edge == domain.y1inv:
                    s = x - 1
                    if s not in el:
                        v[j,s] += 1/n
                        pass
                    else:
                        pass
                    pass
                elif edge == domain.hyp:
                    raise Exception
                else:
                    pass
                pass
            pass
        
        # Generate list of all nodes on exterior boundary
        ext = []
        for l in gl:
            (domain, edge) = l
            for e in domain.edges:
                if e == domain.diag:
                    pass
                else:
                    ext += self.get_edge_coords(domain, e)
                    pass
                pass
            pass
        
        ext.sort()
        
        # Eliminate exterior boundary nodes from the gluing vector v - TODO
        i = 0
        for line in v:
            indices = np.nonzero(line)
            for j in indices[0]:
                if j in ext:
                    v[i,j] = 0
                    pass
                else:
                    pass
                pass
            i += 1
            pass
        
        # Apply gluing to the Laplacian matrix
        
        for l in gl:
            (domain, edge) = l
            i = gl.index(l)
            for x in edge_coords[i]:
                j = edge_coords[i].index(x)
                
                nzs = np.nonzero(L[:,x])
                for row in nzs[0]:
                    
                    e = L[row, x]
                    L[row] += e*v[j]
                    L[row, x] = 0
                    pass
                el.append(x)
                L[x] = 0
                L[:,x] = 0
                
            self.edges_with_bc_appl.append((domain,edge))
            pass
        
        el.sort()
        
        self.L = L
        
        self.eliminated_vars = el
        
        self.glueings.append(gl)
        
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
            pass
        
        self.sL = L
        
        return None
    
    def check_if_lapl_gen(self):
        # Check if the Laplacian matrix has been generated and if not, generate it
        
        if self.lapl_gen == True:
            pass
        else:
            self.lapl()
            pass
        
        return None
    
    def apply_dirichlet(self,domain,edge):
        # Apply Dirichlet boundary conditions to a particular edge in the complex
        
        self.check_if_lapl_gen()
        
        N = domain.N
        L = self.L
        cells = self.cells
        el = self.eliminated_vars
        
        nodes_list = list(domain.G.nodes())
        edge_coords = [nodes_list.index(x) for x in edge]
        
        i = cells.index(domain)
        for cell in cells:
            if cells.index(cell) < i:
                t = np.array(edge_coords)
                t += cell.G.number_of_nodes()
                edge_coords = t.tolist()
                pass
            else:
                pass
            pass
        
        for x in edge_coords:
            L[x] = 0
            L[:,x] = 0
            pass
        
        self.L = L
        
        el += edge_coords
        el.sort()
        self.eliminated_vars = el
        
        for node in edge:
            domain.removed_nodes.append(node)
            pass
        
        self.edges_with_bc_appl.append((domain,edge))
        
        return None
    
    def apply_neumann(self,domain,edge):
        # Apply neumann boundary conditions to a particular edge in the complex
        
        self.check_if_lapl_gen()
        
        L = self.L
        cells = self.cells
        el = self.eliminated_vars
        
        nodes_list = list(domain.G.nodes())
        edge_coords = [nodes_list.index(x) for x in edge]
        
        # Identify the correct indices for nodes in the domain based on the block structure of L
        for cell in cells:
            i = cells.index(domain)
            if cells.index(cell) < i:
                t = np.array(edge_coords)
                t += cell.G.number_of_nodes()
                edge_coords = t.tolist()
                pass
            else:
                pass
            pass
        
        N = domain.N
        
        # Modify the Laplacian to account for the boundary conditions
        for i in edge_coords:
            
            if edge == domain.x0:
                j = i + N
                k = i + 2*N
                pass
            if edge == domain.y0:
                j = i + 1
                k = i + 2
                pass
            if edge == domain.x1:
                j = i - N
                k = i - 2*N
                pass
            if edge == domain.y1:
                j = i - 1
                k = i - 2
                pass
            else:
                pass
            
            # First order finite difference
            L[j,j] -= 1
            
            L[i] = 0
            L[:,i] = 0
            el.append(i)
            pass
        
        self.L = L
        
        self.eliminated_vars = el
        
        self.edges_with_bc_appl.append((domain,edge))
        
        return None
    
    def apply_robin(self,domain,edge,a,b):
        # Apply neumann boundary conditions to a particular edge in the complex
        
        self.check_if_lapl_gen()
        
        L = self.L
        cells = self.cells
        el = self.eliminated_vars
        
        
        N = domain.N
        
        # Scaling factor
        h = (np.pi)/(N-1)
        
        # Quotient defining the Robin boundary conditions
        K = (a*h)/b
        
        nodes_list = list(domain.G.nodes())
        edge_coords = [nodes_list.index(x) for x in edge]
        
        # Identify the correct indices for nodes in the domain based on the block structure of L
        for cell in cells:
            i = cells.index(domain)
            if cells.index(cell) < i:
                t = np.array(edge_coords)
                t += cell.G.number_of_nodes()
                edge_coords = t.tolist()
                pass
            else:
                pass
            pass
        
        # Modify the Laplacian to account for the boundary conditions
        for i in edge_coords:
            
            if edge == domain.x0:
                j = i + N
                k = i + 2*N
                pass
            if edge == domain.y0:
                j = i + 1
                k = i + 2
                pass
            if edge == domain.x1:
                j = i - N
                k = i - 2*N
                pass
            if edge == domain.y1:
                j = i - 1
                k = i - 2
                pass
            else:
                pass
            
            # First order finite difference
            s = (-1)/(K - 1)
            L[j,j] -= s
            
            L[i] = 0
            L[:,i] = 0
            el.append(i)
            pass
        
        self.L = L
        
        self.eliminated_vars = el
        
        self.edges_with_bc_appl.append((domain,edge))
        
        return None
    
    def exterior_bc(self,condition):
        # Apply an boundary condition to every edge on the exterior - TODO
        
        cells = self.cells
        
        for cell in cells:
            for edge in cell.edges:
                if isinstance(cell, squareFDM.Domain) and edge != cell.diag:
                    if (cell, edge) not in self.edges_with_bc_appl:
                        if condition == "dirichlet":
                            self.apply_dirichlet(cell,edge)
                            pass
                        elif condition == "neumann":
                            self.apply_neumann(cell,edge)
                            pass
                        elif condition == "robin":
                            self.apply_robin(cell,edge,1,1)
                            pass
                        else:
                            pass
                        pass
                    else:
                        pass
                    pass
                else:
                    pass
                pass
            pass
        
        return None
    
    def diagonal_bc(self,condition):
        # Apply an boundary condition to every edge on the exterior - TODO
        
        cells = self.cells
        
        for cell in cells:
            if isinstance(cell, squareFDM.Domain) and cell.split == True:
                if (cell, cell.diag) not in self.edges_with_bc_appl:
                    if condition == "dirichlet":
                        self.apply_dirichlet(cell,cell.diag)
                        pass
                    elif condition == "neumann":
                        #self.apply_neumann_diag(cell,cell.diag) TODO
                        pass
                    elif condition == "robin":
                        #self.apply_robin_diag(cell,cell.diag,1,1) TODO
                        pass
                    else:
                        pass
                    
                    self.edges_with_bc_appl.append((cell,cell.diag))
                    pass
                else:
                    pass
                pass
            else:
                pass
            pass
        
        return None
    
    def lapl_spectrum(self,h,dps):
        # Calculate the spectrum of a matrix M
        # Scaling factor h
        # Decimal places to round to is dps
        
        self.simplify_lapl()
        matrix = self.sL

        matrix = matrix.astype('float32')
        matrix.tocsr()
        eigs = sp.sparse.linalg.eigs(matrix,which='SM',k=30,return_eigenvectors=False)

        spec = (h)**(-2) * np.abs(eigs)
        spec.sort()
        return np.round(spec,dps)
    
    def lapl_solve(self,h,dps):
        # Calculate the spectrum of a matrix M
        # Scaling factor h
        # Decimal places to round to is dps
        
        self.simplify_lapl()
        matrix = self.sL
        
        matrix = matrix.astype('float32')
        matrix.tocsr()
        eigvals, eigvecs = sp.sparse.linalg.eigs(matrix,which='SM',k=5,return_eigenvectors=True)
        
        spec = (h)**(-2) * np.abs(eigvals)
        spec.sort()
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



def plot_state(D, state):
    # Plot a state that solves the Schrodinger equation on domain D
    
    N = D.N
    #pos = nx.spring_layout(D.G)
    
    coords = []
    for node in D.G:
        if node not in D.removed_nodes:
            x = node[0]/N
            y = node[1]/N
            coords.append([x,y])
            pass
        else:
            pass
        pass
    
    coords = np.array(coords)
    
    
    # Define coordinates
    
    # Plot the surface
    xs = coords[:,0]
    ys = coords[:,1]
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    state = np.real(state)
    ax.plot_trisurf(xs, ys, state, vmin=state.min() * 2, cmap=cm.autumn)
    
    
    return plt



if __name__ == "__main__":
    # Main
    
    N = 50 # The length of a wire
    
    # Scaling factor
    h = (np.pi)/(N-1)
    
    pass



"""
    def diag_elim_dirichlet(self,cell,edge):
        # Apply the elimation of a fourth variable when Dirichlet conditions are applied to a diagonal boundary
        
        self.check_if_lapl_gen()
        
        L = self.L
        el = self.eliminated_vars
        cells = self.cells
        
        edge_coords = self.get_edge_coords(cell,cell.diag)
        
        N = cell.N
        replaced_dict = {}
        
        for loc in edge_coords:
            a = loc - N
            b = loc - 1
            c = loc + 1
            d = loc + N
            vec = np.zeros(L.shape[0])
            if a in replaced_dict.keys():
                vec -= np.array(replaced_dict[a])
                pass
            else:
                vec[a] = - 1
                pass
            if b in replaced_dict.keys():
                vec -= np.array(replaced_dict[b])
                pass
            else:
                vec[b] = - 1
                pass
            if c in replaced_dict.keys():
                vec -= np.array(replaced_dict[c])
                pass
            else:
                vec[c] = - 1
                pass
            nzs = np.nonzero(L[:,d])
            for row in nzs[0]:
                e = L[row, d]
                L[row] += e*vec
                L[row, d] = 0
                pass
            L[d] = 0
            el.append(d)
            replaced_dict[d] = vec
        
        self.L = L
        
        self.eliminated_vars = el
        
        self.diag_replaced_dict = replaced_dict
        
        self.edges_with_bc_appl.append((cell,edge))
        
        return None
    """
    
"""
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
    
    #print(Network.L)
    M = Network.L.toarray()
    #print(M)
    N = Network.simplify_lapl()
    new = Network.sL
    #print(new.toarray())
    spectrum = Network.lapl_spectrum(h,2)
    print(spectrum)
    print("Scaling factor:")
    print(5/spectrum[0])
    print(4*spectrum)
    
    #spectrum, states = Network.lapl_solve(h,2)
    #print(spectrum)
    #n = 4
    #for cell in cells:
    #        plot = plot_state(cell, states[:,n])
    #        plot.show()
    #        pass
    
    """