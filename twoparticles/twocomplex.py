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
        self.gluings = []
        self.find_edges()
        self.edges_with_bc_appl = []
        self.solved = False
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
    
    def find_edges(self):
        # Find the edges in the complex

        cells = self.cells

        free_edges = []
        diagonal = []

        for cell in cells:
            if isinstance(cell, squareFDM.Domain):
                for edge in cell.edges:
                    if edge == cell.diag:
                        diagonal.append(edge)
                        pass
                    else:
                        free_edges.append(edge)
                        pass
                    pass
                pass
            else:
                pass
            pass

        self.free_edges = free_edges
        self.diagonal = diagonal

        return None
    
    def exterior_bc(self, condition):
        # Set boundary conditions on the exterior of the complex

        self.exterior_bc = condition

        return None
    
    def diagonal_bc(self, condition):
        # Set boundary conditions on the diagonal of domains in the complex

        self.diagonal_bc = condition

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
            self.apply_gluing(gluing)
            pass

        # Apply diagonal boundary conditions to the Laplacian matrix
        if self.diagonal_bc == "dirichlet":
            for edge in self.diagonal:
                if edge.domain.split == True:
                    self.apply_dirichlet(edge)
                    pass
                else:
                    pass
                pass
            pass
        elif self.diagnonal_bc == "neumann":
            pass
        elif self.diagnonal_bc == "robin":
            pass
        else:
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

        # Set the non-eliminated coordinates on each cell in the complex
        self.set_non_elim_cell_coords()

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
    
    def get_edge_coords(self,boundary):
        # Get the locations in the complex of the locations of nodes on a particular edge
        
        cells = self.cells

        domain = boundary.domain
        edge = boundary.edge
        
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
    
    def get_cell_coords(self, domain):
        # Get the locations in the complex of the locations of nodes on a particular cell

        cells = self.cells

        cell_coords = list(range(0,domain.G.number_of_nodes()))

        for cell in cells:
            i = cells.index(domain)
            if cells.index(cell) < i:
                t = np.array(cell_coords)
                t += cell.G.number_of_nodes()
                cell_coords = t.tolist()
                pass
            else:
                pass
            pass

        return cell_coords
    
    def get_non_elim_cell_coords(self, domain):
        # Get the locations in the complex of the non-eliminated nodes on a particular cell

        coords = self.get_cell_coords(domain)

        el = self.eliminated_vars

        el_set = set(el)

        # Remove eliminated nodes from the coordinate list
        non_elim_nodes = [element for element in coords if element not in el_set]

        return non_elim_nodes
    
    def set_non_elim_cell_coords(self):
        # Set the non-eliminated coordinates on each cell in the complex

        cells = self.cells

        for cell in cells:
            non_elim_nodes = self.get_non_elim_cell_coords(cell)
            cell.non_elim_nodes = non_elim_nodes
            pass

        return None
    
    def apply_gluing(self,gl):
        # Glue together a list of endpoints
        
        self.check_if_lapl_gen()
        
        L = self.L
        cells = self.cells
        el = self.eliminated_vars
        
        locs = []
        
        blocks = []
        edge_coords = []
        
        for boundary in gl:
            domain = boundary.domain
            edge = boundary.edge
            nodes_list = list(domain.G.nodes())
            locs.append([nodes_list.index(x) for x in nodes_list])
            blocks.append(cells.index(domain))
            edge_coords.append([nodes_list.index(x) for x in edge])
            pass
        
        for cell in cells:
            for block in list(set(blocks)):
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
        
        for boundary in gl:
            i = gl.index(boundary)
            domain = boundary.domain
            edge = boundary.edge
            for x in edge_coords[i]:
                j = edge_coords[i].index(x)
                #
                if boundary == domain.x0:# or edge == domain.x0inv:
                    s = x + domain.N
                    if s not in el:
                        v[j,s] += 1/n
                        pass
                    else:
                        pass
                    pass
                elif boundary == domain.y0:# or edge == domain.y0inv:
                    s = x + 1
                    if s not in el:
                        v[j,s] += 1/n
                        pass
                    else:
                        pass
                    pass
                elif boundary == domain.x1:# or edge == domain.x1inv:
                    s = x - domain.N
                    if s not in el:
                        v[j,s] += 1/n
                        pass
                    else:
                        pass
                    pass
                elif boundary == domain.y1:# or edge == domain.y1inv:
                    s = x - 1
                    if s not in el:
                        v[j,s] += 1/n
                        pass
                    else:
                        pass
                    pass
                elif boundary == domain.diag:
                    raise Exception
                else:
                    pass
                pass
            pass
        
        # Generate list of all nodes on exterior boundary
        ext = []
        for boundary in gl:
            domain = boundary.domain
            edge = boundary.edge
            for e in domain.edges:
                if e == domain.diag:
                    pass
                else:
                    ext += self.get_edge_coords(e)
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
        
        for boundary in gl:
            i = gl.index(boundary)
            domain = boundary.domain
            edge = boundary.edge
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
                pass
            
            self.edges_with_bc_appl.append(boundary)
            self.free_edges.remove(boundary)
            pass
        
        el.sort()
        
        self.L = L
        
        self.eliminated_vars = el
        
        #self.glueings.append(gl)
        
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
    
    def apply_dirichlet(self,boundary):
        # Apply Dirichlet boundary conditions to a particular edge in the complex
        
        self.check_if_lapl_gen()

        domain = boundary.domain
        edge = boundary.edge
        
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

    def lapl_spectrum(self,h,dps,N_eigs):
        # Calculate the spectrum of a matrix M
        # Scaling factor h
        # Decimal places to round to is dps
        
        matrix = self.sL

        matrix = matrix.astype('float32')
        matrix.tocsr()
        eigs = sp.sparse.linalg.eigs(matrix,which='SM',k=N_eigs,return_eigenvectors=False)

        spec = (h)**(-2) * np.abs(eigs)

        self.spectrum = np.round(spec,dps)

        spec.sort()
        return np.round(spec,dps)
    
    def lapl_solve(self,h,dps,N_eigs):
        # Calculate the spectrum of a matrix M
        # Scaling factor h
        # Decimal places to round to is dps
        
        matrix = self.sL
        
        matrix = matrix.astype('float32')
        matrix.tocsr()
        eigvals, eigvecs = sp.sparse.linalg.eigs(matrix,which='SM',k=N_eigs,return_eigenvectors=True)
        
        spec = (h)**(-2) * np.abs(eigvals)
        #spec.sort()

        self.spectrum = np.round(spec,dps)
        self.states = eigvecs

        self.solved = True

        self.identify_states_with_cells()

        return np.round(spec,dps), eigvecs
    
    def print_eqs(self):
        # Print out the system of equations described by the Laplacian matrix
        # to check that the boundary conditions have been correctly implemented
        
        # Generate a symbol vector of function values
    
        L = self.L
        
        dim = L.shape[0]

        N = self.cells[0].N
        
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
            if c % N**2 == 0:
                # Wait for user input every 100 equations
                input("Press Enter to continue...")
                pass
            pass
        
        return None
    
    def identify_states_with_cells(self):
        # Identify the components of the solved eigenstates with each cell in the complex

        cells = self.cells

        if self.solved == True:
            pass
        else:
            h = (np.pi)/(cells[0].N - 1)
            self.lapl_solve(h,2,dps=2)
            pass

        states = self.states

        i = 0
        for cell in cells:

            self.construct_coords(cell)
            l = cell.num_non_elim
            if cell.split == True:
                cell.eigenstates = states[i:i+l,:]
                pass
            elif cell.split == False:
                cell.eigenstates = states[i:i+l,:]
                pass
            else:
                raise Exception
                pass
            i += l
            pass


        return None
    
    def construct_coords(self,cell):
        # Construct the non-eliminated coordinate system on each cell

        N = cell.N

        x = np.linspace(0,1,N-2,endpoint=False)
        y = np.linspace(0,1,N-2,endpoint=False)

        xx, yy = np.meshgrid(x,y)
        grid = np.array((xx.ravel(), yy.ravel())).T

        if cell.split == True:
            # Remove the diagonal from the coordinate grid
            diag = np.stack((x,x),axis=1)
            indices = np.where((grid==diag[:,None]).all(-1))[1]
            coords = np.delete(grid,indices,axis=0)
            pass
        else:
            coords = grid
            pass
        
        cell.non_elim_coords = coords
        cell.num_non_elim = len(coords)

        return None
    
    def plot_states(self, n):
        # Plot the states of the system

        if self.solved == True:
            pass
        else:
            self.lapl_solve(2,dps=2)
            pass
        
        cells = self.cells
        states = self.states
        
        c = 0
        for cell in cells:
            i = cells.index(cell)

            self.construct_coords(cell)

            coords = cell.non_elim_coords
            length = cell.num_non_elim

            #state = states[c:length+c,n]
            state = cell.eigenstates[:,n]
            state = np.real(state)

            ax = plt.figure().add_subplot(projection='3d')

            ax.plot_trisurf(coords[:,0],coords[:,1], state, linewidth=0.2, antialiased=True, cmap=cm.inferno)
            #ax.scatter(coords[:,0],coords[:,1],state,c=state,cmap=cm.autumn)
            indices = cell.indices

            plt.title("D"+str(indices))
            ax.set_xlabel("x_e"+str(indices[0]))
            ax.set_ylabel("y_e"+str(indices[1]))
            ax.set_zlabel(r'$\psi$'+str(indices)+"(x,y)",rotation=90)
            ax.ticklabel_format(axis='z',style='plain',useOffset=True)
            ax.zaxis.labelpad=10
            plt.show()
            c+=length
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
    
    pass
