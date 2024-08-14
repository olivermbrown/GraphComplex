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

import cells as cls

class ConfigurationSpace:

    def __init__(self,cells):
        """
        Initialize the configuration space with a list of cells.

        Args:
            self: The ConfigurationSpace object
            cells (list): A list of cells in the configuration space
        """

        self.cells = cells
        self.gluings = []
        self.gluings_with_branch_cut = []
        self.lapl_gen = False
        self.edges_with_bc_appl = []
        self.eliminated_vars = []
        self.solved = False

        # Find the corners and edges from each cell in the configuration space.
        self.corners = self.find_corner_indices()
        self.free_edges, self.diagonal = self.get_edges()
        #self.find_edge_indices()
        
        return None

    def find_corner_indices(self):
        """
        Find the indices of every corner in the configuration space.

        Args:
            self: The ConfigurationSpace object
        """

        # Get the list of cells in the configuration space and initialize variables for the function.
        cells = self.cells
        corners = []
        c = 0

        # Iterate over each cell in the configuration space and find the indices of the corners in each cell, appending them to the list of corners in the configuration space.
        for cell in cells:
            i = cells.index(cell)
            nodes_list = list(cell.G.nodes())
            cell_corners = np.array([nodes_list.index(cs) for cs in cell.corners])
            cell_corners += c
            corners.append(cell_corners.tolist())
            c += cell.G.number_of_nodes()
            pass

        corner_indices = reduce(lambda x,y: x+y, corners)

        return corner_indices

    def get_edges(self):
        """
        Find the edges of each cell in the configuration space.

        Args:
            self: The ConfigurationSpace object
        """

        cells = self.cells

        free_edges = []
        diagonal = []

        for cell in cells:
            if isinstance(cell, cls.SquareCell):
                for edge in cell.edges:
                    free_edges.append(edge)
                    edge.find_edge_indices_in_configuration_space(self)
                    pass
                pass
            elif isinstance(cell, cls.TriangleCell):
                for edge in cell.edges:
                    if edge == cell.diag:
                        diagonal.append(edge)
                        pass
                    else:
                        free_edges.append(edge)
                        pass
                    edge.find_edge_indices_in_configuration_space(self)
                    pass
                pass
            else:
                pass
            pass

        return free_edges, diagonal

    def exterior_bc(self, condition):
        """
        Set boundary conditions on the exterior of the configuration space.

        Args:
            self: The ConfigurationSpace object
            condition (str): The type of boundary condition to apply to the exterior of the configuration space.
        """

        self.exterior_bc = condition

        return None
    
    def diagonal_bc(self, condition):
        """
        Set boundary conditions on the diagonal of cells in the configuration space.

        Args:
            self: The ConfigurationSpace object
            condition (str): The type of boundary condition to apply to the diagonal of cells in the configuration space.
        """

        self.diagonal_bc = condition

        return None
    
    def glue(self, gluing):
        """
        Define a gluing of a set of edges in the configuration space.

        Args:
            self: The ConfigurationSpace object
            gluing (list): A list of edges to glue together in the configuration space.        
        """

        self.gluings.append(gluing)

        return None

    def glue_with_branch_cut(self, gluing, phase):
        """
        Define a gluing of a set of edges in the configuration space with a branch cut applied to the final edge in the gluing.

        Args:
            self: The ConfigurationSpace object
            gluing (list): A list of edges to glue together in the configuration space with a branch cut applied to the final edge in the gluing.     
            phase (complex): The phase factor to apply to the final edge in the gluing.
        """

        self.gluings_with_branch_cut.append((gluing, phase))

        return None

    def gen_lapl(self):
        # Generate the Laplacian matrix on the glued complex

        # Generate the Laplacian matrix on the unglued complex
        self.lapl()

        # Apply the gluing map to the Laplacian matrix
        for gluing in self.gluings:
            self.apply_gluing(gluing)
            pass

        for (gluing, phase) in self.gluings_with_branch_cut:
            self.apply_gluing_with_branch_cut(gluing, phase)
            pass

        # Apply diagonal boundary conditions to the Laplacian matrix
        # TODO: Implement diagonal boundary conditions
        if self.diagonal_bc == "dirichlet":
            for edge in self.diagonal:
                self.apply_dirichlet(edge)
                pass
            pass
        elif self.exterior_bc == "neumann":
            pass
        elif self.exterior_bc == "robin":
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
        #self.set_non_elim_cell_coords()

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
            edge_coords.append(boundary.edge_indices)
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
                    if type(domain) == cls.SquareCell:
                        s = x - domain.N
                        pass
                    elif type(domain) == cls.TriangleCell:
                        s = x - domain.N + 1
                        pass
                    else:
                        raise Exception
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
            for e in domain.edges:
                ext += boundary.edge_indices
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
    
    def apply_gluing_with_branch_cut(self, gl, phase):
        """
        Glue together edges in the configuration space with a branch cut applied to the final edge in the gluing.

        Args:
            self: The ConfigurationSpace object
            gl (list): A list of edges to glue together in the configuration space with a branch cut applied to the final edge in the gluing.
        """
        
        self.check_if_lapl_gen()
        
        L = self.L
        cells = self.cells
        el = self.eliminated_vars
        
        locs = []
        
        blocks = []
        edge_coords = []

        # Set phase factor for branch cut

        for boundary in gl:
            edge_coords.append(boundary.edge_indices)
            pass
        
        dim = L.shape[0]
        
        L = L.astype('complex')
        
        # Construct vector to implement the gluing condition
        
        v = np.zeros((np.array(edge_coords).shape[1],dim),dtype=complex)
        
        n = len(gl)
        
        for boundary in gl:
            i = gl.index(boundary)
            domain = boundary.domain
            edge = boundary.edge
            for x in edge_coords[i]:
                j = edge_coords[i].index(x)
                if boundary == domain.x0:# or edge == domain.x0inv:
                    s = x + domain.N
                    if s not in el:
                        if i == n-1:
                            v[j,s] += phase*1/n
                            pass
                        else:
                            v[j,s] += 1/n
                            pass
                        pass
                    else:
                        pass
                    pass
                elif boundary == domain.y0:# or edge == domain.y0inv:
                    s = x + 1
                    if s not in el:
                        if i == n-1:
                            v[j,s] += phase*1/n
                            pass
                        else:
                            v[j,s] += 1/n
                            pass
                        pass
                    else:
                        pass
                    pass
                elif boundary == domain.x1:# or edge == domain.x1inv:
                    s = x - domain.N
                    if s not in el:
                        if i == n-1:
                            v[j,s] += phase*1/n
                            pass
                        else:
                            v[j,s] += 1/n
                            pass
                        pass
                    else:
                        pass
                    pass
                elif boundary == domain.y1:# or edge == domain.y1inv:
                    s = x - 1
                    if s not in el:
                        if i == n-1:
                            v[j,s] += phase*1/n
                            pass
                        else:
                            v[j,s] += 1/n
                            pass
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
            for e in domain.edges:
                ext += boundary.edge_indices
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
                    if i == n-1:
                        L[row] += np.conjugate(phase)*e*v[j]
                        pass
                    else:
                        L[row] += e*v[j]
                        pass
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
                
        return None

    def apply_dirichlet(self, boundary):
        """
        Apply Dirichlet boundary conditions in the Laplacian matrix to a particular edge in the configuration space.

        Args:
            self: The ConfigurationSpace object
            boundary (Edge): The edge in the configuration space to apply Dirichlet boundary conditions to.
        """
        
        self.check_if_lapl_gen()

        domain = boundary.domain
        edge = boundary.edge
        
        N = domain.N
        L = self.L
        cells = self.cells
        el = self.eliminated_vars

        edge_coords = boundary.edge_indices
        
        # nodes_list = list(domain.G.nodes())
        # edge_coords = [nodes_list.index(x) for x in edge]
        
        # i = cells.index(domain)
        # for cell in cells:
        #     if cells.index(cell) < i:
        #         t = np.array(edge_coords)
        #         t += cell.G.number_of_nodes()
        #         edge_coords = t.tolist()
        #         pass
        #     else:
        #         pass
        #     pass
        
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

    def lapl_spectrum(self,h,dps,N_eigs):
        """
        Calculate the spectrum of the Laplacian matrix.

        Args:
            self: The ConfigurationSpace object
            h (float): The scaling factor for the Laplacian matrix
            dps (int): The number of decimal places to round the spectrum to
            N_eigs (int): The number of eigenvalues to calculate
        """
        
        self.simplify_lapl()
        matrix = self.sL

        matrix = matrix.astype('complex64')
        matrix.tocsr()
        eigs = sp.sparse.linalg.eigs(matrix,which='SM',k=N_eigs,return_eigenvectors=False)

        spec = (h)**(-2) * np.abs(eigs)
        spec.sort()
        return np.round(spec,dps)
    
    def lapl_solve(self,h,dps,N_eigs):
        # Calculate the spectrum of a matrix M
        # Scaling factor h
        # Decimal places to round to is dps
        
        self.simplify_lapl()
        matrix = self.sL
        
        matrix = matrix.astype('complex64')
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
            if isinstance(cell, cls.SquareCell):
                cell.eigenstates = states[i:i+l,:]
                pass
            elif isinstance(cell, cls.TriangleCell):
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

        if isinstance(cell, cls.TriangleCell):
            # Remove the upper triangle from the coordinate grid.
            g = np.flip(grid,axis=1)

            indices = []
            for j in range(0,N-2):
                indices += list(range(j*(N-2)+(j),(j+1)*(N-2)))
                pass
            
            c = np.delete(g,indices,axis=0)
            coords = np.flip(c,axis=1)
            pass
        elif isinstance(cell, cls.SquareCell):
            coords = grid
            pass
        else:
            raise Exception
            pass
        
        cell.non_elim_coords = coords
        cell.num_non_elim = len(coords)

        return None
    
    def plot_states(self, n,):
        # Plot the states of the system

        if self.solved == True:
            pass
        else:
            h = (np.pi)/(self.cells[0].N - 1)
            self.lapl_solve(h,dps=2,N_eigs=10)
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



if __name__ == "__main__":
    # Main

    N = 20

    h = (np.pi)/(N-1)
    
    D11 = cls.TriangleCell(N)
    #D12 = cls.SquareCell(N)
    D21 = cls.SquareCell(N)
    D22 = cls.TriangleCell(N)

    D11.indices = (1,1)
    D21.indices = (2,1)
    D22.indices = (2,2)

    C = ConfigurationSpace([D11,D21,D22])
    #C = ConfigurationSpace([D21])

    gluing1 = [D11.x1,D21.x0]
    gluing2 = [D21.y1,D22.y0]
    C.glue(gluing1)
    C.glue(gluing2)

    C.exterior_bc("dirichlet")
    C.diagonal_bc("dirichlet")
    C.gen_lapl()
    lapl = C.sL
    #print(lapl.toarray())
    #C.print_eqs()

    #spec = C.lapl_spectrum(h,2,10)
    #print(spec)

    C.plot_states(0,10)
    spec = C.spectrum
    spec.sort()
    print(spec)



    pass
