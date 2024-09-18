"""
Code to solve the Laplacian on a CW 2-Complex formed from 2 particles on a graph
"""

import numpy as np
import scipy as sp
import sympy as sym
import networkx as nx
from functools import reduce
import multiprocessing as mp
import time
from itertools import chain

import matplotlib.pyplot as plt
import matplotlib.tri as tri
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
        self.el_dict = {}
        self.robin_constant = 0

        # Find the corners and edges from each cell in the configuration space.
        self.corners = self.find_corner_indices()
        self.free_edges, self.diagonal = self.get_edges()
        #self.find_edge_indices()

        self.diagonal_boundary_condition = None
        self.exterior_boundary_condition = None
        
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
            if type(cell) == cls.SquareCell:
                for edge in cell.edges:
                    free_edges.append(edge)
                    edge.find_edge_indices_in_configuration_space(self)
                    pass
                pass
            elif type(cell) == cls.TriangleCell:
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

        self.exterior_boundary_condition = condition

        return None
    
    def diagonal_bc(self, condition):
        """
        Set boundary conditions on the diagonal of cells in the configuration space.

        Args:
            self: The ConfigurationSpace object
            condition (str): The type of boundary condition to apply to the diagonal of cells in the configuration space.
        """

        self.diagonal_boundary_condition = condition

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
        """
        Generate the Laplacian matrix for the configuration space.

        Args:
            self: The ConfigurationSpace object
        """

        # Generate the Laplacian matrix on the unglued complex
        self.lapl()

        # Record the start time
        start_time = time.time()

        # Apply the gluing map to the Laplacian matrix
        for gluing in self.gluings:
            self.apply_gluing(gluing)
            pass

        for (gluing, phase) in self.gluings_with_branch_cut:
            self.apply_gluing_with_branch_cut(gluing, phase)
            pass

        # Record the end time
        end_time = time.time()

        # Calculate the time taken
        elapsed_time = end_time - start_time

        print(f"Time taken for gluing: {elapsed_time} seconds")

        # Apply diagonal boundary conditions to the Laplacian matrix
        # TODO: Implement diagonal boundary conditions
        if self.diagonal_boundary_condition == "dirichlet":
            for edge in self.diagonal:
                self.apply_dirichlet(edge)
                pass
            pass
        elif self.diagonal_boundary_condition == "neumann":
            for edge in self.diagonal:
                self.apply_neumann(edge)
                pass
            pass
        elif self.diagonal_boundary_condition == "robin":
            for edge in self.diagonal:
                self.apply_robin(edge)
                pass
            pass
        else:
            pass

        # Apply exterior boundary conditions to the Laplacian matrix
        if self.exterior_boundary_condition == "dirichlet":
            for edge in self.free_edges:
                self.apply_dirichlet(edge)
                pass
            pass
        elif self.exterior_boundary_condition == "neumann":
            for edge in self.free_edges:
                self.apply_neumann(edge)
            pass
        elif self.exterior_boundary_condition == "robin":
            for edge in self.free_edges:
                self.apply_robin(edge)
                pass
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
        
        L[corners] = 0
        L[:,corners] = 0
        
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
    
    def get_cell_coords(self, domain):

        cells = self.cells

        if type(domain) == cls.SquareCell:
            N = domain.N
            coords = np.arange(N**2)
            pass
        elif type(domain) == cls.TriangleCell:
            N = domain.N
            drop_indices = []
            for i in range (1,N):
                drop_indices += [d for d in range(i*N,i*(N+1))]
                pass

            coords = np.arange(N**2 - len(drop_indices))

            pass
        else:
            raise Exception("Unknown domain type")
        
        max = cells.index(domain)

        for i in range(0,max):
            coords += cells[i].G.number_of_nodes()
            pass

        return coords

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

         # Convert el to a sorted set for fast lookups
        el = sorted(set(el))
        
        # Get the total number of rows and columns
        num_rows, num_cols = L.shape
        
        # Create boolean masks for rows and columns to keep
        rows_to_keep = np.ones(num_rows, dtype=bool)
        cols_to_keep = np.ones(num_cols, dtype=bool)
        
        # Mark rows and columns to be removed
        rows_to_keep[el] = False
        cols_to_keep[el] = False
        
        # Use boolean masks to select the rows and columns to keep
        L = L[rows_to_keep, :]
        L = L[:, cols_to_keep]
        
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
        el = self.eliminated_vars
        
        edge_coords = []

        for boundary in gl:
            edge_coords.append(boundary.edge_indices)
            pass
        
        dim = L.shape[0]
        
        L = L.astype('float32')
        
        # Construct vector to implement the gluing condition
        
        v = np.zeros((np.array(edge_coords).shape[1],dim))
        
        n = len(gl)

        # Convert el to a set for faster lookups
        el_set = set(el)

        # NEW CODE START
        results = np.empty((0, 2), dtype=int)
        for i, boundary in enumerate(gl):
            x = np.array(edge_coords[i])
            domain = boundary.domain
            N = domain.N
            if boundary == domain.x0:
                if type(domain) == cls.SquareCell:
                    s = x + N
                    pass
                elif type(domain) == cls.TriangleCell:
                    s = x + N - 1
                    pass
                else:
                    raise Exception("Unknown domain type")
                pass
            elif boundary == domain.y0:
                s = x + 1
                pass
            elif boundary == domain.x1:
                s = x - N
                pass
            elif boundary == domain.y1:
                if type(domain) == cls.SquareCell:
                    s = x - 1
                    pass
                elif type(domain) == cls.TriangleCell:
                    s = x - 1
                    pass
                else:
                    raise Exception("Unknown domain type")
                pass
            elif boundary == domain.diag:
                raise Exception("Invalid boundary: domain.diag")
            else:
                raise Exception("Unknown boundary")
            # Remove any elements of s that are in el_set
            mask = np.isin(s, list(el_set), invert=True)
            s = s[mask]
            # Create an array of integers of length len(edge_coords[i]) with the value i
            j = np.arange(len(edge_coords[i]))
            # Append the results to the list
            new_results = np.column_stack((j, s))
            # np.concatenate((j[:, None], s[:, None]), axis=1)
            results = np.vstack((results, new_results))
            pass
        # NEW CODE END
        
        # Create an array of the same length as indices
        values = np.ones(results.shape[0]) / n

        # Use advanced indexing to update the array `v`
        np.add.at(v, (results[:, 0], results[:, 1]), values)

        # Convert v to a sparse matrix
        v = sp.sparse.csr_matrix(v)
        
        # Generate list of all nodes on exterior boundary
        ch = chain.from_iterable([boundary.edge_indices for boundary in gl])
        ext = list(ch)
        ext.sort()
        
        # Eliminate exterior boundary nodes from the gluing vector v
        ext_set = set(ext)  # Convert ext to a set for faster membership checks

        # Iterate over each line in v
        for i, line in enumerate(v):
            # Get indices of non-zero elements
            non_zero_indices = np.nonzero(line)[0]
            
            # Create a boolean mask for indices that are in ext_set
            mask = np.isin(non_zero_indices, ext_set)
            
            # Set the corresponding entries to zero
            line[non_zero_indices[mask]] = 0
            pass
        
        # Apply gluing to the Laplacian matrix
        
        for i, boundary in enumerate(gl):
            domain = boundary.domain
            edge = boundary.edge

            # NEW CODE START
            x = edge_coords[i]
            j = np.arange(len(x))
            nzs = L[:, x].nonzero()[0]
            nzs = np.array(list(set(nzs)))
            if nzs.size > 0:
                # Get the corresponding values in L for these rows and column x
                values = L[nzs[:,np.newaxis], x]
        
                # Update L rows with the corresponding values times v[j]
                L[nzs] += values @ v[j]
        
                # Set the column x in these rows to zero
                L[nzs[:,np.newaxis], x] = 0
                pass
            else:
                pass

            eliminated_dict = {coord: v[row_idx, :] for row_idx, coord in enumerate(x)}
            #NEW CODE END

            # Mark column and row x as zero
            el += edge_coords[i]
            L[edge_coords[i]] = 0
            L[:, [edge_coords[i]]] = 0

            boundary.eliminated = True
            
            self.edges_with_bc_appl.append(boundary)
            self.free_edges.remove(boundary)
            self.el_dict.update(eliminated_dict)
            pass
        
        el.sort()
        
        self.L = L
        
        self.eliminated_vars = el
        
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
        el = self.eliminated_vars

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
        
        # Convert el to a set for faster lookups
        el_set = set(el)

        # NEW CODE START
        for i, boundary in enumerate(gl):
            x = np.array(edge_coords[i])
            domain = boundary.domain
            N = domain.N
            if boundary == domain.x0:
                if type(domain) == cls.SquareCell:
                    s = x + N
                    pass
                elif type(domain) == cls.TriangleCell:
                    s = x + N - 1
                    pass
                else:
                    raise Exception("Unknown domain type")
                pass
            elif boundary == domain.y0:
                s = x + 1
                pass
            elif boundary == domain.x1:
                s = x - N
                pass
            elif boundary == domain.y1:
                if type(domain) == cls.SquareCell:
                    s = x - 1
                    pass
                elif type(domain) == cls.TriangleCell:
                    s = x - 1
                    pass
                else:
                    raise Exception("Unknown domain type")
                pass
            elif boundary == domain.diag:
                raise Exception("Invalid boundary: domain.diag")
            else:
                raise Exception("Unknown boundary")
            # Remove any elements of s that are in el_set
            mask = np.isin(s, list(el_set), invert=True)
            s = s[mask]
            # Create an array of integers of length len(edge_coords[i]) with the value i
            j = np.arange(len(edge_coords[i]))
            # Append the results to the list
            idx = np.column_stack((j, s))
            if i == n-1:
                v[idx[:,0],idx[:,1]] += phase*1/n
                pass
            else:
                v[idx[:,0],idx[:,1]] += 1/n
                pass
            pass
        # NEW CODE END

        # Generate list of all nodes on exterior boundary
        ch = chain.from_iterable([boundary.edge_indices for boundary in gl])
        ext = list(ch)
        ext.sort()
        
        # Eliminate exterior boundary nodes from the gluing vector v
        ext_set = set(ext)  # Convert ext to a set for faster membership checks

        # Iterate over each line in v
        for i, line in enumerate(v):
            # Get indices of non-zero elements
            non_zero_indices = np.nonzero(line)[0]
            
            # Create a boolean mask for indices that are in ext_set
            mask = np.isin(non_zero_indices, ext_set)
            
            # Set the corresponding entries to zero
            line[non_zero_indices[mask]] = 0
            pass

        # Precompute conjugate of phase as it is used frequently
        conjugate_phase = np.conjugate(phase)
        
        # Apply gluing to the Laplacian matrix
        
        for i, boundary in enumerate(gl):
            domain = boundary.domain

            # NEW CODE START
            x = edge_coords[i]
            j = np.arange(len(x))
            nzs = L[:, x].nonzero()[0]
            nzs = np.array(list(set(nzs)))
            if nzs.size > 0:
                # Get the corresponding values in L for these rows and column x
                values = L[nzs[:,np.newaxis],x]#.toarray().flatten()

                if i == n-1:
                    # Update rows using the conjugate of phase
                    #L[nzs] += conjugate_phase * values @ v[j]
                    v_phase = conjugate_phase * v
                    pass
                else:
                    # Standard update without phase conjugate
                    #L[nzs] += values @ v[j]
                    v_phase = v
                    pass

                L[nzs] += values @ v_phase[j]
        
                # Set the column x in these rows to zero
                L[nzs[:,np.newaxis], x] = 0
                pass
            else:
                pass
            # NEW CODE END

            # Convert v_phase to a sparse matrix
            v_phase = sp.sparse.csr_matrix(v_phase)

            eliminated_dict = {coord: v_phase[row_idx, :] for row_idx, coord in enumerate(x)}

            # Mark column and row x as zero
            el += edge_coords[i]
            L[edge_coords[i]] = 0
            L[:, [edge_coords[i]]] = 0

            boundary.eliminated = True
            
            self.edges_with_bc_appl.append(boundary)
            self.free_edges.remove(boundary)
            self.el_dict.update(eliminated_dict)
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
        
        dirichlet_ids = np.array(edge_coords)

        L[dirichlet_ids] = 0
        L[:,dirichlet_ids] = 0
        
        self.L = L
        
        el += edge_coords
        el.sort()
        self.eliminated_vars = el
        
        #for node in edge:
        #    domain.removed_nodes.append(node)
        #    pass

        # Quicker version of the code with using a loop
        domain.removed_nodes += edge

        boundary.eliminated = True
        
        self.edges_with_bc_appl.append((domain,edge))
        
        return None
    
    def apply_neumann(self, boundary):
        """
        Apply Neumann boundary conditions in the Laplacian matrix to a particular edge in the configuration space.

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
        
        self.L = L
        
        #el += edge_coords
        #el.sort()
        #self.eliminated_vars = el

        # Quicker version of the code with using a loop
        domain.removed_nodes += edge
        
        self.edges_with_bc_appl.append((domain,edge))
        
        return None
    
    def apply_robin(self, boundary):
        """
        Apply Robin boundary conditions in the Laplacian matrix to a particular edge in the configuration space.

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
        c = self.robin_constant

        edge_coords = boundary.edge_indices

        h = (np.pi)/(N-1)

        ones = np.ones(len(edge_coords))

        L[edge_coords,edge_coords] += h * c * ones
        
        self.L = L
        
        #el += edge_coords
        #el.sort()
        #self.eliminated_vars = el

        # Quicker version of the code with using a loop
        domain.removed_nodes += edge
        
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
    
    def lapl_solve(self,h,N_eigs):
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

        self.spectrum = spec
        self.states = eigvecs

        self.solved = True

        self.identify_states_with_cells()

        return np.round(spec,2), eigvecs
    
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

        # Delete diagonal elements if the diagonal boundary condition is neumann or robin
        # delete_nodes = []
        # l = 0
        # if self.diagonal_boundary_condition == "neumann" or self.diagonal_boundary_condition == "robin":
        #     for cell in cells:
        #         coords = np.array(self.get_non_elim_cell_coords(cell))
        #         if type(cell) == cls.SquareCell:
        #             pass
        #         elif type(cell) == cls.TriangleCell:
        #             diag = np.array(cell.diag.edge_indices)
        #             # Find the indices of the diagonal elements in the list coords
        #             indices = np.nonzero(diag[:,None] == coords)[1]
        #             delete_nodes.append(indices+l)
        #             pass
        #         else:
        #             raise Exception("Unknown cell type")
                
        #         l += len(coords)

        #         pass

        #     pass

        # states = np.delete(states,delete_nodes,axis=0)

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
            i += l
            pass

        return None
    
    def construct_coords(self,cell):
        # Construct the non-eliminated coordinate system on each cell

        if cell.cell_coords_constructed == True:
            return None
        else:
            pass

        N = cell.N

        x = np.linspace(0,1,N,endpoint=True)
        y = np.linspace(0,1,N,endpoint=True)

        xx, yy = np.meshgrid(x,y)
        grid = np.array((xx.ravel(), yy.ravel())).T

        if isinstance(cell, cls.TriangleCell):
            # Remove the lower triangle from the coordinate grid.
            g = np.flip(grid,axis=1)

            indices = []
            for j in range(0,N):
                indices += list(range((j)*(N),j*(N)+(j)))
                pass
            
            coords = np.delete(g,indices,axis=0)

            # Remove the nodes that have been eliminated
            delete_nodes = []
            nodes_list = list(cell.G.nodes())
            for e in cell.edges:
                edge_coords = [nodes_list.index(x) for x in e.edge]
                if e.eliminated == True:
                    delete_nodes += edge_coords
                    pass
                else:
                    pass
                pass

            # Add corner nodes to the list of nodes to delete
            corner_coords = [nodes_list.index(x) for x in cell.corners]
            delete_nodes += corner_coords

            coords = np.delete(coords,delete_nodes,axis=0)
            pass
        elif isinstance(cell, cls.SquareCell):
            coords = np.flip(grid,axis=1)

            # Remove the nodes that have been eliminated
            delete_nodes = []
            nodes_list = list(cell.G.nodes())
            for e in cell.edges:
                edge_coords = [nodes_list.index(x) for x in e.edge]
                if e.eliminated == True:
                    delete_nodes += edge_coords
                    pass
                else:
                    pass
                pass

            # Add corner nodes to the list of nodes to delete
            corner_coords = [nodes_list.index(x) for x in cell.corners]
            delete_nodes += corner_coords

            coords = np.delete(coords,delete_nodes,axis=0)
            pass
        else:
            raise Exception
        
        cell.non_elim_coords = coords
        cell.num_non_elim = len(coords)

        cell.cell_coords_constructed = True

        return None
    
    def plot_states(self, n, return_data=False, show_plots=True, solve_lapl=False,plotting_method="surface", realimag="real"):
        # Plot the states of the system

        if solve_lapl == False:
            pass
        elif solve_lapl == True:
            h = (np.pi)/(self.cells[0].N - 1)
            self.lapl_solve(h,dps=2,N_eigs=10)
            pass
        else:
            raise Exception("Invalid value for solve_lapl")
        
        cells = self.cells

        data = []
        
        c = 0
        for cell in cells:
            i = cells.index(cell)

            N = cell.N

            self.construct_coords(cell)

            coords = cell.non_elim_coords
            length = cell.num_non_elim
            indices = cell.indices

            state = cell.eigenstates[:,n]

            if realimag == "real":
                state = np.real(state)
                pass
            elif realimag == "imag":
                state = np.imag(state)
                pass
            elif realimag == "abs":
                state = np.abs(state)
                pass
            else:
                raise Exception("Invalid value for realimag")

            if plotting_method == "surface":
                ax = plt.figure().add_subplot(projection='3d')
                ax.plot_trisurf(coords[:,0],coords[:,1], state, linewidth=0.2, antialiased=True, cmap="Spectral")
                ax.set_zlabel(r'$\psi$'+str(indices)+"(x,y)",rotation=90)
                ax.ticklabel_format(axis='z',style='plain',useOffset=True)
                ax.zaxis.labelpad=10
                plt.title("D"+str(indices))
                pass
            elif plotting_method == "contour":
                triang = tri.Triangulation(coords[:, 0], coords[:, 1])
                fig, ax = plt.subplots()
                contour_plot = ax.tricontourf(triang, state, cmap="cividis")
                contour_lines = ax.tricontour(triang, state, colors='black', linewidths=0.4)
                cbar = fig.colorbar(contour_plot, ax=ax)
                ax.set_title(r'$\psi$'+str(indices)+"(x,y)")
                if type(cell) == cls.TriangleCell:
                    xs = np.linspace(coords[0,0],coords[-1,0],100)
                    if cell.diag.eliminated == True:
                        ys = xs + coords[0,0]
                        pass
                    else:
                        ys = xs
                        pass                    
                    line = ax.plot(xs,ys,color='black',linewidth=0.5)
                    pass
                pass
            else:
                raise Exception("Invalid plotting method")
            
            ax.set_xlabel("x_e"+str(indices[0]))
            ax.set_ylabel("y_e"+str(indices[1]))

            if show_plots == True:
                plt.show()
                pass
            else:
                pass

            plt.clf()
            plt.cla()
            plt.close()
            
            c+=length

            if return_data == True:
                data.append((coords[:,0],coords[:,1],state))
                pass
            else:
                pass
            pass
        
        if return_data == True:
            return data
        else:
            pass

        return None
    
    def save_states(self, filename):
        # Save the states of the system to a file

        directory_path = "C:/Users/vq22287/OneDrive - University of Bristol/Documents/Anyons on graphs/Eigenstates data/"

        filename = directory_path + filename

        cells = self.cells

        c = 0
        for cell in cells:
            self.construct_coords(cell)

            coords = cell.non_elim_coords
            length = cell.num_non_elim

            states = cell.eigenstates
            states = np.real(states)

            #data = np.array((coords[:,0],coords[:,1],states[:])).T
            data = np.hstack((coords, states))

            indices = cell.indices

            np.savetxt(filename+"D"+str(indices)+".csv",data,delimiter=",")
            c+=length
            pass

        return None
    
    def load_states(self, filename):
        # Load the states of the system from a file

        directory_path = "C:/Users/vq22287/OneDrive - University of Bristol/Documents/Anyons on graphs/Eigenstates data/"

        filename = directory_path + filename

        cells = self.cells

        c = 0
        for cell in cells:
            self.construct_coords(cell)

            coords = cell.non_elim_coords
            length = cell.num_non_elim

            data = np.loadtxt(filename+"D"+str(cell.indices)+".csv",delimiter=",")
            states = data[:,2:]

            cell.eigenstates = states

            c+=length
            pass

        return None
    
    def save_eigenvalues(self, filename):
        # Save the eigenvalues of the system to a file

        directory_path = "C:/Users/vq22287/OneDrive - University of Bristol/Documents/Anyons on graphs/Eigenstates data/"

        filename = directory_path + filename

        spec = self.spectrum

        np.savetxt(filename+".csv",spec,delimiter=",")

        return None
    
    def load_eigenvalues(self, filename):
        # Load the eigenvalues of the system from a file

        directory_path = "C:/Users/vq22287/OneDrive - University of Bristol/Documents/Anyons on graphs/Eigenstates data/"

        filename = directory_path + filename

        spec = np.loadtxt(filename+".csv",delimiter=",")

        self.spectrum = spec

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

    C.exterior_boundary_condition("dirichlet")
    C.diagonal_boundary_condition("dirichlet")
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







# def apply_neumann(self, boundary):
    #     """
    #     Apply Neumann boundary conditions in the Laplacian matrix to a particular edge in the configuration space.

    #     Args:
    #         self: The ConfigurationSpace object
    #         boundary (Edge): The edge in the configuration space to apply Dirichlet boundary conditions to.
    #     """
        
    #     self.check_if_lapl_gen()

    #     domain = boundary.domain
    #     edge = boundary.edge
        
    #     N = domain.N
    #     L = self.L
    #     cells = self.cells
    #     el = self.eliminated_vars
    #     el_dict = self.el_dict

    #     edge_coords = boundary.edge_indices
        
    #     edge_indices = np.array(edge_coords)

    #     if type(domain) == cls.SquareCell:
    #         if boundary == domain.x0:
    #             # Apply Neumann boundary conditions to the x0 boundary of the square cell
    #             neumann_ids = edge_indices + N
    #             pass
    #         elif boundary == domain.y0:
    #             # Apply Neumann boundary conditions to the y0 boundary of the square cell
    #             neumann_ids = edge_indices + 1
    #             pass
    #         elif boundary == domain.x1:
    #             # Apply Neumann boundary conditions to the x1 boundary of the square cell
    #             neumann_ids = edge_indices - N
    #             pass
    #         elif boundary == domain.y1:
    #             # Apply Neumann boundary conditions to the y1 boundary of the square cell
    #             neumann_ids = edge_indices - 1
    #             pass
    #         else:
    #             raise Exception("Unknown boundary")
    #         pass
    #     elif type(domain) == cls.TriangleCell:
    #         if boundary == domain.x0:
    #             # Apply Neumann boundary conditions to the x0 boundary of the triangle cell
    #             neumann_ids = edge_indices + (domain.N - 1)
    #             pass
    #         elif boundary == domain.y1:
    #             # Apply Neumann boundary conditions to the y1 boundary of the triangle cell
    #             neumann_ids = edge_indices - 1
    #             pass
    #         elif boundary == domain.diag:
    #             # Apply Neumann boundary conditions to the diagonal boundary of the triangle cell
    #             pass
    #         else:
    #             raise Exception("Unknown boundary")
    #         pass
    #     else:
    #         raise Exception("Unknown domain type")

        
    #     if type(domain) == cls.TriangleCell and boundary == domain.diag:
    #         # Generate a list of all nodes adjacent to the diagonal
            
    #         idx0 = edge_indices[0] - ( N - 1 )
    #         off_diagonal = np.array([idx0])
    #         off_diagonal = np.append(off_diagonal, edge_indices + 1)

    #         # Subtract 1 to the diagonal elements corresponding to the Neumann boundary conditions

    #         # Get the current diagonal elements of L
    #         diagonal = L.diagonal()

    #         # Increment the diagonal elements at positions in neumann_ids
    #         diagonal[off_diagonal[1:-1]] -= 2

    #         # Set the updated diagonal back to the matrix
    #         L.setdiag(diagonal)

    #         # Remove any elements of off_diagonal that are in el
    #         if off_diagonal[0] in el:
    #             # g = el_dict[off_diagonal[0]].toarray().flatten()

    #             # locs = g.nonzero()[0]

    #             # L[locs] += g

    #             # ratio = 1/g[edge_indices[0]]

    #             # print(ratio)

    #             # vec = g * ratio

    #             # vec = vec * 1/(ratio-1)

    #             # vec[edge_indices[0]] = 0

    #             # vec = sp.sparse.csr_matrix(vec)

    #             # print(vec)

    #             # # Remove diagonal[0] from locs
    #             # locs = locs[locs != edge_indices[0]]

    #             # mat = np.tile(vec.toarray(), (len(locs), 1))

    #             # print(mat.shape)

    #             # L[locs] -= mat

    #             col = L[:, edge_indices[0]].toarray()#.flatten()

    #             col[edge_indices[0]] = 0

    #             L[:,off_diagonal[1]] += col

    #             L[edge_indices[0]] = 0

    #             pass
    #         elif off_diagonal[-1] in el:
    #             # g = el_dict[off_diagonal[-1]].toarray().flatten()

    #             # locs = g.nonzero()[0]

    #             # L[locs] += g

    #             # ratio = 1/g[edge_indices[-1]]

    #             # print(ratio)

    #             # vec = g * ratio

    #             # vec = vec * 1/(ratio-1)

    #             # vec[edge_indices[-1]] = 0

    #             # vec = sp.sparse.csr_matrix(vec)

    #             # print(vec)

    #             # # Remove diagonal[-1] from locs
    #             # locs = locs[locs != edge_indices[-1]]

    #             # print(L[locs].shape)
    #             # print(vec.shape)

    #             # mat = np.tile(vec.toarray(), (len(locs), 1))

    #             # print(mat.shape)

    #             # L[locs] -= mat

    #             col = L[:, edge_indices[-1]].toarray()#.flatten()

    #             col[edge_indices[-1]] = 0

    #             L[:,off_diagonal[-2]] += col

    #             L[edge_indices[-1]] = 0

    #             pass
    #         else:
    #             pass

    #         pass
    #     elif type(domain) == cls.TriangleCell or type(domain) == cls.SquareCell:
    #         # Subtract 1 to the diagonal elements corresponding to the Neumann boundary conditions

    #         # Get the current diagonal elements of L
    #         diagonal = L.diagonal()

    #         # Increment the diagonal elements at positions in neumann_ids
    #         diagonal[neumann_ids] -= 1

    #         # Set the updated diagonal back to the matrix
    #         L.setdiag(diagonal)
    #         pass
    #     else:
    #         raise Exception("Unknown domain type")
        

    #     L[edge_indices] = 0
    #     L[:,edge_indices] = 0
        
    #     self.L = L
        
    #     el += edge_coords
    #     el.sort()
    #     self.eliminated_vars = el

    #     # Quicker version of the code with using a loop
    #     domain.removed_nodes += edge
        
    #     self.edges_with_bc_appl.append((domain,edge))
        
    #     return None

    # def apply_neumann_old(self, boundary):
    #     """
    #     Apply Neumann boundary conditions in the Laplacian matrix to a particular edge in the configuration space.

    #     Args:
    #         self: The ConfigurationSpace object
    #         boundary (Edge): The edge in the configuration space to apply Dirichlet boundary conditions to.
    #     """
        
    #     self.check_if_lapl_gen()

    #     domain = boundary.domain
    #     edge = boundary.edge
        
    #     N = domain.N
    #     L = self.L
    #     cells = self.cells
    #     el = self.eliminated_vars
    #     el_dict = self.el_dict

    #     edge_coords = boundary.edge_indices
        
    #     edge_indices = np.array(edge_coords)

    #     if type(domain) == cls.SquareCell:
    #         if boundary == domain.x0:
    #             # Apply Neumann boundary conditions to the x0 boundary of the square cell
    #             neumann_ids = edge_indices + N
    #             pass
    #         elif boundary == domain.y0:
    #             # Apply Neumann boundary conditions to the y0 boundary of the square cell
    #             neumann_ids = edge_indices + 1
    #             pass
    #         elif boundary == domain.x1:
    #             # Apply Neumann boundary conditions to the x1 boundary of the square cell
    #             neumann_ids = edge_indices - N
    #             pass
    #         elif boundary == domain.y1:
    #             # Apply Neumann boundary conditions to the y1 boundary of the square cell
    #             neumann_ids = edge_indices - 1
    #             pass
    #         else:
    #             raise Exception("Unknown boundary")
    #         pass
    #     elif type(domain) == cls.TriangleCell:
    #         if boundary == domain.x0:
    #             # Apply Neumann boundary conditions to the x0 boundary of the triangle cell
    #             neumann_ids = edge_indices + (domain.N - 1)
    #             pass
    #         elif boundary == domain.y1:
    #             # Apply Neumann boundary conditions to the y1 boundary of the triangle cell
    #             neumann_ids = edge_indices - 1
    #             pass
    #         elif boundary == domain.diag:
    #             # Apply Neumann boundary conditions to the diagonal boundary of the triangle cell
    #             pass
    #         else:
    #             raise Exception("Unknown boundary")
    #         pass
    #     else:
    #         raise Exception("Unknown domain type")

        
    #     if type(domain) == cls.TriangleCell and boundary == domain.diag:
    #         # Generate a list of all nodes adjacent to the diagonal
            
    #         idx0 = edge_indices[0] - ( N - 1 )
    #         off_diagonal = np.array([idx0])
    #         off_diagonal = np.append(off_diagonal, edge_indices + 1)

    #         # Remove any elements of off_diagonal that are in el
    #         #off_diagonal = [si for si in off_diagonal if si not in el]
            
    #         # Create a list of pairs of adjacent nodes along the off diagonal
    #         off_diagonal_pairs =  np.column_stack((off_diagonal[:-1], off_diagonal[1:]))

    #         # Convert L to complex
    #         L = L.astype('complex')

    #         # Replace all instances of an amplitude on the diagonal with the average of the adjacent amplitudes on the off-diagonal
    #         for i, idx in enumerate(edge_indices):

    #             # Find the rows where the idx column of L is non-zero
    #             nzs = L[:, idx].nonzero()[0]
    #             nzs = np.array(list(set(nzs)))

    #             if i == 0:
    #                 if idx - (N - 1) in el:
    #                     g = el_dict[idx-(N-1)].toarray().flatten()
    #                     end = off_diagonal[1]

    #                     # Create zero vector of length L.shape[0]
    #                     vec = np.zeros(L.shape[0], dtype='complex128')

    #                     vec[end] = 1/2

    #                     vec += 1/2 * g

    #                     ratio = 1/(1 - vec[idx])
    #                     vec = vec * ratio
    #                     vec[idx] = 0

    #                     vec = sp.sparse.csr_matrix(vec)

    #                     # Update L to account for boundary conditions
    #                     L[end] -= vec
    #                     pass
    #                 else:
    #                     if nzs.size > 0:
    #                         # Find the next component along the off-diagonal
    #                         loc = off_diagonal[1]

    #                         # Update the row in L with the corresponding value
    #                         L[loc, loc] -= 1/2
    #                         pass
    #                     else:
    #                         raise Exception("No non-zero elements in column")
    #                     pass
    #                 pass
    #             elif i == len(edge_indices) - 1:
    #                 # Check if any components of nzs are in the eliminated variables
    #                 if idx + 1 in el:
    #                     g = el_dict[idx+1].toarray().flatten()
    #                     end = off_diagonal[-2]

    #                     # Create zero vector of length L.shape[0]
    #                     vec = np.zeros(L.shape[0], dtype='complex128')

    #                     vec[end] = 1/2

    #                     vec += 1/2 * g

    #                     ratio = 1/(1 - vec[idx])
    #                     vec = vec * ratio
    #                     vec[idx] = 0

    #                     vec = sp.sparse.csr_matrix(vec)

    #                     # Update L to account for boundary conditions
    #                     L[end] -= vec
    #                     pass
    #                 else:
    #                     if nzs.size > 0:
    #                         # Find the next component along the off-diagonal
    #                         loc = off_diagonal[-2]

    #                         # Update the row in L with the corresponding value
    #                         L[loc, loc] -= 1/2
    #                         pass
    #                     else:
    #                         raise Exception("No non-zero elements in column")
    #                     pass
    #                 pass
    #             elif 1 <= i < len(edge_indices) - 1:
    #                 # Find the rows where the ith column of L is non-zero
    #                 if nzs.size > 0:
    #                     # Get the indices of the adjacent nodes on the off-diagonal
    #                     pair = np.array(off_diagonal_pairs[i])

    #                     # Create numpy matrix of ones with dimension nzs size x pair size
    #                     ones = np.ones((nzs.size, pair.size))

    #                     # Update the rows in L with the corresponding values times the ones matrix
    #                     L[nzs[:,np.newaxis], pair] -= 1/2 * ones
    #                     pass
    #                 else:
    #                     pass
    #                 pass
    #             else:
    #                 pass
    #             pass
    #         pass
    #     elif type(domain) == cls.TriangleCell or type(domain) == cls.SquareCell:
    #         # Add 1 to the diagonal elements corresponding to the Neumann boundary conditions

    #         # Get the current diagonal elements of L
    #         diagonal = L.diagonal()

    #         # Increment the diagonal elements at positions in neumann_ids
    #         diagonal[neumann_ids] -= 1

    #         # Set the updated diagonal back to the matrix
    #         L.setdiag(diagonal)
    #         pass
    #     else:
    #         raise Exception("Unknown domain type")
        

    #     L[edge_indices] = 0
    #     L[:,edge_indices] = 0
        
    #     self.L = L
        
    #     el += edge_coords
    #     el.sort()
    #     self.eliminated_vars = el

    #     # Quicker version of the code with using a loop
    #     domain.removed_nodes += edge
        
    #     self.edges_with_bc_appl.append((domain,edge))
        
    #     return None

# def construct_coords_old(self,cell):
    #     # Construct the non-eliminated coordinate system on each cell

    #     if cell.cell_coords_constructed == True:
    #         return None
    #     else:
    #         pass

    #     N = cell.N

    #     x = np.linspace(1/N,1-1/N,N-2,endpoint=True)
    #     y = np.linspace(1/N,1-1/N,N-2,endpoint=True)

    #     xx, yy = np.meshgrid(x,y)
    #     grid = np.array((xx.ravel(), yy.ravel())).T

    #     if isinstance(cell, cls.TriangleCell):
    #         # Remove the lower triangle from the coordinate grid.
    #         g = np.flip(grid,axis=1)

    #         indices = []
    #         for j in range(0,N-2):
    #             indices += list(range((j)*(N-2),j*(N-2)+(j)+1))
    #             pass
            
    #         coords = np.delete(g,indices,axis=0)
    #         pass
    #     elif isinstance(cell, cls.SquareCell):
    #         coords = np.flip(grid,axis=1)
    #         pass
    #     else:
    #         raise Exception
        
    #     cell.non_elim_coords = coords
    #     cell.num_non_elim = len(coords)

    #     cell.cell_coords_constructed = True

    #     return None