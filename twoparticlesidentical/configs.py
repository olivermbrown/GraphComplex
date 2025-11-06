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
import pandas as pd

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

        self.plot_dim = (1,1)
        self.plots_off = []
        self.figuresize = (8,6)

        self.diagonal_boundary_condition = None
        self.exterior_boundary_condition = None

        self.alpha = 0

        mins = []
        for cell in cells:
            mins.append(np.min((cell.N, cell.N2)))
        pass
        self.Nmin = np.min(mins)
        
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
            N1 = domain.N
            N2 = domain.N2
            coords = np.arange(N1*N2)
            pass
        elif type(domain) == cls.TriangleCell:
            N1 = domain.N
            N2 = domain.N2

            # Find the indices of the nodes along the diagonal
            diag_ids = []
            for x in np.linspace(0,(N1*N2 - 1), N2):
                if x == 0:
                    continue
                else:
                    pass
                # Round x to the nearest integer
                idx = int(round(x))
                diag_ids.append(idx)
                pass

            # Find all nodes on the lower right half of the triangle
            drop_indices = []
            for i in range (1,N2+1):
                drop_indices += [d for d in range(i*N1,diag_ids[i-1])]
                pass

            coords = np.arange(N1*N2 - len(drop_indices))

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
    
    def apply_neumann_strong(self, boundary):
        """
        Apply strong Neumann boundary conditions in the Laplacian matrix to a particular edge in the configuration space.

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

        off_diagonal = np.array(edge_coords[:-1]) + 1

        neumann_ids = np.concatenate((np.array(edge_coords), off_diagonal))

        L[neumann_ids] = 0
        L[:,neumann_ids] = 0
        
        self.L = L
        
        el += list(neumann_ids)
        el.sort()
        self.eliminated_vars = el

        boundary.eliminated = True
        boundary.strongly_eliminated = True

        # Quicker version of the code with using a loop
        domain.removed_nodes += list(neumann_ids)
        
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
    
    def apply_robin_strong(self, boundary):
        """
        Apply strong Robin boundary conditions in the Laplacian matrix to a particular edge in the configuration space.

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

        L[edge_coords,edge_coords] = 1/(1+np.exp(1j*np.pi*self.alpha)) * L[edge_coords,edge_coords]
        
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
        
        #matrix = matrix.astype('complex64')
        matrix = matrix.astype('complex128')
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
        N2 = cell.N2

        min = self.Nmin
        L1 = N/min
        L2 = N2/min

        x = np.linspace(0,L1,N,endpoint=True)
        y = np.linspace(0,L2,N2,endpoint=True)

        xx, yy = np.meshgrid(x,y)
        grid = np.array((xx.ravel(), yy.ravel())).T

        if isinstance(cell, cls.TriangleCell):
            # Remove the lower triangle from the coordinate grid.
            g = np.flip(grid,axis=1)

            indices = []
            #for j in range(0,N):
            #    indices += list(range((j)*(N),j*(N)+(j)))
            #    pass
            
            indices = cell.dropped_indices
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
                if e.strongly_eliminated == True and e == cell.diag:
                    off_diagonal = np.array(edge_coords[:-1]) + 1
                    delete_nodes += list(off_diagonal)
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
    
    def plot_states(self, n, return_data=False, show_plots=True, solve_lapl=False, plotting_method="surface", realimag="abs", N_levels=20):
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

        # Calculate contour grading
        if plotting_method == "contour":
            mins = []
            maxs = []
            for cell in cells:
                if realimag == "phase":
                    #min = np.min(np.angle(cell.eigenstates[:,n]))
                    #max = np.max(np.angle(cell.eigenstates[:,n]))
                    min = -np.pi
                    max = np.pi
                    pass
                elif realimag == "abs":
                    min = np.min(abs(cell.eigenstates[:,n]))
                    max = np.max(abs(cell.eigenstates[:,n]))
                    pass
                elif realimag == "real":
                    min = np.min(np.real(cell.eigenstates[:,n]))
                    max = np.max(np.real(cell.eigenstates[:,n]))
                    pass
                elif realimag == "imag":
                    min = np.min(np.imag(cell.eigenstates[:,n]))
                    max = np.max(np.imag(cell.eigenstates[:,n]))
                    pass
                else:
                    min = np.min(cell.eigenstates[:,n])
                    max = np.max(cell.eigenstates[:,n])
                    pass
                mins.append(min)
                maxs.append(max)
                pass
            min = np.min(mins)
            max = np.max(maxs)
            # Generate a range of levels that will be consistent for both plots
            levels = np.linspace(min, max, N_levels)  # 10 contour levels
            pass
        else:
            pass

        data = []

        if plotting_method == "surface":
            pass
        elif plotting_method == "contour":
            dim1, dim2 = self.plot_dim
            fig, ax = plt.subplots(dim1, dim2, figsize=self.figuresize)
            # Choose colormap
            if realimag == "phase":
                plot_cmp = "gray"
                pass
            else:
                plot_cmp = "cividis"
                pass
            pass
        else:
            raise Exception("Invalid plotting method")
        
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
                state = abs(state)
                if plotting_method == "contour":
                    levels = abs(levels)
                    # Order levels as min to max
                    levels.sort()
                    pass
                else:
                    pass
                pass
            elif realimag == "phase":
                state = np.angle(state)
                pass
            else:
                raise Exception("Invalid value for realimag")

            if plotting_method == "surface":
                ax = plt.figure().add_subplot(projection='3d')
                ax.plot_trisurf(coords[:,0],coords[:,1], state, linewidth=0.2, antialiased=True, cmap="Spectral")
                # Adjust aspect ratio
                #min = np.min((cell.N, cell.N2))
                L1 = cell.N/self.Nmin
                L2 = cell.N2/self.Nmin
                ax.set_box_aspect([L2, L1, 1])  # (X, Y, Z) aspect ratio
                ax.set_zlabel(r'$\psi$'+str(indices)+"(x,y)",rotation=90)
                ax.ticklabel_format(axis='z',style='plain',useOffset=True)
                ax.zaxis.labelpad=10
                plt.title("D"+str(indices))
                pass
            elif plotting_method == "contour":
                triang = tri.Triangulation(coords[:, 0], coords[:, 1])
                #dim1, dim2 = self.plot_dim
                #fig, ax = plt.subplots(dim1, dim2)
                sb1, sb2 = cell.plot_loc
                if self.plot_dim == (1,1):
                    contour_plot = ax.tricontourf(triang, state, levels=levels, cmap=plot_cmp)
                    contour_lines = ax.tricontour(triang, state, levels=levels, colors='black', linewidths=0.4)
                    pass
                else:
                    contour_plot = ax[sb1,sb2].tricontourf(triang, state, levels=levels, cmap=plot_cmp)
                    contour_lines = ax[sb1,sb2].tricontour(triang, state, levels=levels, colors='black', linewidths=0.4)
                    pass
                #cbar = fig.colorbar(contour_plot, ax=ax)
                id1, id2 = indices
                strid = str(id1)+","+str(id2)
                title = r'$\psi_{{{}}}$'.format(strid)
                if self.plot_dim == (1,1):
                    ax.set_title(title, fontsize=12)
                    pass
                else:
                    ax[sb1,sb2].set_title(title, fontsize=20)
                    pass
                if type(cell) == cls.TriangleCell:
                    xs = np.linspace(coords[0,0],coords[-1,0],100)
                    if cell.diag.eliminated == True and cell.diag.strongly_eliminated == True:
                        ys = xs + 2*coords[0,0]
                        pass
                    elif cell.diag.eliminated == True:
                        ys = xs + coords[0,0]
                        pass
                    else:
                        ys = xs
                        pass 
                    if self.plot_dim == (1,1):
                        line = ax.plot(xs,ys,color='black',linewidth=0.5)
                        pass
                    else:
                        line = ax[sb1,sb2].plot(xs,ys,color='black',linewidth=0.5)
                        pass
                    pass
                pass
            else:
                raise Exception("Invalid plotting method")
            
            id1 , id2 = indices
            xl = fr'$x_{{e_{{{id1}}}}}$'
            yl = fr'$y_{{e_{{{id2}}}}}$'
            if plotting_method == "surface" or (plotting_method == "contour" and self.plot_dim == (1,1)):
                ax.set_xlabel(xl)
                ax.set_ylabel(yl)
                pass
            elif plotting_method == "contour":
                if cell.use_x_labels == True:
                    ax[sb1,sb2].set_xlabel(xl,fontsize=16)
                    # Make x ticks larger
                    ax[sb1,sb2].tick_params(axis='x', labelsize=12)
                    pass
                elif cell.use_x_labels == False:
                    ax[sb1,sb2].get_xaxis().set_visible(False)
                    pass
                else:
                    pass
                if cell.use_y_labels == True:
                    ax[sb1,sb2].set_ylabel(yl,fontsize=16)
                    # Make y ticks larger
                    ax[sb1,sb2].tick_params(axis='y', labelsize=12)
                    # Rotate y-axis label
                    ax[sb1,sb2].yaxis.label.set_rotation(0)
                    ax[sb1, sb2].yaxis.set_label_coords(-0.5, 0.5)
                    pass
                elif cell.use_y_labels == False:
                    ax[sb1,sb2].get_yaxis().set_visible(False)
                    pass
                else:
                    pass
                pass
            else:
                raise Exception("Invalid plotting method")

            if show_plots == True and plotting_method == "surface":
                plt.show()
                pass
            else:
                pass

            #plt.clf()
            #plt.cla()
            #plt.close()
            
            c+=length

            if return_data == True:
                data.append((coords[:,0],coords[:,1],state))
                pass
            else:
                pass
            pass

        if show_plots == True:

            if plotting_method == "surface":
                #plt.show()
                pass
            elif plotting_method == "contour":
                if self.plot_dim == (1,1):
                    cbar = fig.colorbar(contour_plot, ax=ax)
                    plt.show()
                    pass
                else:
                    for po in self.plots_off:
                        po1, po2 = po
                        ax[po1,po2].axis('off')
                        pass
                    #fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.3, hspace=0.4)
                    cbar = fig.colorbar(contour_plot, ax=ax)
                    cbar.ax.tick_params(labelsize=12)
                    #plt.savefig("foo.pdf", bbox_inches="tight")
                    plt.show()
                    pass
                pass
            else:
                raise Exception("Invalid plotting method")
            pass
        else:
            pass

        plt.close('all')
        plt.clf()
        
        if return_data == True:
            return data
        else:
            pass

        return None
    
    def save_states(self, filename):
        # Save the states of the system to a file

        #directory_path = "C:/Users/vq22287/OneDrive - University of Bristol/Documents/Anyons on graphs/Eigenstates data/"

        #filename = directory_path + filename

        cells = self.cells

        c = 0
        for cell in cells:
            self.construct_coords(cell)

            coords = cell.non_elim_coords
            length = cell.num_non_elim

            states = cell.eigenstates
            #states = np.real(states)

            #data = np.array((coords[:,0],coords[:,1],states[:])).T
            data = np.hstack((coords, states))

            indices = cell.indices

            #np.savetxt(filename+"D"+str(indices)+".csv",data,delimiter=",")
            np.savetxt(filename+"D"+str(indices)+".csv",data.view(float),delimiter=",")
            c+=length
            pass

        return None
    
    def load_states(self, filename,override_directory_path=False):
        # Load the states of the system from a file

        if override_directory_path == True:
            directory_path = ""
            pass
        else:
            directory_path = "C:/Users/vq22287/OneDrive - University of Bristol/Documents/Anyons on graphs/Eigenstates data/"
            pass

        filename = directory_path + filename

        cells = self.cells

        # Set the appropriate flags
        if self.solved == False:
            if self.diagonal_boundary_condition == "dirichlet":
                diag = self.diagonal
                for edge in diag:
                    edge.eliminated = True
                    pass
                pass
            elif self.diagonal_boundary_condition == "neumann" or self.diagonal_boundary_condition == "robin":
                pass
            else:
                raise Exception("Unknown diagonal boundary condition")
            
            if self.exterior_boundary_condition == "dirichlet":
                free_edges = self.free_edges
                for edge in free_edges:
                    edge.eliminated = True
                    pass
                pass
            elif self.exterior_boundary_condition == "neumann" or self.exterior_boundary_condition == "robin":
                pass
            else:
                raise Exception("Unknown exterior boundary condition")
            
            for g in self.gluings:
                for edge in g:
                    edge.eliminated = True
                    pass
                pass

            for g in self.gluings_with_branch_cut:
                for edge in g[0]:
                    edge.eliminated = True
                    pass
                pass

            self.solved = True
            pass
        else:
            pass

        c = 0
        for cell in cells:
            self.construct_coords(cell)

            coords = cell.non_elim_coords
            length = cell.num_non_elim

            #data = np.loadtxt(filename+"D"+str(cell.indices)+".csv",delimiter=",")
            data = np.loadtxt(filename+"D"+str(cell.indices)+".csv",delimiter=",").view(complex)
            states = data[:,2:]

            cell.eigenstates = states

            c+=length
            pass

        return None
    
    def save_eigenvalues(self, filename):
        # Save the eigenvalues of the system to a file

        #directory_path = "C:/Users/vq22287/OneDrive - University of Bristol/Documents/Anyons on graphs/Eigenstates data/"

        #filename = directory_path + filename

        spec = self.spectrum

        np.savetxt(filename+".csv",spec,delimiter=",")

        return None
    
    def load_eigenvalues(self, filename, override_directory_path=False):
        # Load the eigenvalues of the system from a file

        if override_directory_path == True:
            directory_path = ""
            pass
        else:
            directory_path = "C:/Users/vq22287/OneDrive - University of Bristol/Documents/Anyons on graphs/Eigenstates data/"
            pass
        
        filename = directory_path + filename

        spec = np.loadtxt(filename+".csv",delimiter=",")

        self.spectrum = spec

        return None



if __name__ == "__main__":
    # Main

    N = 50

    h = (np.pi)/(N-1)
    
    D11 = cls.TriangleCell(N)
    D12 = cls.SquareCell(N)
    D22 = cls.TriangleCell(N)

    D11.indices = (1,1)
    D12.indices = (1,2)
    D22.indices = (2,2)

    C = ConfigurationSpace([D11,D12,D22])

    gluing1 = [D11.y1,D12.y0]
    gluing2 = [D12.x1,D22.x0]
    C.glue(gluing1)
    C.glue(gluing2)

    #C.exterior_boundary_condition("dirichlet")
    #C.diagonal_boundary_condition("dirichlet")
    C.gen_lapl()
    lapl = C.sL
    #print(lapl.toarray())
    #C.print_eqs()

    #spec = C.lapl_spectrum(h,2,10)
    #print(spec)

    #C.plot_states(0,10)
    #spec = C.spectrum
    #spec.sort()
    #print(spec)



    pass