"""
Code to implement the finite difference method on a square
"""

import numpy as np
import scipy as sp
import networkx as nx

class SquareCell:

    def __init__(self, length1, length2=None):
        """
        Initialize the SquareCell object with the given parameters.

        Args:
            self: The SquareCell object.
            length: The side length of the square domain.
        """

        # Creae a grid graph of the given dimensions.
        self.N = length1
        if length2 is None:
            length2 = length1
            pass
        else:
            pass
        self.N2 = length2
        self.G = nx.grid_graph(dim=(length1,length2))

        # Find the corners and edges of the domain.
        self.corners = self.find_corners()
        self.edges = self.find_edges()

        # Initialize the remaining attributes.
        self.removed_nodes = []
        self.non_elim_coords = []
        self.non_elim_nodes = []
        self.num_non_elim = None
        self.indices = None
        self.eigenstates = None
        self.cell_coords_constructed = False
        self.plot_loc = (0,0)
        self.use_x_labels = True
        self.use_y_labels = True

        return None

    def find_corners(self):
        """
        Find the nodes in the graph corresponding to the corners of the domain.

        Args:
            self: The SquareCell object.
        
        Returns:
            corners: A list of the nodes corresponding to the corners of the domain.
        """

        N1 = self.N
        N2 = self.N2
        nodes_list = list(self.G.nodes())
        corner_indices = [0,N1-1,(N1*N2)-N1,(N1*N2)-1]
        corner_indices.sort()

        corners = []

        for i in corner_indices:
            corners.append(nodes_list[i])
            pass

        return corners

    def find_edges(self):
        """
        Find the nodes corresponding to each edge of the domain.

        Args:
            self: The SquareCell object.
        
        Returns:
            edges: A list of the nodes corresponding to each edge of the domain.
        """
        
        N1 = self.N
        N2 = self.N2
        nodes_list = list(self.G.nodes())
        
        # Find the nodes on the first y-axis
        x0_nodes = []
        for x in range(1,N1-1):
            x0_nodes.append(nodes_list[x])
            pass
        
        # Find the nodes on the first x-axis
        y0_nodes = []
        for x in range(N1,(N1*N2)-N1,N1):
            y0_nodes.append(nodes_list[x])
            pass
        
        # Find the nodes on the opposite y-axis
        x1_nodes = []
        for x in range((N1*N2)-N1+1,(N1*N2)-1):
            x1_nodes.append(nodes_list[x])
            pass
        
        # Find the nodes on the opposite x-axis
        y1_nodes = []
        for x in range(2*N1-1,N1*N2-1,N1):
            y1_nodes.append(nodes_list[x])
            pass

        # Create edge objects
        x0 = Edge(self,x0_nodes)
        y0 = Edge(self,y0_nodes)
        x1 = Edge(self,x1_nodes)
        y1 = Edge(self,y1_nodes)
        
        # Append the edges to a list
        edges = []
        edges.append(x0)
        edges.append(y0)
        edges.append(x1)
        edges.append(y1)
        
        # Assign the edges as class variables
        self.x0 = x0
        self.y0 = y0
        self.x1 = x1
        self.y1 = y1
        
        return edges

    def lapl(self):
        """
        Generate the Laplacian for the square domain.

        Args:
            self: The SquareCell object.
        """
        
        N = self.N
        G = self.G
        corners = self.corners

        # Calculate Laplacian for Graph G
        L = nx.laplacian_matrix(G)
        
        # Get list of nodes
        nodes_list = list(G.nodes())
        
        removed_nodes = []
        # Delete corners from the Laplacian matrix
        
        self.removed_nodes = removed_nodes
        
        self.L = L
        
        return None

class TriangleCell(SquareCell):

    def __init__(self, length1, length2=None):
        """
        Initialize the TriangleCell object with the given parameters.

        Args:
            self: The TriangleCell object.
            length: The side length of the square domain.
        """

        # Create a grid graph of the given dimensions.

        super().__init__(length1,length2)

        # Find the diagonal of the domain.
        diag, diag_ids = self.find_diag()
        self.edges.append(diag)

        # Delete nodes beyond the diagonal from the domain in order to form a triangular cell.
        self.G = self.form_triangle(diag_ids)

        return None

    def find_diag(self):
        """
        Find the nodes in the graph corresponding to the diagonal of the domain.

        Args:
            self: The TriangleCell object.
        
        Returns:
            diag: The nodes corresponding to the diagonal of the domain.
        """

        # Initialize variables.
        N1 = self.N
        N2 = self.N2
        nodes_list = list(self.G.nodes())

        if N1 % 2 == 0 and N2 % 2 == 0:
            pass
        else:
            raise ValueError("The side lengths of the square must be even.")

        # Find the nodes along the diagonal x=y.
        diag_nodes = []
        # Get the x-y coordinates of the diagonal nodes        

        diag_ids = []
        for x in np.linspace(0,(N1*N2 - 1), N2):#range(0,(N1*N2 - 1), N2 - 1):
            if x == 0:
                continue
            elif x == (N1*N2 - 1):
                continue
            else:
                pass
            # Round x to the nearest integer
            idx = int(round(x))
            diag_ids.append(idx)
            diag_nodes.append(nodes_list[idx])
            pass

        # Create an edge object for the diagonal.
        diag = Edge(self,diag_nodes)

        # Create a class variable for the diagonal.
        self.diag = diag

        return diag, diag_ids

    def form_triangle(self,diag_ids):
        """
        Form a triangular domain from the square domain.

        Args:
            self: The TriangleCell object.
        """

        N1 = self.N
        N2 = self.N2
        G = self.G

        # Find the nodes in the graph corresponding to the corners of the domain.
        corners = self.corners

        # Calculate Laplacian for Graph G
        L = nx.laplacian_matrix(G)

        # Get list of nodes
        nodes_list = list(G.nodes())

        # Create list of indices to drop from the Laplacian matrix.
        drop_indices = []
        diag_ids.append(N1*N2 - 1)
        for i in range (1,N2):
            diag_idx = diag_ids[i-1]
            drop_indices += [d for d in range(i*N1,diag_idx)]
            pass

        # Identify the nodes to drop from the graph.
        drop_nodes = [nodes_list[i] for i in drop_indices]

        # Drop excess nodes from the grid graph.
        G.remove_nodes_from(drop_nodes)

        # Record deleted nodes.
        self.dropped_indices = drop_indices

        # Delete corners and extra nodes from the objects's data.
        self.corners.remove(nodes_list[N1*(N2-1)])

        self.edges.remove(self.y0)
        self.y0 = None
        self.edges.remove(self.x1)
        self.x1 = None
        

        return G
    
class Edge:
    """
    A class containing data about the edges of a two-particle domain.
    """

    def __init__(self, domain, edge_nodes):
        """
        Initialize the Edge object with the given parameters.

        Args:
            self: The Edge object.
            domain: The Domain object to which the edge belongs.
            edge_nodes: The nodes corresponding to the edge of the domain.
        """

        self.domain = domain
        self.edge = edge_nodes
        self.eliminated = False
        self.strongly_eliminated = False

        return None

    def find_edge_indices_in_configuration_space(self, cfgspace):
        """
        Find the indices of the nodes corresponding to the edge in the configuration space.

        Args:
            self: The Edge object.
            cfgspace: The configuration space of the system.
        """

        cells = cfgspace.cells

        edge = self.edge
        domain = self.domain

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
        
        self.edge_indices = edge_coords

        return None



if __name__=="__main__":
    # Main
    
    N = 6 # The side length of the square
    
    # Scaling factor
    h = (np.pi)/(N-1)
    
    D = TriangleCell(N,N-2)
    print(D.G)
    print(D.G.nodes())
    print(D.corners)
    for e in D.edges:
        print(e.edge)
        pass

    D.lapl()
    print(D.L.toarray())
