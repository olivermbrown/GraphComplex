"""
Code to implement the finite difference method on a square
"""

import numpy as np
import scipy as sp
import networkx as nx

class SquareCell:

    def __init__(self, length):
        """
        Initialize the SquareCell object with the given parameters.

        Args:
            self: The SquareCell object.
            length: The side length of the square domain.
        """

        # Creae a grid graph of the given dimensions.
        self.N = length
        self.G = nx.grid_graph(dim=(length,length))

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

        return None

    def find_corners(self):
        """
        Find the nodes in the graph corresponding to the corners of the domain.

        Args:
            self: The SquareCell object.
        
        Returns:
            corners: A list of the nodes corresponding to the corners of the domain.
        """

        N = self.N
        nodes_list = list(self.G.nodes())
        corner_indices = [0,N-1,(N**2)-N,(N**2)-1]
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
        
        N = self.N
        nodes_list = list(self.G.nodes())
        
        # Find the nodes on the first y-axis
        x0_nodes = []
        for x in range(1,N-1):
            x0_nodes.append(nodes_list[x])
            pass
        
        # Find the nodes on the first x-axis
        y0_nodes = []
        for x in range(N,(N**2)-N,N):
            y0_nodes.append(nodes_list[x])
            pass
        
        # Find the nodes on the opposite y-axis
        x1_nodes = []
        for x in range((N**2)-N+1,(N**2)-1):
            x1_nodes.append(nodes_list[x])
            pass
        
        # Find the nodes on the opposite x-axis
        y1_nodes = []
        for x in range(2*N-1,N**2-1,N):
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

    def __init__(self, length):
        """
        Initialize the TriangleCell object with the given parameters.

        Args:
            self: The TriangleCell object.
            length: The side length of the square domain.
        """

        # Create a grid graph of the given dimensions.

        super().__init__(length)

        # Find the diagonal of the domain.
        diag = self.find_diag()
        self.edges.append(diag)

        # Delete nodes beyond the diagonal from the domain in order to form a triangular cell.
        self.G = self.form_triangle()

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
        N = self.N
        nodes_list = list(self.G.nodes())

        # Find the nodes along the diagonal x=y.
        diag_nodes = []
        for x in range(N+1,(N**2 - 1), N + 1):
            diag_nodes.append(nodes_list[x])
            pass

        # Create an edge object for the diagonal.
        diag = Edge(self,diag_nodes)

        # Create a class variable for the diagonal.
        self.diag = diag

        return diag

    def form_triangle(self):
        """
        Form a triangular domain from the square domain.

        Args:
            self: The TriangleCell object.
        """

        N = self.N
        G = self.G

        # Find the nodes in the graph corresponding to the corners of the domain.
        corners = self.corners

        # Calculate Laplacian for Graph G
        L = nx.laplacian_matrix(G)

        # Get list of nodes
        nodes_list = list(G.nodes())

        # Create list of indices to drop from the Laplacian matrix.
        drop_indices = []
        for i in range (N):
            drop_indices += [d for d in range(i*N+i+1,(i+1)*N)]
            pass

        # Identify the nodes to drop from the graph.
        drop_nodes = [nodes_list[i] for i in drop_indices]

        # Drop excess nodes from the grid graph.
        G.remove_nodes_from(drop_nodes)

        # Record deleted nodes.
        self.dropped_indices = drop_indices

        # Delete corners and extra nodes from the objects's data.
        self.corners.remove(nodes_list[N-1])

        self.edges.remove(self.x0)
        self.x0 = None
        self.edges.remove(self.y1)
        self.y1 = None
        

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
    
    N = 5 # The side length of the square
    
    # Scaling factor
    h = (np.pi)/(N-1)
    
    D = TriangleCell(N)
    print(D.G)
    print(D.G.nodes())
    print(D.corners)
    for e in D.edges:
        print(e.edge)
        pass

    D.lapl()
    print(D.L.toarray())
