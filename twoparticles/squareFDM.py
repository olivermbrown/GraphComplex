"""
Code to implement the finite difference method on a square
"""

import numpy as np
import scipy as sp
import networkx as nx



class Domain:
    # A class for discretisation of square domains, which can be glued together
    
    def __init__(self, length):
        self.N = length
        self.G = nx.grid_graph(dim=(length,length))
        self.split = False
        self.find_corners()
        self.edges = []
        self.find_edges()
        self.removed_nodes = []
        return None
    
    def find_corners(self):
        # Find the nodes in the graph of the corners of the domain
        N = self.N
        nodes_list = list(self.G.nodes())
        corner_indices = [0,N-1,(N**2)-N,(N**2)-1]
        corner_indices.sort()
        self.corners = []
        for i in corner_indices:
            self.corners.append(nodes_list[i])
        return None
    
    def find_edges(self):
        # Find the nodes corresponding to each edge of the domain
        
        N = self.N
        nodes_list = list(self.G.nodes())
        
        # Find the nodes on the first x-axis
        x0 = []
        for x in range(1,N-1):
            x0.append(nodes_list[x])
            pass
        
        # Find the nodes on the first y-axis
        y0 = []
        for x in range(N,(N**2)-N,N):
            y0.append(nodes_list[x])
            pass
        
        # Find the nodes on the opposite x-axis
        x1 = []
        for x in range((N**2)-N+1,(N**2)-1):
            x1.append(nodes_list[x])
            pass
        
        # Find the nodes on the opposite y-axis
        y1 = []
        for x in range(2*N-1,N**2-1,N):
            y1.append(nodes_list[x])
            pass
        
        # Find the nodes along the diagonal x=y
        diag = []
        for x in range(N+1,(N**2 - 1), N + 1):
            diag.append(nodes_list[x])
            pass
        
        x0inv = x0[::-1]
        y0inv = y0[::-1]
        x1inv = x1[::-1]
        y1inv = y1[::-1]
        
        self.edges.append(x0)
        self.edges.append(y0)
        self.edges.append(x1)
        self.edges.append(y1)
        self.edges.append(diag)
        
        # Assign the lists as class variables
        self.x0 = x0
        self.y0 = y0
        self.x1 = x1
        self.y1 = y1
        self.diag = diag
        self.x0inv = x0inv
        self.y0inv = y0inv
        self.x1inv = x1inv
        self.y1inv = y1inv
        
        return None
    
    def split_domain(self,condition):
        
        self.split = True
        self.split_condition = condition
        
        return None
    
    def lapl(self):
        # Generate the Laplacian for the square domain
        
        N = self.N
        G = self.G
        corners = self.corners

        # Calculate Laplacian for Graph G
        L = nx.laplacian_matrix(G)
        
        # Get list of nodes
        nodes_list = list(G.nodes())
        
        removed_nodes = []
        # Delete corners from the Laplacian matrix
        """
        for n in corners:
            i = nodes_list.index(n)
            c = N**2 - L.shape[0]
            i -= c
            L = sp.sparse.vstack([L[:i,:],L[i+1:,:]],format='csr')
            L = sp.sparse.hstack([L[:,:i],L[:,i+1:]],format='csr')
            removed_nodes.append(n)
            pass
        """
        
        self.removed_nodes = removed_nodes
        
        self.L = L
        
        return None
    
    def apply_dirichlet(self,edge):
        # Apply Dirichlet boundary conditions to a particular edge
        
        N = self.N
        G = self.G
        L = self.L
        
        # Identify which edge has been selected
        B = edge
        
        # Get list of nodes
        nodes_list = list(G.nodes())
        
        removed_nodes = self.removed_nodes
        # Apply Dirichlet boundary condition to this edge in the Laplacian matrix
        for n in B:
            i = nodes_list.index(n)
            c=0
            for r in removed_nodes:
                if nodes_list.index(r) > i:
                    pass
                else:
                    c += 1
                pass
            i -= c
            L = sp.sparse.vstack([L[:i,:],L[i+1:,:]],format='csr')
            L = sp.sparse.hstack([L[:,:i],L[:,i+1:]],format='csr')
            removed_nodes.append(n)
            pass
        
        self.removed_nodes = removed_nodes
        
        self.L = L
        
        return None
    
    def apply_neumann(self,edge):
        # Apply Neumann boundary conditions to a particular edge
        
        N = self.N
        G = self.G
        L = self.L
        
        # Identify which edge has been selected
        B = edge
        
        # Get list of nodes
        nodes_list = list(G.nodes())
        
        removed_nodes = self.removed_nodes
        
        # Apply Neumann boundary condition to this edge in the Laplacian matrix
        
        for n in B:
            i = nodes_list.index(n)
            
            if edge == self.x0:
                j = i + N
                k = i + 2*N
                pass
            if edge == self.y0:
                j = i + 1
                k = i + 2
                pass
            if edge == self.x1:
                j = i - N
                k = i - 2*N
                pass
            if edge == self.y1:
                j = i - 1
                k = i - 2
                pass
            pass
        
            c=0
            d=0
            e=0
            
            for r in removed_nodes:
                if nodes_list.index(r) > i:
                    pass
                else:
                    c += 1
                    pass
                
                if nodes_list.index(r) > j:
                    pass
                else:
                    d += 1
                    pass
                
                if nodes_list.index(r) > k:
                    pass
                else:
                    e += 1
                    pass
                pass
            
            i -= c
            j -= d
            k -= e
            
            # First order finite difference
            L[j,j] -= 1
            
            # Second order finite difference
            #L = L.astype('float32')
            #L[j,j] -= 4/3
            #L[j,k] += 1/3
            
            L = sp.sparse.vstack([L[:i,:],L[i+1:,:]],format='csr')
            L = sp.sparse.hstack([L[:,:i],L[:,i+1:]],format='csr')
            removed_nodes.append(n)
            
            
        
        self.removed_nodes = removed_nodes
        
        self.L = L
        
        return None
    
    def apply_robin(self,edge,a,b):
        # Apply Robin boundary conditions to a particular edge
        
        N = self.N
        G = self.G
        L = self.L
        
        # Scaling factor
        h = (np.pi)/(N-1)
        
        # Quotient defining the Robin boundary conditions
        K = (a*h)/b
        
        # Identify which edge has been selected
        B = edge
        
        # Get list of nodes
        nodes_list = list(G.nodes())
        
        removed_nodes = self.removed_nodes
        
        # Apply Robin boundary condition to this edge in the Laplacian matrix
        
        for n in B:
            i = nodes_list.index(n)
            
            if edge == self.x0:
                j = i + N
                k = i + 2*N
                pass
            if edge == self.y0:
                j = i + 1
                k = i + 2
                pass
            if edge == self.x1:
                j = i - N
                k = i - 2*N
                pass
            if edge == self.y1:
                j = i - 1
                k = i - 2
                pass
            pass
        
            c=0
            d=0
            e=0
            
            for r in removed_nodes:
                if nodes_list.index(r) > i:
                    pass
                else:
                    c += 1
                    pass
                
                if nodes_list.index(r) > j:
                    pass
                else:
                    d += 1
                    pass
                
                if nodes_list.index(r) > k:
                    pass
                else:
                    e += 1
                    pass
                pass
            
            i -= c
            j -= d
            k -= e
            
            L = L.astype('float32')
            
            # First order finite difference
            s = (-1)/(K - 1)
            L[j,j] -= s
            
            # Second order finite difference
            #s = (-2)/(K - 3/2)
            #t = 1/(2*K - 3)
            #L[j,j] -= s
            #L[j,k] -= t
            
            L = sp.sparse.vstack([L[:i,:],L[i+1:,:]],format='csr')
            L = sp.sparse.hstack([L[:,:i],L[:,i+1:]],format='csr')
            removed_nodes.append(n)
            
            
        
        self.removed_nodes = removed_nodes
        
        self.L = L
        
        return None
    
class Edge:
    # A class for the edge of a two-particle domain

    def __init__(self, domain, edge_nodes):
        self.domain = domain
        self.edge = edge_nodes

        return None


def lapl_spectrum(matrix,h,dps):
    # Calculate the spectrum of a matrix M
    # Scaling factor h
    # Decimal places to round to is dps
    
    matrix = matrix.astype('float32')
    eigs = sp.sparse.linalg.eigs(matrix,which='SM',k=10,return_eigenvectors=False)
    
    spec = (h)**(-2) * np.abs(eigs)
    spec.sort()
    return np.round(spec,dps)

if __name__=="__main__":
    # Main
    
    N = 1000 # The side length of the square
    
    # Scaling factor
    h = (np.pi)/(N-1)
    
    D = Domain(N)
    G = D.G
    print(G)
    
    D.lapl()
    D.apply_neumann(D.x0)
    D.apply_neumann(D.y0)
    D.apply_neumann(D.x1)
    D.apply_neumann(D.y1)
    print(D.L)
    print(D.L.toarray())
    
    
    # Calculate the spectrum of the Laplacian
    
    spectrum = lapl_spectrum(D.L,h,2)
    print(spectrum)










"""
def dirichlet_lapl(G):
    # Calculate the Laplacian matrix with Dirichlet boundary conditions for a grid G
    
    # Find the size of the graph
    N = int(np.sqrt(nx.number_of_nodes(G)))
    
    # Calculate Laplacian for Graph G
    L = nx.laplacian_matrix(G)
    
    # Generate list of boundary "end nodes"
    nodes_list = list(G.nodes())
    
    # Generate a list of the indices of nodes which are on the boundary of the square
    end_indices = []
    for x in range(0,N):
        end_indices.append(x)
    for x in range(N,(N**2)-N,N):
        end_indices.append(x)
    for x in range(2*N-1,N**2-1,N):
        end_indices.append(x)
    for x in range((N**2)-N,N**2):
        end_indices.append(x)
    
    end_indices.sort()
    end_nodes = []
    for i in end_indices:
        end_nodes.append(nodes_list[i])
    
    #end_nodes = [x for x in G.nodes() if G.degree(x)==2]
    
    # Modify Laplacian for Dirichlet boundary conditions at each "end node"
    c = 0
    for n in end_nodes:
        i = nodes_list.index(n)
        #L[i,:] = 0
        #L[:,i] = 0
        c = N**2 - L.shape[0]
        i -= c
        L = sp.sparse.vstack([L[:i,:],L[i+1:,:]],format='csr')
        L = sp.sparse.hstack([L[:,:i],L[:,i+1:]],format='csr')
        #c += 1
        
    c = 0   
    
    L.eliminate_zeros()
    
    return L
"""

"""
def neumann_lapl(G):
    # Calculate the Laplacian matrix with Neumann boundary conditions for a grid G
    
    # Find the size of the graph
    N = int(np.sqrt(nx.number_of_nodes(G)))
    
    # Calculate Laplacian for Graph G
    L = nx.laplacian_matrix(G)
    
    # Generate list of boundary "end nodes"
    nodes_list = list(G.nodes())
    
    # Generate a list of the indices of nodes which are on the boundary of the square
    end_indices = []
    for x in range(0,N):
        end_indices.append(x)
    for x in range(N,(N**2)-N,N):
        end_indices.append(x)
    for x in range(2*N-1,N**2-1,N):
        end_indices.append(x)
    for x in range((N**2)-N,N**2):
        end_indices.append(x)
    
    end_indices.sort()
    end_nodes = []
    for i in end_indices:
        end_nodes.append(nodes_list[i])
    
    #end_nodes = [x for x in G.nodes() if G.degree(x)==2]
    
    # Modify Laplacian for Neumann boundary conditions at each "end node"
    # Begin by deleting end nodes from the matrix
    c = 0
    for n in end_nodes:
        i = nodes_list.index(n)
        #L[i,:] = 0
        #L[:,i] = 0
        i -= c
        L = sp.sparse.vstack([L[:i,:],L[i+1:,:]],format='csr')
        L = sp.sparse.hstack([L[:,:i],L[:,i+1:]],format='csr')
        c += 1
    c = 0
    
    # Modifies rows of the L matrix according to the central difference approximation
    
    L.eliminate_zeros()
    
    return L
"""