"""
Code to implement the finite difference method on a square
"""

import numpy as np
import scipy as sp
import networkx as nx

import matplotlib.pyplot as plt
from matplotlib import cm


class TriangleDomain:
    # A class for discretisation of square domains, which can be glued together
    
    def __init__(self, length):
        self.N = length
        self.G = nx.grid_graph(dim=(length,length))
        self.find_corners()
        self.edges = []
        self.find_edges()
        self.find_hypotenuse()
        self.remove_half()
        return None
    
    def find_corners(self):
        # Find the nodes in the graph of the corners of the domain
        N = self.N
        nodes_list = list(self.G.nodes())
        corner_indices = [0,N-1,(N**2)-N]
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
        """
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
        """
        # Assign the lists as class variables
        self.x0 = x0
        self.y0 = y0
        
        self.edges.append(x0)
        self.edges.append(y0)
        #self.x1 = x1
        #self.y1 = y1
        
        return None
    
    def find_hypotenuse(self):
        # Add the hyptenuse to the list of edges
        
        N = self.N
        nodes_list = list(self.G.nodes())
        
        # Find the nodes on the hypotenuse
        hyp = []
        for x in range (2*N-2,N**2-N,N-1):
            hyp.append(nodes_list[x])
            pass
        
        self.hyp = hyp
        
        self.edges.append(hyp)
        
        return None
    
    def remove_half(self):
        # Remove the second half of the triangle from the graph
        
        N = self.N
        nodes_list = list(self.G.nodes())
        
        nodes_to_remove = []
        
        i = 0
        for node in self.hyp:
            start = nodes_list.index(node) + 1
            end = nodes_list.index(self.y0[i]) + N
            for n in range (start, end):
                nodes_to_remove.append(nodes_list[n])
                pass
            i += 1
            pass
        
        for m in range((N**2)-N+1,(N**2)):
            nodes_to_remove.append(nodes_list[m])
            pass
        
        for rm in nodes_to_remove:
            self.G.remove_node(rm)
            pass
        
        return None
    
    def lapl(self):
        # Generate the Laplacian for the triangular domain
        
        N = self.N
        G = self.G
        corners = self.corners

        # Calculate Laplacian for Graph G
        L = nx.laplacian_matrix(G)
        
        # Store the original dimension of the Laplacian matrix
        Dim = L.shape[0]
        
        # Get list of nodes
        nodes_list = list(G.nodes())
        
        removed_nodes = []
        # Delete corners from the Laplacian matrix
        """
        for n in corners:
            i = nodes_list.index(n)
            c = Dim - L.shape[0]
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
        if edge == self.hyp:
            self.apply_dirichlet(edge)
            #self.apply_neumann_hyp(self,edge)
            pass
        else:
            self.apply_neumann_exterior(edge)
            pass
        return None
    
    def apply_neumann_hyp(self,hyp):
        # Apply Neumann boundary conditions to the hypotenuse #TODO
        
        N = self.N
        G = self.G
        L = self.L
        
        # Identify which edge has been selected
        #B = hyp
        
        # Get list of nodes
        nodes_list = list(G.nodes())
        
        removed_nodes = self.removed_nodes
        
        # Apply Neumann boundary condition to this edge in the Laplacian matrix
        
        for n in hyp:
            i = nodes_list.index(n)
            """
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
            """
            j = i - 1
            k = i - N
            
            
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
            
            # TODO
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
    
    def apply_neumann_exterior(self,edge):
        # Apply Neumann boundary conditions to an exterior edge
        
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
            """
            if edge == self.x1:
                j = i - N
                k = i - 2*N
                pass
            if edge == self.y1:
                j = i - 1
                k = i - 2
                pass
            """
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
        if edge == self.hyp:
            self.apply_dirichlet(edge)
            #self.apply_robin_hyp(self,edge)
            pass
        else:
            self.apply_robin_exterior(edge,a,b)
            pass
        return None
    
    def apply_robin_exterior(self,edge,a,b):
        # Apply Robin boundary conditions to an exterior edge
        
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

def lapl_spectrum(matrix,h,dps):
    # Calculate the spectrum of a matrix M
    # Scaling factor h
    # Decimal places to round to is dps
    
    matrix = matrix.astype('float32')
    eigs = sp.sparse.linalg.eigs(matrix,which='SM',k=10,return_eigenvectors=False)
    
    spec = (h)**(-2) * np.abs(eigs)
    spec.sort()
    return np.round(spec,dps)

def lapl_solve(matrix,h,dps):
    # Calculate the spectrum of a matrix M
    # Scaling factor h
    # Decimal places to round to is dps
    
    matrix = matrix.astype('float32')
    eigvals, eigvecs = sp.sparse.linalg.eigs(matrix,which='SM',k=10,return_eigenvectors=True)
    
    spec = (h)**(-2) * np.abs(eigvals)
    spec.sort()
    return np.round(spec,dps), eigvecs

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

if __name__=="__main__":
    # Main
    
    N = 150 # The side length of the square
    
    # Scaling factor
    h = (np.pi)/(N-1)
    
    D = TriangleDomain(N)
    G = D.G
    print(G)
    
    D.lapl()
    D.apply_neumann(D.x0)
    D.apply_neumann(D.y0)
    D.apply_dirichlet(D.hyp)
    #D.apply_neumann(D.y1)
    print(D.L)
    print(D.L.toarray())
    
    
    # Calculate the spectrum of the Laplacian
    
    #spectrum = lapl_spectrum(D.L,h,2)
    #print(spectrum)
    
    spectrum, states = lapl_solve(D.L,h,2)
    print(spectrum)
    plot = plot_state(D, states[:,0])
    plot.show()
    







