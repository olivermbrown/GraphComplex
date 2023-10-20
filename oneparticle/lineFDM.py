"""
Code to implement the finite difference method on a line
"""

import numpy as np
import scipy as sp
import networkx as nx

import matplotlib.pyplot as plt
from matplotlib import cm

class Line:
    # A class for discretisation of linear domains, which can be glued together
    
    def __init__(self, length):
        self.N = length
        self.G = nx.grid_graph(dim=(length,))
        self.find_endpoints()
        return None
    
    def find_endpoints(self):
        # Find the nodes in the graph at the endpoints of the line
        N = self.N
        nodes_list = list(self.G.nodes())
        
        endpoints = []
        self.start = End(self,nodes_list[0])#nodes_list[0]
        endpoints.append(self.start)
        self.end = End(self,nodes_list[N-1])#nodes_list[N-1]
        endpoints.append(self.end)
        self.endpoints = endpoints
        return None
    
    def lapl(self):
        # Generate the Laplacian matrix for the line
        
        N = self.N
        G = self.G
        ends = self.endpoints

        # Calculate Laplacian for Graph G
        L = nx.laplacian_matrix(G)
        removed_nodes = []
        """
        # Store the original dimension of the Laplacian matrix
        Dim = L.shape[0]
        
        # Get list of nodes
        nodes_list = list(G.nodes())
        
        
        # Delete corners from the Laplacian matrix
        for n in ends:
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
    
    def apply_dirichlet(self, end):
        # Apply Dirichlet boundary conditions to an end of the line
        
        N = self.N
        G = self.G
        L = self.L
        
        # Get list of nodes
        nodes_list = list(G.nodes())
        removed_nodes = self.removed_nodes
        
        # Apply Dirichlet boundary condition to this edge in the Laplacian matrix

        i = nodes_list.index(end)
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
        removed_nodes.append(end)
        
        self.removed_nodes = removed_nodes
        
        self.L = L
        
        return None
    
    def apply_neumann(self,end):
        # Apply Dirichlet boundary condition to this edge in the Laplacian matrix
        
        N = self.N
        G = self.G
        L = self.L
        
        # Get list of nodes
        nodes_list = list(G.nodes())
        removed_nodes = self.removed_nodes
        
        # Modify Laplacian for Neumann boundary conditions at each "end node"
        i = nodes_list.index(end)
        
        if end == 0:
            j = i+1
            k = i+2
            pass
        if end == N-1:
            j = i-1
            k = i-2
            pass
        else:
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
        
        """
        if i == 0:
            L[i+1,i+2] -= 1
            pass
        else:
            L[i-1,i-2] -= 1
            pass
        """
        L = sp.sparse.vstack([L[:i,:],L[i+1:,:]],format='csr')
        L = sp.sparse.hstack([L[:,:i],L[:,i+1:]],format='csr')
        removed_nodes.append(end)
        
        self.removed_nodes = removed_nodes
        
        self.L = L
        
        return None
    
    def apply_robin(self,end,a,b,h):
        # Apply Robin boundary condition to this edge in the Laplacian matrix
        
        N = self.N
        G = self.G
        L = self.L
        
        L = L.astype('float32')
        
        # Quotient defining the Robin boundary conditions
        K = (2*a*h)/b
        
        # Get list of nodes
        nodes_list = list(G.nodes())
        removed_nodes = self.removed_nodes
        
        # Modify Laplacian for Neumann boundary conditions at each "end node"
        i = nodes_list.index(end)
        
        if end == 0:
            j = i+1
            k = i+2
            pass
        if end == N-1:
            j = i-1
            k = i-2
            pass
        else:
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
        s = (-1)/(K - 1)
        L[j,j] -= s
        
        """
        if i == 0:
            L[i+1,i+2] -= 1
            pass
        else:
            L[i-1,i-2] -= 1
            pass
        """
        L = sp.sparse.vstack([L[:i,:],L[i+1:,:]],format='csr')
        L = sp.sparse.hstack([L[:,:i],L[:,i+1:]],format='csr')
        removed_nodes.append(end)
        
        self.removed_nodes = removed_nodes
        
        self.L = L
        
        return None

def lapl_spectrum(matrix,h,dps):
    # Calculate the spectrum of a matrix M
    # Scaling factor h
    # Decimal places to round to is dps
    
    matrix = matrix.astype('float32')
    eigs = sp.sparse.linalg.eigs(matrix,which='SM',k=8,return_eigenvectors=False)
    
    spec = (h)**(-1) * np.sqrt(np.abs(eigs))
    spec.sort()
    return np.round(spec,dps)

def lapl_solve(matrix,h,dps):
    # Calculate the spectrum of a matrix M
    # Scaling factor h
    # Decimal places to round to is dps
    
    matrix = matrix.astype('float32')
    eigvals, eigvecs = sp.sparse.linalg.eigs(matrix,which='SM',k=8,return_eigenvectors=True)
    
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
            x = node/N
            #y = node[1]/N
            coords.append([x])
            pass
        else:
            pass
        pass
    
    coords = np.array(coords)
    
    
    # Define coordinates
    
    
    # Plot the surface
    xs = coords[:,0]
    #ys = coords[:,1]
    fig, ax = plt.subplots()
    state = np.real(state)
    ax.plot(xs, state)
    
    
    return plt

class End:
    # A class for the ends of a line

    def __init__(self, line, node):
        self.line = line
        self.end = node

        return None

if __name__=="__main__":
    # Main
    
    N = 100 # The length of the line-chain
    
    # Scaling factor
    h = (np.pi)/(N-1)
    
    #G = nx.grid_graph(dim=(N,))
    #print(G)
    
    #L = dirichlet_lapl(G)
    #L = neumann_lapl(G)
    #L = robin_lapl(G,1,1,h)
    D = Line(N)
    G = D.G
    print(G)
    
    D.lapl()
    D.apply_robin(D.start,1,1,h)
    D.apply_robin(D.end,1,1,h)
    L = D.L
    
    print(L)
    print(L.toarray())
    
    # Calculate the spectrum of the Laplacian
    
    #spectrum = lapl_spectrum(L,h,2)
    #print(spectrum)
    
    spectrum, states = lapl_solve(D.L,h,2)
    print(spectrum)
    plot = plot_state(D, states[:,4])
    plot.show()



"""
def dirichlet_lapl(G):
    # Calculate the Laplacian matrix with Dirichlet boundary conditions for a grid G
    
    # Calculate Laplacian for Graph G
    L = nx.laplacian_matrix(G)
    
    # Generate list of boundary "end nodes"
    nodes_list = list(G.nodes())
    end_nodes = [x for x in G.nodes() if G.degree(x)==1]
    
    # Modify Laplacian for Dirichlet boundary conditions at each "end node"
    for n in end_nodes:
        i = nodes_list.index(n)
        L[i,:] = 0
        L[:,i] = 0
    
    L.eliminate_zeros()
    
    return L

def neumann_lapl(G):
    # Calculate the Laplacian matrix with Neumann boundary conditions for a grid G
    
    # Calculate Laplacian for Graph G
    L = nx.laplacian_matrix(G)
    
    # Generate list of boundary "end nodes"
    nodes_list = list(G.nodes())
    end_nodes = [x for x in G.nodes() if G.degree(x)==1]
    
    # Modify Laplacian for Neumann boundary conditions at each "end node"
    for n in end_nodes:
        i = nodes_list.index(n)
        L[i,:] = 0
        L[:,i] = 0
        if i == 0:
            L[i+1,i+2] -= 1
        else:
            L[i-1,i-2] -= 1
    
    L.eliminate_zeros()
    
    return L

def robin_lapl(G,a,b,h):
    # Calculate the Laplacian matrix with Neumann boundary conditions for a grid G
    
    # Quotient defining the Robin boundary conditions
    K = (2*a*h)/b
    
    # Calculate Laplacian for Graph G
    L = nx.laplacian_matrix(G)
    L = L.astype('float32')
    
    # Generate list of boundary "end nodes"
    nodes_list = list(G.nodes())
    end_nodes = [x for x in G.nodes() if G.degree(x)==1]
    
    # Modify Laplacian for Neumann boundary conditions at each "end node"
    for n in end_nodes:
        i = nodes_list.index(n)
        L[i,:] = 0
        L[:,i] = 0
        if i == 0:
            L[i+1,i+1] -= K
            L[i+1,i+2] -= 1
        else:
            L[i-1,i-1] += K
            L[i-1,i-2] -= 1
    
    L.eliminate_zeros()
    
    return L
"""