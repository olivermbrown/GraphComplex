import networkx as nx
import matplotlib.pyplot as plt

def create_hexagonal_lattice():
    G = nx.Graph()
    
    # Define the positions of the 6 sites
    positions = {
        0: (0, 0),
        1: (1, 0),
        2: (1.5, 0.866),
        3: (1, 1.732),
        4: (0, 1.732),
        5: (-0.5, 0.866)
    }
    
    # Add nodes with their positions
    for i in range(6):
        G.add_node(i, pos=positions[i])
    
    # Add edges to form the hexagon
    edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
    G.add_edges_from(edges)
    
    # Add periodic boundary conditions (opposite vertices)
    G.add_edge(0, 3)
    G.add_edge(1, 4)
    G.add_edge(2, 5)
    
    return G

# Create the lattice
G = create_hexagonal_lattice()

# Draw the graph
pos = nx.get_node_attributes(G, 'pos')
plt.figure(figsize=(8, 6))
nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=500, font_size=12, font_weight='bold')
edge_labels = {(u, v): '' for (u, v) in G.edges()}
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
plt.title("6-site Hexagonal Lattice with Periodic Boundary Conditions")
plt.axis('off')
plt.tight_layout()
plt.show()
# Calculate adjacency matrix of G
A = nx.adjacency_matrix(G).toarray()
print(A)
