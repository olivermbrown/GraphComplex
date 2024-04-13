"""
Code for finding the configuration space of two particles on a network of wires
"""

import numpy as np
import scipy as sp
import sympy as sym
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import squareFDM
import twocomplex

class FundamentalDomain(squareFDM.Domain):

    def __init__(self, domain):

        self.split = domain.split

        self.indices = domain.indices[0:2]

        self.non_elim_nodes = domain.non_elim_nodes

        self.non_elim_coords = domain.non_elim_coords

        self.associated_cells = []

        return None
    
    def associate(self, cell):

        self.associated_cells.append(cell)

        return None

class ConfigurationSpace:
    # A class for the configuration space of two particles on a network of wires

    def __init__(self, complx):

        self.covering_cells = complx.cells

        self.find_fundamental_domains()

        self.identify_cells()

        return None
    
    def find_fundamental_domains(self):
        # Define domains in the configuration space and identify cells in the complex with these domains
        # TODO

        covering_cells = self.covering_cells

        fundamental_domains = []

        for cell in covering_cells:
            if cell.split == True:
                if cell.indices[2] == 0:
                    fundamental_domains.append(FundamentalDomain(cell))
                    pass
                else:
                    pass
            elif cell.split == False:
                if cell.indices[0] < cell.indices[1] and cell.indices[2] == 0:
                    fundamental_domains.append(FundamentalDomain(cell))
                    pass
                pass
            else:
                raise Exception
                pass
            pass

        self.domains = fundamental_domains

        return None
    
    def identify_cells(self):

        covering_cells = self.covering_cells

        domains = self.domains

        for domain in domains:
            for cell in covering_cells:
                if domain.indices[0:2] == cell.indices[0:2]:
                    domain.associate(cell)
                    pass
                elif domain.indices[0:2] == cell.indices[1::-1]:
                    domain.associate(cell)
                    pass
                else:
                    pass
                pass
            pass

        print("Cells identified with domains")

        return None
    
    def bosonic_projection(self):
        # Project eigenstates on the covering complex onto the bosonic subspace in the indistinguishable configuration space
        
        domains = self.domains

        for domain in domains:
            if domain.split == True:
                
                for cell in domain.associated_cells:

                    domain.bosonic_states = sum([cell.eigenstates for cell in domain.associated_cells])

                    pass
                
                non_elim_coords = domain.non_elim_coords

                # Need to try to speed this code up:
                    
                # Create a list of coordinates with the axes flipped
                col1 = non_elim_coords[:,0]
                col2 = non_elim_coords[:,1]
                rev = np.c_[col2,col1]

                paired_indices = []

                for c in non_elim_coords:
                    if c in rev:
                        index1 = non_elim_coords.tolist().index(c.tolist())
                        index2 = rev.tolist().index(c.tolist())
                        paired_indices.append([index1,index2])
                        pass
                    else:
                        pass
                    pass

                bosonic_states_transposed = np.transpose(domain.bosonic_states)

                for state in bosonic_states_transposed:
                    projected_state = state.copy()
                    for indices in paired_indices:
                        projected_state[indices[0]] = state[indices[0]] + state[indices[1]]
                        pass
                    bosonic_states_transposed[bosonic_states_transposed.tolist().index(state.tolist())] = projected_state
                    pass

                domain.bosonic_states = np.transpose(bosonic_states_transposed)
                
                pass
            elif domain.split == False:

                domain.bosonic_states = sum([cell.eigenstates for cell in domain.associated_cells])

                pass
            else:
                raise Exception
            pass

        return None
    
    def plot_projected_wavefunction(self,n):
        # Plot the projected n-th wavefunction

        domains = self.domains

        for domain in domains:
            i = domains.index(domain)

            coords = domain.non_elim_coords

            state = domain.bosonic_states[:,n]
            state = np.real(state)

            ax = plt.figure().add_subplot(projection='3d')

            ax.plot_trisurf(coords[:,0],coords[:,1], state, linewidth=0.2, antialiased=True, cmap=cm.inferno)
            #ax.scatter(coords[:,0],coords[:,1],state,c=state,cmap=cm.autumn)
            indices = domain.indices

            plt.title("D"+str(indices))
            ax.set_xlabel("x_e"+str(indices[0]))
            ax.set_ylabel("y_e"+str(indices[1]))
            ax.zaxis.set_rotate_label(False)
            ax.set_zlabel(r'$\psi$'+str(indices)+"(x,y)",rotation=90)
            plt.show()
            pass

        return None
