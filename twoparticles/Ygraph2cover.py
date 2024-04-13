"""
Code to solve the Laplacian on the distinguishable configuration space of two particles on the Y graph
"""

import numpy as np
import scipy as sp
import networkx as nx
from functools import reduce

import matplotlib.pyplot as plt
from matplotlib import cm

import squareFDM
import twocomplex
import configspace
"""
def twoparticlesYgraph2cover(N):
    # Evaluate the spectrum of the Laplacian on the double cover of the distinguishable configuration space of two particles on the Y graph

    # Scaling factor
    h = (np.pi)/(N-1)

    # Define the indices of the graph
    wires = np.array([0,1,2])

    # Create the configurations
    configurations = np.zeros((len(wires),len(wires),2))
    for i in wires:
        for j in wires:
            configurations[i,j] = np.array([i,j])
            pass
        pass

    domains = []

    # Create the domains
    configurations = np.reshape(configurations,(len(wires)**2,2))
    for config in configurations:
        D = squareFDM.Domain(N)
        if config[0] == config[1]:
            D.split_domain()
            pass
        else:
            pass
        D.indices = tuple(config)
        domains.append(D)
        pass

    # Reshape the domains
    #configurations = np.reshape(configurations,(len(wires),len(wires)),order='F')
    domains = np.reshape(domains,(len(wires),len(wires)))

    # Create the double cover
    DC = np.array([domains,domains])

    """"""
    Old code
    # Create domains
    D11 = squareFDM.Domain(N)
    D22 = squareFDM.Domain(N)
    D33 = squareFDM.Domain(N)
    D12 = squareFDM.Domain(N)
    D13 = squareFDM.Domain(N)
    D23 = squareFDM.Domain(N)
    D21 = squareFDM.Domain(N)
    D31 = squareFDM.Domain(N)
    D32 = squareFDM.Domain(N)
    D11.split_domain()
    D22.split_domain()
    D33.split_domain()
    D11.indices = (1,1)
    D22.indices = (2,2)
    D33.indices = (3,3)
    D12.indices = (1,2)
    D13.indices = (1,3)
    D23.indices = (2,3)
    D21.indices = (2,1)
    D31.indices = (3,1)
    D32.indices = (3,2)
    
    
    cells = [D11,D22,D33,D12,D13,D23,D21,D31,D32]
    """"""

    # Define gluings for the domains

    gluings = []

    gluings.append([DC[0,0,0].x0, DC[1,0,2].x0, DC[0,0,1].x0])
    gluings.append([DC[0,1,1].y0, DC[0,0,1].y0, DC[0,2,1].y0])
    gluings.append([DC[0,2,2].x0, DC[0,2,1].x0, DC[0,2,0].x0])
    gluings.append([DC[0,0,0].y0, DC[0,2,0].y0, DC[0,1,0].y0])
    gluings.append([DC[0,1,1].x0, DC[0,1,0].x0, DC[0,1,2].x0])
    gluings.append([DC[0,2,2].y0, DC[0,1,2].y0, DC[0,0,2].y0])

    gluings.append([DC[1,0,0].x0, DC[0,0,2].x0, DC[1,0,1].x0])
    gluings.append([DC[1,1,1].y0, DC[1,0,1].y0, DC[1,2,1].y0])
    gluings.append([DC[1,2,2].x0, DC[1,2,1].x0, DC[1,2,0].x0])
    gluings.append([DC[1,0,0].y0, DC[1,2,0].y0, DC[1,1,0].y0])
    gluings.append([DC[1,1,1].x0, DC[1,1,0].x0, DC[1,1,2].x0])
    gluings.append([DC[1,2,2].y0, DC[1,1,2].y0, DC[1,0,2].y0])

    """"""
    Old code
    gluingA = [D11.x0, D13.x0, D12.x0]
    gluingB = [D22.y0, D12.y0, D32.y0]
    gluingC = [D33.x0, D32.x0, D31.x0]
    gluingD = [D11.y0, D31.y0, D21.y0]
    gluingE = [D22.x0, D21.x0, D23.x0]
    gluingF = [D33.y0, D23.y0, D13.y0]
    
    Network.glue(gluingA)
    Network.glue(gluingB)
    Network.glue(gluingC)
    Network.glue(gluingD)
    Network.glue(gluingE)
    Network.glue(gluingF)
    """"""

    cells = list(DC.flatten())

    # Create two complex from domains
    Network = twocomplex.Complex(cells)

    Network.gluings = gluings

    # Apply boundary conditions to the exterior and diagonal of the complex
    Network.diagonal_bc("dirichlet")
    Network.exterior_bc("dirichlet")
    
    # Generate the discrete Laplacian
    Network.gen_lapl()
    
    return Network
"""

def bosonic_projection(Complex, n):
        # Project the nth eigenstate onto the bosonic subspace

        if Complex.solved == True:
            pass
        else:
            Complex.lapl_solve(2,dps=2)
            pass

        cells = Complex.cells
        states = Complex.states

        # Find coordinates in the complex of each cell
        state_components = []
        c = 0
        for cell in cells:
            
            Complex.construct_coords(cell)

            length = len(cell.non_elim_coords)

            state_components.append(states[c:length+c,n])

            c += length
            pass

        #state_components = np.vstack(state_components)

        # Construct the bosonic state
        bosonic_state = np.zeros(len(states[n]))

        for cell in cells:
            i = cells.index(cell)
            if cell.split == True:
                for cell2 in cells:
                    j = cells.index(cell2)
                    if cell.indices[0:2] == cell2.indices[0:2] and cell != cell2 and cell2.split == True:
                        state_components[i] += state_components[j]
                        pass
                    elif cell.indices[0:2] == tuple(reversed(cell2.indices[0:2])) and cell != cell2 and cell2.split == True:
                        state_components[i] += state_components[j]
                        pass
                    else:
                        pass
                    pass
                pass
            else:
                for cell2 in cells:
                    j = cells.index(cell2)
                    if cell.indices[0:2] == cell2.indices[0:2] and cell != cell2 and cell2.split == False:
                        state_components[i] += state_components[j]
                        pass
                    elif cell.indices[0:2] == tuple(reversed(cell2.indices[0:2])) and cell != cell2 and cell2.split == False:
                        state_components[i] += state_components[j]
                        pass
                    else:
                        pass
                    pass
                pass
            pass

        # Find cells corresponding to the same domains in the indistinguishable configuration space
        # and sum the states on those cells

        bosonic_states = []

        c = 0
        for cell in cells:
            i = cells.index(cell)

            Complex.construct_coords(cell)

            coords = cell.non_elim_coords
            length = cell.num_non_elim

            #state = state_components[c:length+c]
            state = state_components[i]
            state = np.real(state)

            ax = plt.figure().add_subplot(projection='3d')

            ax.plot_trisurf(coords[:,0],coords[:,1], state, linewidth=0.2, antialiased=True, cmap=cm.inferno)
            #ax.scatter(coords[:,0],coords[:,1],state,c=state,cmap=cm.autumn)
            indices = cell.indices

            plt.title("D"+str(indices))
            ax.set_xlabel("x_e"+str(indices[0]))
            ax.set_ylabel("y_e"+str(indices[1]))
            ax.zaxis.set_rotate_label(False)
            ax.set_zlabel(r'$\psi$'+str(indices)+"(x,y)",rotation=90)
            plt.show()
            c+=length
            pass

        return None

def twoparticlesYgraph2cover(N):
    # Evaluate the spectrum of the Laplacian on the distinguishable configuration space of two particles on the Y graph

    # Scaling factor
    h = (np.pi)/(N-1)

    # Create domains
    D11_0 = squareFDM.Domain(N)
    D22_0 = squareFDM.Domain(N)
    D33_0 = squareFDM.Domain(N)
    D12_0 = squareFDM.Domain(N)
    D13_0 = squareFDM.Domain(N)
    D23_0 = squareFDM.Domain(N)
    D21_0 = squareFDM.Domain(N)
    D31_0 = squareFDM.Domain(N)
    D32_0 = squareFDM.Domain(N)
    D11_0.split_domain()
    D22_0.split_domain()
    D33_0.split_domain()
    D11_0.indices = (1,1,0)
    D22_0.indices = (2,2,0)
    D33_0.indices = (3,3,0)
    D12_0.indices = (1,2,0)
    D13_0.indices = (1,3,0)
    D23_0.indices = (2,3,0)
    D21_0.indices = (2,1,0)
    D31_0.indices = (3,1,0)
    D32_0.indices = (3,2,0)
    D11_1 = squareFDM.Domain(N)
    D22_1 = squareFDM.Domain(N)
    D33_1 = squareFDM.Domain(N)
    D12_1 = squareFDM.Domain(N)
    D13_1 = squareFDM.Domain(N)
    D23_1 = squareFDM.Domain(N)
    D21_1 = squareFDM.Domain(N)
    D31_1 = squareFDM.Domain(N)
    D32_1 = squareFDM.Domain(N)
    D11_1.split_domain()
    D22_1.split_domain()
    D33_1.split_domain()
    D11_1.indices = (1,1,1)
    D22_1.indices = (2,2,1)
    D33_1.indices = (3,3,1)
    D12_1.indices = (1,2,1)
    D13_1.indices = (1,3,1)
    D23_1.indices = (2,3,1)
    D21_1.indices = (2,1,1)
    D31_1.indices = (3,1,1)
    D32_1.indices = (3,2,1)
    
    cells = [D11_0,D22_0,D33_0,D12_0,D13_0,D23_0,D21_0,D31_0,D32_0,D11_1,D22_1,D33_1,D12_1,D13_1,D23_1,D21_1,D31_1,D32_1]
    
    # Create two complex from domains
    Network = twocomplex.Complex(cells)
    
    # Glue the domains together
    gluingA = [D11_0.x0, D13_1.x0, D12_0.x0]
    gluingB = [D22_0.y0, D12_0.y0, D32_0.y0]
    gluingC = [D33_0.x0, D32_0.x0, D31_0.x0]
    gluingD = [D11_0.y0, D31_0.y0, D21_0.y0]
    gluingE = [D22_0.x0, D21_0.x0, D23_0.x0]
    gluingF = [D33_0.y0, D23_0.y0, D13_0.y0]
    gluingG = [D11_1.x0, D13_0.x0, D12_1.x0]
    gluingH = [D22_1.y0, D12_1.y0, D32_1.y0]
    gluingI = [D33_1.x0, D32_1.x0, D31_1.x0]
    gluingJ = [D11_1.y0, D31_1.y0, D21_1.y0]
    gluingK = [D22_1.x0, D21_1.x0, D23_1.x0]
    gluingL = [D33_1.y0, D23_1.y0, D13_1.y0]
    
    Network.glue(gluingA)
    Network.glue(gluingB)
    Network.glue(gluingC)
    Network.glue(gluingD)
    Network.glue(gluingE)
    Network.glue(gluingF)
    Network.glue(gluingG)
    Network.glue(gluingH)
    Network.glue(gluingI)
    Network.glue(gluingJ)
    Network.glue(gluingK)
    Network.glue(gluingL)

    # Apply boundary conditions to the exterior and diagonal of the complex
    Network.diagonal_bc("dirichlet")
    Network.exterior_bc("dirichlet")
    
    # Generate the discrete Laplacian
    Network.gen_lapl()
    
    return Network

if __name__ == "__main__":
    # Main
    
    N = 20 # The length of a wire

    # Scaling factor
    h = (np.pi)/(N-1)

    Ygraph = twoparticlesYgraph2cover(N)

    #Ygraph.print_eqs()

    #spectrum = Ygraph.lapl_spectrum(h,2,200)
    #print(spectrum)

    Ygraph.lapl_solve(h,2,20)
    print("Unsorted eigenvalues = " + str(Ygraph.spectrum))
    print("Sorted eigenvalues = " + str(np.sort(Ygraph.spectrum)))

    # Plot the states
    #Ygraph.plot_states(0)

    # Project the states onto the bosonic subspace
    #bosonic_projection(Ygraph,1)
    Configs = configspace.ConfigurationSpace(Ygraph)

    Configs.bosonic_projection()
    Configs.plot_projected_wavefunction(0)
    
    pass
