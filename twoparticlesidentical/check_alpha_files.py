import numpy as np
import os

N = 150
alphas = np.linspace(0, 1, 25)
base = r'C:\Users\vq22287\OneDrive - University of Bristol\Documents\Anyons on graphs\Eigenstates data'

for alpha in alphas:
    states_stub = os.path.join(base, 'cross_hardcore_states', 'cross_N' + str(N) + '_alpha' + str(alpha) + '_')
    eigs_path = os.path.join(base, 'cross_hardcore_eigenvalues', 'cross_N' + str(N) + '_alpha' + str(alpha) + '.csv')
    states_exists = os.path.exists(states_stub + 'D(1, 1).csv')
    eigs_exists = os.path.exists(eigs_path)
    if not (states_exists and eigs_exists):
        print('missing', alpha, states_exists, eigs_exists)
        break
else:
    print('all present')
