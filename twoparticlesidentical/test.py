import numpy as np
import scipy as sp
import configs as cfg
import dumbbell as db

if __name__ == "__main__":
    print("Testing scipy linalg")

    # Generate 1000 x 1000 sparse matrix of zeros
    #A = sp.sparse.csr_matrix((1000,1000))

    # Put some values in the matrix
    #for i in range(100):
    #    A[i,i] = i*10
    #    pass

    # Generate a dumbbell system
    N = 30
    h = (np.pi)/(N-1)
    alpha1 = 0.5
    alpha2 = 0.25

    C = db.DumbbellAnyonsHardcore(N,alpha1,alpha2)

    C.gen_lapl()

    C.simplify_lapl()

    A = C.sL

    A = A.astype(np.complex64)
    A.tocsr()

    # Solve the system
    eigvals, eigvecs = sp.sparse.linalg.eigs(A, k=10, which='SM', return_eigenvectors=True)

    print(eigvals)