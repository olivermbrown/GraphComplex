"""
Module to calculate the number of n-particle domains from arrangements of the n particles on a graph with E wires
"""

import numpy as np
import scipy as sp
import more_itertools as mit
from functools import reduce

def number_of_domains(N, E, arr):
    # Calculate the number of domains from an arrangement arr of n particles on the graph

    try:
        m = len(arr)

        if len(arr)==E:
            pass
        elif len(arr)<=E:
            diff = E - len(arr)
            zeros = tuple(np.zeros(diff,dtype=int))
            arr = arr + zeros
            pass
        else:
            raise Exception
        pass

    except:
        print("Error: The arrangement is not valid")
        return None
    
    else:
        # Split the arrangement into blocks where each site contains the same number of particles

        partition = mit.split_when(arr, lambda x,y: y-x != 0)
        b = tuple(map(len, partition))
        
        # Calculate the number of domains that satisfy the arrangement arr

        A = list(map(np.math.factorial, arr))
        B = list(map(np.math.factorial, b))
        nd = (np.math.factorial(N)*np.math.factorial(E))/(np.prod(np.array(A))*np.prod(np.array(B)))
        return nd
        

if __name__=="__main__":
    N = 3
    E = 3
    arrangement = (3,)
    print(number_of_domains(N,E,arrangement))
    pass