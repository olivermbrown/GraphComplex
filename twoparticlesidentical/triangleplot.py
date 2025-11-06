import sympy
from sympy.plotting import plot3d
import numpy as np

def eigenstate(x, y, a, k1, k2):
    """
    Calculate the eigenstate of a system given two wave vectors k1 and k2.
    
    Parameters:
    k1 (float): The first wave vector.
    k2 (float): The second wave vector.
    
    Returns:
    sympy.Matrix: The eigenstate represented as a matrix.
    """
    # Define the variables x, y
    x, y = sympy.symbols('x y')

    # Define the wave function as a complex function of x and y
    wave_function = sympy.exp(1j * (k1 * x + k2 * y)) + sympy.exp(1j * (k2 - sympy.pi * a + k2 * x + k1 * y))
    
    return wave_function

if __name__ == "__main__":

    # Choose a value for a
    a = 0.0

    # Choose k1 and k2 values
    k1 = a*sympy.pi
    k2 = a*sympy.pi + 2*sympy.pi

    # Define the variables x, y and alpha
    x, y = sympy.symbols('x y')

    # Calculate the eigenstate
    state = sympy.arg(eigenstate(x, y, a, k1, k2))

    # Plot the eigenstate in 3d
    plot3d(state, (x, 0, 1), (y, 0, 1), title='Eigenstate Plot', xlabel='x', ylabel='y')