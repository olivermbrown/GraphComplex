import numpy as np
import sympy as sym

def set_P():
    
    P = np.zeros((6,6))
    
    P[0,0] = 2/3
    P[1,0] = -1/3
    P[2,0] = -1/3
    P[0,1] = -1/3
    P[1,1] = 2/3
    P[1,2] = -1/3
    P[0,2] = -1/3
    P[2,1] = -1/3
    P[2,2] = 2/3
    
    P[3,3] = 1
    P[4,4] = 1
    P[5,5] = 1
    
    return P

def set_T(k):
    
    ek = sym.exp(sym.I*k)
    
    t = sym.diag(ek,ek,ek)
    
    O = sym.zeros(3,3)
    
    T = sym.Matrix([[O, t],[t, O]])
    
    return T

def Scattering_Matrix(P,L,k):
    
    Q = np.eye(6) - P
    
    S = - P - Q * sym.Matrix(L - sym.I*k*sym.eye(6)).inv() * (L + sym.I*k*sym.eye(6)) * Q
    
    return S

if __name__ == "__main__":
    
    k = sym.symbols('k')
    
    P = set_P()
    
    L = np.zeros((6,6))
    
    S = Scattering_Matrix(P,L,k)
    
    T = set_T(k)
    
    U = S * T
    
    V = sym.eye(6) + U
    
    F = V.det()
    
    sym.det
    
    print("P = ")
    print(P)
    
    print("L = ")
    print(L)
    
    print("S = ")
    print(S)
    
    print("U = ")
    print(U)
    
    print("V = ")
    print(V)
    
    print("Determinant F(k) = ")
    print(F)
    
    pass