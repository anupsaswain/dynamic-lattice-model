# 19 June 2024
# Anupsa Swain

import os
import time
import numpy as np
import matplotlib.pyplot as plt
import dynamic_truss_lib as dlib
import scipy.signal as ss
from tf_unit import tf_static_unit, meshplot
import sympy as sp
# from sympy.abc import s, t
from mpmath import invertlaplace

sp.init_printing()
s = sp.Symbol('s')
t = sp.Symbol('t', positive=True)
# a = sp.Symbol('a')

# t, s = sp.symbols('t, s')

def tf_unit(u, K, M, C=0):

    # t, s = sp.Symbol('t, s')
    mid_n = int(len(M)/2)
    M21 = M[mid_n:, :mid_n]
    K21 = K[mid_n:, :mid_n]

    M22 = M[mid_n:,mid_n:]
    K22 = K[mid_n:,mid_n:]

    A21 = (s*s)*M21 + K21
    A22 = (s*s)*M22 + K22

    A21 = sp.Matrix(A21)
    A22 = sp.Matrix(A22)

    q = A21@u

    v = -A22.inv()@q


    # v = list(v)

    v_new = []

    for i in range(len(v)):
        print(i)
        v[i] = sp.simplify(v[i])
        q = v[i]
        b = sp.inverse_laplace_transform(q, s, t) # it's failing here!!!
        v_new.append(b)

    # def v_func(s_var, i):
    #     q = v[i]
    #     a = q.subs(s, s_var)
    #     return a
    
    # for i in range(len(v)):
    #     def v_fun_element(s_var):
    #         return v_func(s_var, i)
    #     v_new.append(v_fun_element) 

    return v_new

# def tf_unit_for_specific_truss(s):
#     return tf_unit(s, u0, truss.K, truss.M)

def inverse_v(v, time):
    v_new = []

    for i in range(len(v)):
        print(i)
        q = v[i]
        fs = lambda s: q(s)
        a = invertlaplace(q, time)
        v_new.append(a)
    return v_new

if __name__ == "__main__":
    t0 = time.time() # start
    L = 10.0
    rho = 2710.0
    E = 7.0e10
    A = 0.25

    u0 = np.array([0.5, 0, 0.5, 0])

    macro_nodes = np.array([[0, L], [0, 0], [L, L], [L, 0]])

    macro_edges = np.array([[0, 1], [1, 2], [2, 3], [3, 0], [1, 3], [0, 2]])

    truss = dlib.PlaneTruss(macro_nodes, macro_edges, E, A, rho)
    truss.compute_mass()
    truss.compute_stiffness()

    v = tf_unit(u0, truss.K, truss.M)

    # v is an array of functions with variable s
    # v = inverse_v(v, 0.1)

    # u_array = np.append(u0, v, axis=0)

    # u_array = np.reshape(u_array, [4, 2])

    # moved_nodes = macro_nodes + u_array

    # meshplot(moved_nodes, macro_edges)

    print(v)
    
    t1 = time.time()

    print("Time taken: " + str(t1-t0) + " seconds")
