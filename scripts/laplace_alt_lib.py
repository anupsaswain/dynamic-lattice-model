import time
import numpy as np
import dynamic_truss_lib as dlib
import sympy as sp
from mpmath import invertlaplace
import matplotlib.pyplot as plt

sp.init_printing()
s = sp.Symbol('s')
t = sp.Symbol('t', positive=True)

def numerical_inverse_laplace(v_exprs, s, t_vals):
    # Convert symbolic expressions to numerical functions
    v_funcs = [sp.lambdify(s, expr, 'mpmath') for expr in v_exprs]
    
    # Define the numerical inverse Laplace transform using Talbot method
    def inverse_laplace_transform_num(f_s, t_val):
        return float(invertlaplace(f_s, t_val))
    
    # Compute numerical inverse Laplace transform for each function
    v_new = []
    for v_func in v_funcs:
        v_num = []
        for t_val in t_vals:
            v_num.append(inverse_laplace_transform_num(v_func, t_val))
        v_new.append(v_num)
    
    return v_new

def laplace_u(u_array, t_array):
    n, m = u_array.shape
    laplace_result = []
    dt = t_array[1] - t_array[0]
    for i in range(n):
        u_laplace = 0
        for j in range(m):
            u_laplace += u_array[i, j] * sp.exp(-s * t_array[j]) * dt
        laplace_result.append(u_laplace)
    return sp.Matrix(laplace_result)

def tf_unit(u,t_vals, K, M, C=0):
    mid_n = int(len(M) / 2)
    M21 = M[mid_n:, :mid_n]
    K21 = K[mid_n:, :mid_n]

    M22 = M[mid_n:, mid_n:]
    K22 = K[mid_n:, mid_n:]

    A21 = (s * s) * M21 + K21
    A22 = (s * s) * M22 + K22

    A21 = sp.Matrix(A21)
    A22 = sp.Matrix(A22)

    # t_vals = np.linspace(0.01, 0.1, 100)  # Define time range

    u = laplace_u(u, t_vals)
    print("laplace done")

    q = A21 @ u

    v = -A22.inv() @ q

    # Use numerical inverse Laplace transform
    v_new = numerical_inverse_laplace(v, s, t_vals)

    return np.array(v_new)

if __name__ == "__main__":
    t0 = time.time()  # start
    L = 10.0
    rho = 2710.0
    E = 7.0e10
    A = 0.25

    f = 1

    t_vals = np.linspace(0.1, 1, 100)
    u0 = np.array([
        np.sin(2*np.pi*f*t_vals),
        np.zeros_like(t_vals),
        np.sin(2*np.pi*f*t_vals),
        np.zeros_like(t_vals)
    ])

    macro_nodes = np.array([[0, L], [0, 0], [L, L], [L, 0]])

    macro_edges = np.array([[0, 1], [1, 2], [2, 3], [3, 0], [1, 3], [0, 2]])

    truss = dlib.PlaneTruss(macro_nodes, macro_edges, E, A, rho)
    truss.compute_mass()
    truss.compute_stiffness()

    v = tf_unit(u0 / 100, truss.K, truss.M)
    t_vals = np.linspace(0.01, 1, 100)

    for i in range(4):
        plt.plot(t_vals, v[i])

    plt.legend(["p1_x", "p1_y", "p2_x", "p2_y"])

    # print(v)

    t1 = time.time()

    print("Time taken: " + str(t1 - t0) + " seconds")
    plt.show()
