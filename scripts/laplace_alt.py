# 02 July 2024
# Anupsa Swain

import os
import time
import numpy as np
import matplotlib.pyplot as plt
import dynamic_truss_lib as dlib
import sympy as sp
from laplace_alt_lib import laplace_u, numerical_inverse_laplace
from scipy.fft import fft, ifft, fftfreq

sp.init_printing()
s = sp.Symbol('s')
t = sp.Symbol('t', positive=True)

def u_t(f, t):
    # t is a timestamp
    # f = 50
    return 0.5*np.sin(f*2*np.pi*t)

def laplace_output(u_input, truss, t_array):
    lt0 = time.time()
    K = truss.K
    M = truss.M 

    ni = len(u_input)

    u_s = laplace_u(u_input, t_array)

    lt1 = time.time() - lt0
    print("Laplace: " + str(lt1) + "s")

    A = (s*s)*M + K
    A = sp.Matrix(A)

    A11 = A[:ni, :ni]
    A12 = A[:ni, ni:-ni]
    A13 = A[:ni, -ni:]

    A21 = A[ni:-ni, :ni]
    A22 = A[ni:-ni, ni:-ni]
    A23 = A[ni:-ni, -ni:]

    A31 = A[-ni:, :ni]
    A32 = A[-ni:, ni:-ni]
    A33 = A[-ni:, -ni:]

    lt2 = time.time() - lt0
    print("Matrix A ready: " + str(lt2) + "s")

    A22 = sp.Matrix(A22)
    print("Length of A22: " + str(len(A22)))

    A22_inv = A22.inv()

    lt3 = time.time() - lt0
    print("A22 inverted: " + str(lt3) + "s")

    B_in = A31 - A32@(A22_inv@A21)
    B_out = A33 - A32@(A22_inv@A23)

    lt4 = time.time() - lt0
    print("B matrices ready: " + str(lt4) + "s")

    v_s = -B_out.inv()@(B_in@u_s)

    lt5 = time.time() - lt0
    print("v_s calculated: " + str(lt5) + "s")

    # v_s = some linalg with * u_s

    v_t = numerical_inverse_laplace(v_s, s, t_array)

    lt6 = time.time() - lt0
    print("Laplace inversion done: " + str(lt6) + "s")

    return np.array(v_t)


if __name__ == "__main__": 
    t0 = time.time()
    L = 10
    rho = 2710
    E = 7e10
    A = 0.25

    layers = 1
    meshsize = 1
    freq = 50

    nodes = np.load('data/dynamic/new_square/l%d_m%.3f_final_trussCord.npz' %(layers, meshsize))['x']
    edges = np.load('data/dynamic/new_square/l%d_m%.3f_final_node_connect.npz' %(layers, meshsize))['x']
    
    edges = edges - 1

    E_truss = E /(layers*2 + 1)
    rho_truss = rho /(layers*2 + 1)

    t_array = np.linspace(0.01, 0.1, 100)
    len_t = len(t_array)

    truss = dlib.PlaneTruss(nodes, edges, E_truss, A, rho_truss)
    truss.compute_mass()
    truss.compute_stiffness()

    u_in = np.zeros([int(2*(layers+1)), len_t])
    for i in range(layers + 1):
        for j in range(len_t):
            u_in[2*i, j] = u_t(freq, t_array[j])

    u_out = laplace_output(u_in, truss, t_array)

    file_name = os.path.join('data/dynamic/new_square/laplace/', 'l%d_m%.3f_laplace_data_f%gHz.npz' %(layers, meshsize, freq))
    np.savez(file_name, x=u_out, t=t_array)

    u2L_array = np.average(u_out, axis=0)

    x = t_array
    y = u2L_array

    N = len(t_array)
    T = t_array[1] - t_array[0]

    yf = fft(y)
    xf = fftfreq(N, T)[:N//2]

    plt.plot(xf[:], 2.0/N* np.abs(yf[:N//2]))
    plt.grid()
    plt.show()