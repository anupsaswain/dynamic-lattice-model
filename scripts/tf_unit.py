# 17 June 2024
# Anupsa Swain

import os
import time
import numpy as np
import matplotlib.pyplot as plt
import dynamic_truss_lib as dlib
import scipy.signal as ss
from dynamic_lattice_model.laplace_alt_lib import tf_unit

def meshplot(trussCord, node_connect):
    fig, ax = plt.subplots()
    i = 0
    for element in node_connect:
        x1, x2 = trussCord[element[0], 0], trussCord[element[1], 0]
        y1, y2 = trussCord[element[0], 1], trussCord[element[1], 1]
        ax.plot([x1,x2], [y1,y2], '-k')
        # plt.annotate(i, xy=((x1+x2)/2, (y1+y2)/2), c='red')
        i = i + 1
    x_scatter = trussCord[:,0]
    y_scatter = trussCord[:,1]
    x_scatter = np.vstack([x_scatter]).flatten()
    y_scatter = np.vstack([y_scatter]).flatten()
    plt.scatter(x_scatter, y_scatter)
    ax.set_aspect('equal')
    # ax.set_xlim([-2, 12])
    # ax.set_ylim([-2, 2])

    for i in range(len(trussCord)):
        plt.annotate(i, xy=(trussCord[i, 0], trussCord[i, 1]), c='blue')

    # for i in range(len(node_connect)):
    #     plt.annotate(i, xy=(node_connect[i, 0], node_connect[i, 1]), c='red')

    plt.show()

def u_t(f, t):
    # t is a timestamp
    # f = 50
    A = 0.5/((f/50)**2)
    return 0.5*np.sin(f*2*np.pi*t)

def tf_static_unit(u, omega, K, M, C=0):
    n = len(M)
    mid_n = int(n/2)
    M21 = M[mid_n:, :mid_n]
    K21 = K[mid_n:,:mid_n]

    M22 = M[mid_n:,mid_n:]
    K22 = K[mid_n:,mid_n:]

    A21 = -omega**2*M21 + K21
    A22 = -omega**2*M22 + K22

    # A21  = ss.TransferFunction(np.array([M21, 0, K21]))
    # A22  = ss.TransferFunction(np.array([M22, 0, K22]))

    # print(A21)

    T = -np.linalg.inv(A22)@A21

    v = T@u

    return v

if __name__ == "__main__":
    t0 = time.time() # start
    L = 10
    rho = 2710
    E = 7e10
    A = 0.25

    layers = 1
    meshsize = 1

    n_cells = L/meshsize

    nodes = np.load('data/dynamic/new_square/l%d_m%.3f_final_trussCord.npz' %(layers, meshsize))['x']
    edges = np.load('data/dynamic/new_square/l%d_m%.3f_final_node_connect.npz' %(layers, meshsize))['x']
    
    edges = edges - 1

    E_truss = E /(layers*2 + 1)
    rho_truss = rho /(layers*2 + 1)

    # meshplot(nodes, edges)

    t_array = np.linspace(0.01, 0.1, 100)

    len_t = len(t_array)

    for freq in [25, 75, 125, 175]:
        u0 = np.zeros([int(2*(layers+1)), len_t])
        for i in range(layers + 1):
            for j in range(len_t):
                u0[2*i, j] = u_t(freq, t_array[j])

        u_array = np.zeros([len(nodes)*2, len_t])

        for i in range(int(len(nodes)/(layers+1) - 1)):
            if i == 0: ui = u0
            macro_nodes = nodes[i:i+ int(2*(layers + 1))]
            macro_edges = np.zeros([0, 2])
            print(i)
            for j in range(layers):
                print(j)
                idx = int((layers+1)*i + j)
                u_array[idx:(idx+2)] = ui[int(2*j): int(2*j + 2)]
                macro_edges = np.append(macro_edges, np.array([[j, j+1], 
                                                                [j, j + layers + 1], 
                                                                [j, j + layers + 2], 
                                                                [j + layers + 1, j + layers+2], 
                                                                [j+1, j+layers+1]]), axis=0)
            macro_edges = np.append(macro_edges,np.array([[layers, 2*layers+1]]), axis=0)
            idx = int((layers+1)*i + layers)
            u_array[idx:(idx+2)] = ui[int(2*layers): int(2*layers + 2)]
            macro_edges = macro_edges.astype('int32')
            # print(macro_edges)
            truss = dlib.PlaneTruss(macro_nodes, macro_edges, E_truss, A, rho_truss)
            truss.compute_mass()
            truss.compute_stiffness()

            # ui = tf_static_unit(ui, 2*np.pi*10, truss.K, truss.M)
            ui = tf_unit(ui, t_array, truss.K, truss.M)

        file_name = os.path.join('data/dynamic/new_square/amplitude/', 'l%d_m%.3f_tf_data_%gHz.npz' %(layers, meshsize, freq))
        np.savez(file_name, x=u_array)

        print(np.shape(ui))
        t1 = time.time()

        print("Time taken: " + str(t1 - t0) + " seconds")

    # moved_nodes = nodes + u_array
    # meshplot(moved_nodes, edges)


