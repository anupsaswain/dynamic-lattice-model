# Author: Anupsa Swain
# Date: 28 May 2024

import os
import numpy as np
import matplotlib.pyplot as plt
from single_beam_lattice_create import meshplot
from dynamic_truss_lib import PlaneTruss
from solve_microtruss import assign_constraint, assign_load, tlib
import time

# def p(t):
#     p = np.zeros(len(nodes)*2)
#     p[-1] = 1e4
#     return p

if __name__ == '__main__':
    t0 = time.time() # start
    L = 10
    rho = 2710
    E = 7e10
    A = 0.25

    for layers in [1]:
        meshsize = .1
        conn_node_tol = 0.5
        H_meshsize = meshsize
        w = H_meshsize*layers

        nodes = np.load('data/dynamic/new_square/l%d_m%.3f_final_trussCord.npz' %(layers, meshsize))
        nodes = nodes['x']

        edges = np.load('data/dynamic/new_square/l%d_m%.3f_final_node_connect.npz' %(layers, meshsize))
        edges = edges['x']

        # nodes and edge check
        # meshplot(nodes, edges)

        edges = edges - 1

        macro_nodes = np.array([[0, 0], [L, 0]])

        macro_truss = PlaneTruss(macro_nodes, np.array([[0, 1]]), E, A, rho)

        A_per = A #/layers

        # horizontal
        E_truss = E /(layers*2 + 1)
        rho_truss = rho /(layers*2 + 1)

        # # vertical
        # E_truss = E*10 /(layers**2.5 + 1)
        # rho_truss = rho*10 /(layers + 1)

        truss = PlaneTruss(nodes, edges, E_truss, A_per, rho_truss)
        static_truss = tlib.PlaneTruss(nodes, edges, E_truss, A)

        macro_constraints = [
            [0, 0, 0.0],
            [0, 1, 0.0],
            [-1, 0, 0.0],
            [-1, 1, 0.0]
        ]

        P = 1e9
        Q = 0
        macro_loads = [
            [0, 0, -P]
        ]

        loads = assign_load(macro_loads, macro_nodes, nodes, w)
        constraints = assign_constraint(macro_constraints, macro_nodes, nodes, w)

        truss.apply_constraints(constraints)
        truss.apply_loads(loads)

        truss.compute_stiffness()
        truss.compute_mass()
        truss.enforce_constraints()

        K = truss.K

        a = np.linalg.eigvals(K)

        # static_truss.apply_constraints(constraints)
        # static_truss.apply_loads(loads)

        # macro_truss.apply_constraints(macro_constraints)
        # macro_truss.apply_loads(macro_loads)

        t1 = time.time() # setup complete

        # qs, qdots, qddots = truss.solve(h=(0.1/100), n_iterations=100)
        # macro_truss.solve(h=0.05, n_iterations=200)
        # static_truss.solve()

        # file_name = os.path.join('data/dynamic/new_square/', 'l%d_m%.3f_antiresonance_data.npz' %(layers, meshsize))
        # np.savez(file_name, x=qs, y=qdots, z=qddots)

    t2 = time.time() # solve complete

    # q0dot = np.zeros([2, truss.n_joints])
    # qs, qdots, qddots = bar_newmark_integration(joints, truss.M, truss.C, truss.K, q0dot, p, h=0.05, n_iterations=200) # need to redefine same thing within truss class definition.

    # print('Reactions: ')
    # print(truss.reactions)

    # print('Internal forces: ')
    # print(truss.member_forces)

    print("Time taken: " + str(t2 - t0) + "s")
    print(np.sort(np.sqrt(a)/2/np.pi)[:20])

    # print("stiffness matrix, K: ")
    # print(truss.K)
    # print("mass matrix, M: ")
    # print(truss.M)
    # print("damping matrix, C: ")
    # print(truss.C)
    # print("Force vector, F:")
    # print(truss.F)


    # truss.plot(mag=25)
    # macro_truss.plot()
    # static_truss.plot()

    # t3 = time.time() # plot complete; also includes time plot was open and stared at







