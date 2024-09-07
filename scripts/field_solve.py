# Date: 07 June 2024
# Author: Anupsa Swain

import os
import numpy as np
import matplotlib.pyplot as plt
from single_beam_lattice_create import meshplot
from dynamic_truss_lib import PlaneTruss
from solve_microtruss import assign_constraint, assign_load, tlib
import time


if __name__ == "__main__":
    # define a field, solve a bar truss with it

    t0 = time.time() # start
    L = 10
    rho = 2710
    E = 7e10
    A = 0.25

    for layers in [2, 4]:
        meshsize = 0.1
        conn_node_tol = 0.4
        H_meshsize = meshsize*(3**.5)/2
        w = H_meshsize*layers

        u_time_array = np.load('data/dynamic/l%d_m%.3f_displacement_data.npz' %(layers, meshsize))['x']

        t = 0.02/5
        # in 0.1 seconds, 201 elements. in t seconds, t/0.1*201
        t_idx = int(t/0.1*201)

        q0 = u_time_array[t_idx]

        nodes = np.load('data/dynamic/l%d_m%.3f_final_trussCord.npz' %(layers, meshsize))
        nodes = nodes['x']

        edges = np.load('data/dynamic/l%d_m%.3f_final_node_connect.npz' %(layers, meshsize))
        edges = edges['x']

        # nodes and edge check
        # meshplot(nodes, edges)

        edges = edges - 1

        macro_nodes = np.array([[0, 0], [L, 0]])

        macro_truss = PlaneTruss(macro_nodes, np.array([[0, 1]]), E, A, rho)

        A_per = A #/layers

        E_truss = E /(2.5*layers + 1)
        rho_truss = rho /(2.5*layers + 1)

        truss = PlaneTruss(nodes, edges, E_truss, A_per, rho_truss)
        static_truss = tlib.PlaneTruss(nodes, edges, E_truss, A)

        macro_constraints = [
            [0, 0, 0.0],
            [0, 1, 0.0],
        ]
        constraints = assign_constraint(macro_constraints, macro_nodes, nodes, w)

        truss.apply_constraints(constraints)

        t1 = time.time() # setup complete

        qs, qdots, qddots = truss.solve(h=(0.1/200), n_iterations=200, q0=q0)


        file_name = os.path.join('data/dynamic/', 'l%d_m%.3f_field_displacement_data.npz' %(layers, meshsize))
        np.savez(file_name, x=qs, y=qdots, z=qddots)

    t2 = time.time() # solve complete

    print("Time taken: " + str(t2 - t0) + "s")