import numpy as np
import matplotlib.pyplot as plt
from single_beam_lattice_create import meshplot
from dynamic_truss_lib import PlaneTruss
from solve_microtruss import assign_constraint, assign_load, tlib
import time

if __name__ == '__main__':
    t0 = time.time() # start
    L = 10
    rho = 2710
    E = 7e10
    A = 0.25

    macro_nodes = np.array([[0, 0], [L, 0], [L/2, L/2], [L/2, -L/2]])

    macro_edges = np.array([[0, 1], [0, 2], [3, 0], [2, 1], [3, 1]])

    macro_truss = PlaneTruss(macro_nodes, macro_edges, E, A, rho)
    truss = tlib.PlaneTruss(macro_nodes, macro_edges, E, A)


    macro_constraints = [
        [0, 0, 0.0],
        [0, 1, 0.0],  
    ]

    P = 1e9
    Q = 0
    macro_loads = [
        [1, 0, P]
    ]


    macro_truss.apply_constraints(macro_constraints)
    macro_truss.apply_loads(macro_loads)

    # print("K before constraints")
    # print(truss.K)

    truss.apply_loads(macro_loads)
    truss.apply_constraints(macro_constraints)

    # print("K after constraints")
    # print(truss.K)

    t1 = time.time() # setup complete

    macro_truss.solve(h=0.05, n_iterations=200)
    truss.solve()

    # print(truss.K)


    t2 = time.time() # solve complete


    # macro_truss.plot(mag=1)
    truss.plot(mag=1)