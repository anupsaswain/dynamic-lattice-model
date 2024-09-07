# Author: Anupsa Swain
# Date: 23 May 2024

import numpy as np
import matplotlib.pyplot as plt

# q is a 12x1 matrix
# M, C, K will be 12x12 matrices. will need setting up. 
# set h = 1 second for now
# everything is solved in SI units

def force(time):
    forces = np.zeros(12)

    forces[11] = np.sin(np.pi*time/6)

    return forces

def newmark_integration(nodes, edges, M, C, K, q0_dot, p, h, n_iterations, gamma=0.5, beta=0.25):
    t = 0
    q0 = np.concatenate(np.array([nodes[i] for i in range(len(nodes))]))
    p0 = p(0) # p is force, function of time; gives force in all 12 dims
    q0_doubledot = np.linalg.inv(M)@(p0 - C@q0_dot - K@q0)
    m, n = np.shape(nodes)
    dim = m*n
    qs = np.zeros([0, dim])
    q_dots = np.zeros([0, dim])
    q_doubledots = np.zeros([0, dim])

    qs = np.append(qs, [q0], axis=0)
    q_dots = np.append(q_dots, [q0_dot], axis=0)
    q_doubledots = np.append(q_doubledots, [q0_doubledot], axis=0)

    for i in range(n_iterations):
        t = t + h
        q_dot_star = q_dots[i] + (1-gamma)*h*q_doubledots[i]
        q_star = qs[i] + h*q_dots[i] + (0.5 - beta)*h*h*q_doubledots[i]

        S = M + h*gamma*C + h*h*beta*K

        qn_doubledot = np.linalg.inv(S)@(p(t) - C@q_dot_star - K@q_star)
        qn_dot = q_dot_star + h*gamma*qn_doubledot
        qn = q_star + h*h*beta*qn_doubledot

        qs = np.append(qs, [qn], axis=0)
        q_dots = np.append(q_dots, [qn_dot], axis=0)
        q_doubledots = np.append(q_doubledots, [qn_doubledot], axis=0)

    return qs, q_dots, q_doubledots



if __name__ == "__main__":
    h = 1 # seconds
    n_iterations = 20
    nodes = np.array([[0, 0],
                                [0, 40],
                                [40, 0],
                                [40,40],
                                [80, 0],
                                [80,40]])
    nodes = nodes + 100*np.ones(np.shape(nodes))

    edges = np.array([[1, 3],
                    [1, 4],
                    [2, 4],
                    [3, 4],
                    [3, 5],
                    [4, 5],
                    [4, 6],
                    [5, 6]])
    
    edges = edges - 1 # no more +1s hereon

    # need to define M, C, K, q0_dot

    

    # qs, q_dots, q_doubledots = newmark_integration(nodes, edges, M, C, K, q0_dot, force, h, n_iterations) # gamma, beta not added here