# Script for solving for deflection in microtrusses
import numpy as np
import matplotlib.pyplot as plt 
import plane_truss_lib as tlib
import time
from numba import jit 
# from scripts.check_results_d2 import meshplot


def meshplot(trussCord, node_connect):
    fig, ax = plt.subplots()
    for element in node_connect:
        x1, x2 = trussCord[element[0]-1, 0], trussCord[element[1]-1, 0]
        y1, y2 = trussCord[element[0]-1, 1], trussCord[element[1]-1, 1]
        ax.plot([x1,x2], [y1,y2], '-k')
    x_scatter = trussCord[:,0]
    y_scatter = trussCord[:,1]
    x_scatter = np.vstack([x_scatter]).flatten()
    y_scatter = np.vstack([y_scatter]).flatten()
    plt.scatter(x_scatter, y_scatter)
    ax.set_aspect('equal')

    for i in range(len(trussCord)):
        plt.annotate(i, xy=(trussCord[i, 0], trussCord[i, 1]), c='blue')

    plt.show()
    # plt.savefig('truss.pdf')

def edge_meshplot(trussCord, node_connect):
    fig, ax = plt.subplots()
    for element in node_connect:
        x1, x2 = trussCord[element[0]-1, 0], trussCord[element[1]-1, 0]
        y1, y2 = trussCord[element[0]-1, 1], trussCord[element[1]-1, 1]
        ax.plot([x1,x2], [y1,y2], '-k')
    x_scatter = trussCord[:,0]
    y_scatter = trussCord[:,1]
    x_scatter = np.vstack([x_scatter]).flatten()
    y_scatter = np.vstack([y_scatter]).flatten()
    plt.scatter(x_scatter, y_scatter)
    ax.set_aspect('equal')

    for i in range(len(node_connect)):
        x = (trussCord[node_connect[i, 0]-1, 0] + trussCord[node_connect[i, 1]-1, 0])/2
        y = (trussCord[node_connect[i, 0]-1, 1] + trussCord[node_connect[i, 1]-1, 1])/2
        plt.annotate(i, xy=(x, y), c='blue')

    plt.show()

# within nodes, the last set are the joints. in the order defined by the macrotruss. 

def assign_constraint(macro_constraints, parent_trussCord, nodes, w):
    constraints = np.zeros([0, 3])
    for i in range(len(macro_constraints)):
        node_num = macro_constraints[i][0]

        node_location = parent_trussCord[node_num] # this is of the form [x, y], a numpy array

        # defining the circle of application
        for j in range(len(nodes)):
            node = nodes[j]

            vector = node - node_location

            if vector[0]**2 + vector[1]**2 - (w)**2 <= 0:
                constraint_element = np.array([[j, macro_constraints[i][1], macro_constraints[i][2]]])
                constraints = np.append(constraints, constraint_element, axis=0)

    constraints_list = constraints.tolist()

    return constraints_list

def assign_load(macro_loads, parent_trussCord, nodes, w):

    # breakpoint()

    loads = np.zeros([0, 3])
    for i in range(len(macro_loads)):
        node_num = macro_loads[i][0]

        node_location = parent_trussCord[node_num] # this is of the form [x, y], a numpy array
        applied_nodes_count = 0
        # defining the circle of application
        for j in range(len(nodes)):
            node = nodes[j]
            # if j > 7030: breakpoint()
            vector = node - node_location

            if vector[0]**2 + vector[1]**2 - (w)**2 <= 0:
                applied_nodes_count = applied_nodes_count +1 
                load_element = np.array([[j, macro_loads[i][1], macro_loads[i][2]]])
                loads = np.append(loads, load_element, axis=0)
        
        loads[-applied_nodes_count:, 2] = loads[-applied_nodes_count:, 2]/applied_nodes_count

    loads_list = loads.tolist()

    # print(np.sum(loads[:, 2])) # to check if forces distributed correctly

    return loads_list

if __name__ == "__main__":

    t0 = time.time()

    parent_trussCord = np.array([[0, 0],
                                    [0, 40],
                                    [40, 0],
                                    [40,40],
                                    [80, 0],
                                    [80,40]])
    
    parent_trussCord = parent_trussCord + 100

    EnConn = np.array([[1, 3],
                    [1, 4],
                    [2, 4],
                    [3, 4],
                    [3, 5],
                    [4, 5],
                    [4, 6],
                    [5, 6]])
    
    parent_edges = EnConn - 1

    # parameters
    layers = 3
    meshsize = 0.1
    H_meshsize = meshsize*3/2        # Height of the hex mesh
    w = layers*H_meshsize + meshsize/2                   # Width 

    # get microtruss nodes and edges data

    nodes_file = np.load('data/squares/l%d_m%.3f_final_trussCord.npz' %(layers, meshsize))
    edges_file = np.load('data/squares/l%d_m%.3f_final_node_connect.npz' %(layers, meshsize))
    joints_file = np.load('data/squares/l%d_m%.3f_joint_nodes_final.npz' %(layers, meshsize))

    nodes = nodes_file['x']
    edges = edges_file['x']
    joints = joints_file['x']

    nodes = np.array(nodes)
    edges = np.array(edges)
    joints = np.array(joints)

    # # # plot check:
    # meshplot(nodes, edges)

    edges = edges - 1

    # fps system
    E = 1e7 
    A = 1.5 

    truss= tlib.PlaneTruss(nodes, edges, E, A)
    macro_truss = tlib.PlaneTruss(parent_trussCord, parent_edges, E, A)

    # assigning constraints and forces

    #node number, horz/ver, displacement
    macro_constraints = [
    [0,0,0.0],
    [0,1,0.0],
    [1,0,0.0],
    [1,1,0.0]
    ]

    macro_loads =[
        [2,1,-2e3], 
        [4, 0, 2e3], 
        [5, 0, 4e3], 
        [5, 1, 6e3]
    ]

    constraints = assign_constraint(macro_constraints, parent_trussCord, nodes, w)
    truss.apply_constraints(constraints)

    loads = assign_load(macro_loads, parent_trussCord, nodes, w)       
    truss.apply_loads(loads)

    # print(w)
    # print(constraints)
    # print(loads)

    t1 = time.time()

    print("Time taken to setup = " + str(t1 - t0) + "s")

    truss.solve()

    t2 = time.time()

    print("Time taken to solve = " + str(t2 - t1) + "s")

    edge_meshplot(nodes, edges + 1)
    truss.plot()

    t3 = time.time()

    print("Time taken to plot = " + str(t3 - t2) + "s")

    macro_truss.apply_constraints(macro_constraints)
    macro_truss.apply_loads(macro_loads)

    macro_truss.solve()

    macro_truss.plot()


# solving using linalg but parsing method? check

# this should be enough and you're hopefully done *party-poppers

