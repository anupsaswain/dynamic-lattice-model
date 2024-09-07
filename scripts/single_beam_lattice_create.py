# Author: Anupsa Swain
# Date: 28 May 2024

import numpy as np
import matplotlib.pyplot as plt

# Generalized code for triangular micro-architectuired truss

import argparse
from argparse import ArgumentParser
import os
import datetime
import time

from lattice_functions import *
# from check_results_d2 import meshplot

def meshplot(trussCord, node_connect):
    fig, ax = plt.subplots()
    i = 0
    for element in node_connect:
        x1, x2 = trussCord[element[0]-1, 0], trussCord[element[1]-1, 0]
        y1, y2 = trussCord[element[0]-1, 1], trussCord[element[1]-1, 1]
        ax.plot([x1,x2], [y1,y2], '-k')
        # plt.annotate(i, xy=((x1+x2)/2, (y1+y2)/2), c='red')
        i = i + 1
    x_scatter = trussCord[:,0]
    y_scatter = trussCord[:,1]
    x_scatter = np.vstack([x_scatter]).flatten()
    y_scatter = np.vstack([y_scatter]).flatten()
    plt.scatter(x_scatter, y_scatter)
    ax.set_aspect('equal')
    ax.set_xlim([-2, 12])
    ax.set_ylim([-2, 2])

    # for i in range(len(trussCord)):
    #     plt.annotate(i, xy=(trussCord[i, 0], trussCord[i, 1]), c='blue')

    # for i in range(len(node_connect)):
    #     plt.annotate(i, xy=(node_connect[i, 0], node_connect[i, 1]), c='red')

    plt.show()
    # plt.savefig('truss.pdf')

def main():
    p = ArgumentParser(description='Micro-architectured triangular truss')
    p.add_argument(
        '-n',
        type=int,
        default=100,
        help='number of pixels for the image to produce along the x-axis')
    args = p.parse_args()

    main_truss_code()


def main_truss_code():
    Tn = time.perf_counter()
    L = 10
    parent_trussCord = np.array([[0, 0],
                                [L, 0]])
    # parent_trussCord = parent_trussCord + 100*np.ones(np.shape(parent_trussCord))
    print('parent_trussCord\n',parent_trussCord)
    EnConn = np.array([[1, 2]])
    
    # parameters
    layers = 3
    meshsize = 1
    propEA = [1e6, 0.005]
    joint_mesh_tol = 0.4
    delete_INmesh_tol = 0.5                  # Two nodes < delete_mesh_tol*mesh, apart is treated as duplicate. Nodes close to parent_trussCord
    delete_OUTmesh_tol = 0.5                  # Two nodes < delete_mesh_tol*mesh, apart is treated as duplicate. Nodes farther to parent_trussCord
    conn_node_tol = 0.3         # used in check_results_d2.py
    upper_threshold = 1.1
    lower_threshold = 1.05
    print('layers', layers)
    print('meshsize', meshsize)
    print('delete_mesh_tol',delete_INmesh_tol)
    print('delete_mesh_tol',delete_OUTmesh_tol)


    H_meshsize = meshsize       # Height of the triangular mesh
    w = layers*H_meshsize                   # Width
    trussCord = np.zeros([0,2])             # Initialization of micro-architecture coordinate matrix
    points_per_strut = [0]                  # Number of nodes per strut. Initializing for computationalpurpose only.
    


# Mesh Generation
    for i in range(len(EnConn)):
        slope, xi, xf, yi, yf = calc_slope(i,EnConn,parent_trussCord)
        dis = dis_btw_points([xi,yi],[xf,yf])
        HzStrut_cord, numPoints = gen_HzStrut_cord(dis,meshsize,w,layers)
        transform_temp = transform_HzStrut_cord(HzStrut_cord,slope,layers)
        transform_temp = centroid_correction(transform_temp,[xi,yi],[xf,yf])
        points_per_strut.append(len(transform_temp))
        file_name = os.path.join('data/dynamic/new_square/', 'l%d_m%.3f_trussCord%d.npz' %(layers, meshsize, i+1))
        np.savez(file_name, x=transform_temp)        # Trusscord array for Individual elements
        trussCord = np.append(trussCord, transform_temp,axis=0)
# Generation of truss nodes COMPLETE!


    centre_line_const = gen_parent_line_eq(EnConn,parent_trussCord)
    layers_line_const = gen_layers_line(centre_line_const, trussCord, layers, points_per_strut)

# Neglecting idle node connected to only one element
    NeConn = gen_NeConn(EnConn) 
    modf_NeConn, node_num, counts = get_modf_NeConn(NeConn)
    print('modf_NeConn',modf_NeConn)
    print('node_num', node_num)
    print('counts',counts)
    joint_nodes = create_vertex(modf_NeConn, node_num, counts, layers_line_const, centre_line_const, layers)

# separating joint nodes and creating dummy nodes for deletion of duplicate trussCord
    start_row = 0
# joint_nodes_extra = this is a dummy nodes to delete extra nodes
# joint_nodes_final = set of all joint nodes
    joint_nodes_extra = np.zeros([0,2])
    joint_nodes_final = np.zeros([0,2])
    row = layers + 1
    for i, item in enumerate(node_num):
        end_row = start_row + (layers+1)*(layers+1)
        joint_nodes_temp = joint_nodes[start_row:end_row,:]
        joint_nodes_temp = joint_nodes_temp.tolist()
        joint_nodes_temp.sort(key=lambda joint_nodes_temp:joint_nodes_temp[0])
        joint_nodes_temp.sort(key=lambda joint_nodes_temp:joint_nodes_temp[1])
        joint_nodes_temp = np.array(joint_nodes_temp)
        globals()[f"joint_nodes{item}"] = joint_nodes_temp
        joint_nodes_final = np.append(joint_nodes_final, joint_nodes_temp, axis=0)
        row_jump = 0
        for j in range(int(layers*layers)):
            if (j%layers) == 0 and j>0:
                row_jump += 1
            if (j%layers) == layers-1:
                a = joint_nodes_temp[row_jump+j,:]
                b = joint_nodes_temp[row_jump+j+1,:]
                c = joint_nodes_temp[row_jump+row+j,:]
                d = joint_nodes_temp[row_jump+row+j+1,:]
                x_new, y_new = 0.5*(b[0]+d[0]), 0.5*(b[1]+d[1])
                joint_nodes_extra = np.append(joint_nodes_extra, [[x_new, y_new]], axis=0)
            if row_jump == layers-1:
                a = joint_nodes_temp[row_jump+j,:]
                b = joint_nodes_temp[row_jump+j+1,:]
                c = joint_nodes_temp[row_jump+row+j,:]
                d = joint_nodes_temp[row_jump+row+j+1,:]
                x_new, y_new = 0.5*(c[0]+d[0]), 0.5*(c[1]+d[1])
                joint_nodes_extra = np.append(joint_nodes_extra, [[x_new, y_new]], axis=0)
            
            a = joint_nodes_temp[row_jump+j,:]
            b = joint_nodes_temp[row_jump+j+1,:]
            c = joint_nodes_temp[row_jump+row+j,:]
            d = joint_nodes_temp[row_jump+row+j+1,:]
            x_new, y_new = 0.25*(a[0]+b[0]+c[0]+d[0]), 0.25*(a[1]+b[1]+c[1]+d[1])
            joint_nodes_extra = np.append(joint_nodes_extra, [[x_new, y_new]], axis=0)

            x_new, y_new = 0.5*(a[0]+b[0]), 0.5*(a[1]+b[1])
            joint_nodes_extra = np.append(joint_nodes_extra, [[x_new, y_new]], axis=0)

            x_new, y_new = 0.5*(a[0]+c[0]), 0.5*(a[1]+c[1])
            joint_nodes_extra = np.append(joint_nodes_extra, [[x_new, y_new]], axis=0)
        print('size_joint_nodes_extra', np.shape(joint_nodes_extra))
        start_row = end_row

    file_name = os.path.join('data/dynamic/new_square/', 'l%d_m%.3f_joint_nodes_final.npz' %(layers, meshsize))
    np.savez(file_name, x=joint_nodes_final)

    file_name = os.path.join('data/dynamic/new_square/', 'l%d_m%.3f_joint_nodes_extra.npz' %(layers, meshsize))
    np.savez(file_name, x=joint_nodes_extra)
# Generation of joint nodes, extra joint nodes COMPLETE!


# Remove redundant nodes 
# Organizing mesh
    joint_nodes = np.append(joint_nodes_final, joint_nodes_extra, axis=0)
    for i in range(len(EnConn)):      # scan over parentstrut elements
        trussCord_temp = np.load('data/dynamic/new_square/l%d_m%.3f_trussCord%d.npz' %(layers, meshsize, i+1))
        trussCord_temp = trussCord_temp['x']
# deleting intersection nodes with joint nodes
        trussCord_temp = delete_dup_vertex(trussCord_temp, joint_nodes, meshsize, joint_mesh_tol)
        file_name = os.path.join('data/dynamic/new_square/', 'l%d_m%.3f_trussCord%d.npz' %(layers, meshsize, i+1))
        np.savez(file_name, x=trussCord_temp)


# deleting intersection nodes with macro-strut nodes
    now = datetime.datetime.now()
    print(now)
    for i in range(len(node_num)):
        node = node_num[i]
        row = np.where(modf_NeConn[:,0] == node)
        elements = modf_NeConn[row[0],1]
        for j in range(len(elements)):
            element_1 = elements[j]
            slope_element_1 = centre_line_const[element_1-1,0]
            for k in range(j+1,len(elements)):
                element_2 = elements[k]
                slope_element_2 = centre_line_const[element_2-1,0]
                if slope_element_1 != slope_element_2:
                    element_1_nodes = np.load('data/dynamic/new_square/l%d_m%.3f_trussCord%d.npz' %(layers, meshsize, element_1))
                    element_1_nodes = element_1_nodes['x']
                    element_2_nodes = np.load('data/dynamic/new_square/l%d_m%.3f_trussCord%d.npz' %(layers, meshsize, element_2))
                    element_2_nodes = element_2_nodes['x']
                    if slope_element_1 == '0.0' or slope_element_1 == 'infinity':
                        print('element_1',element_1)
                        print('element_2',element_2)
                        element_1_nodes, element_2_nodes = get_trussCord(element_1_nodes,element_2_nodes,meshsize,delete_INmesh_tol,delete_OUTmesh_tol,parent_trussCord, w, upper_threshold, lower_threshold)

                        file_name = os.path.join('data/dynamic/new_square/', 'l%d_m%.3f_trussCord%d.npz' %(layers, meshsize, element_1))
                        np.savez(file_name, x=element_1_nodes)
                        file_name = os.path.join('data/dynamic/new_square/', 'l%d_m%.3f_trussCord%d.npz' %(layers, meshsize, element_2))
                        np.savez(file_name, x=element_2_nodes)

                        now = datetime.datetime.now()
                        print(now)
                    else:
                        print('element_1 initial',element_1)
                        print('element_2 initial',element_2)
                        temp = element_1
                        element_1 = element_2
                        element_2 = temp
                        print('element_1',element_1)
                        print('element_2',element_2)

                        temp_nodes = element_1_nodes
                        element_1_nodes = element_2_nodes
                        element_2_nodes = temp_nodes

                        element_1_nodes, element_2_nodes = get_trussCord(element_1_nodes,element_2_nodes,meshsize,delete_INmesh_tol,delete_OUTmesh_tol,parent_trussCord, w, upper_threshold, lower_threshold)

                        file_name = os.path.join('data/dynamic/new_square/', 'l%d_m%.3f_trussCord%d.npz' %(layers, meshsize, element_1))
                        np.savez(file_name, x=element_1_nodes)
                        file_name = os.path.join('data/dynamic/new_square/', 'l%d_m%.3f_trussCord%d.npz' %(layers, meshsize, element_2))
                        np.savez(file_name, x=element_2_nodes)

                        now = datetime.datetime.now()
                        print(now)
                        temp = element_2
                        element_2 = element_1
                        element_1 = temp



    file_name = os.path.join('data/dynamic/new_square/', 'l%d_m%.3f_trussCord_final.npz' %(layers, meshsize))
    np.savez(file_name, x=trussCord)


    trussCord = np.zeros([0,2])
    for i in range(len(node_num)):
        trussCord_file = np.load('data/dynamic/new_square/l%d_m%.3f_trussCord%d.npz' %(layers, meshsize, i+1))
        trussCord = np.append(trussCord, trussCord_file['x'], axis=0)
    joint_nodes_file = np.load('data/dynamic/new_square/l%d_m%.3f_joint_nodes_final.npz' %(layers, meshsize))
    trussCord = np.append(trussCord, joint_nodes_file['x'], axis=0)
    file_name = os.path.join('data/dynamic/new_square/', 'l%d_m%.3f_trussCord_final.npz' %(layers, meshsize))
    np.savez(file_name, x=trussCord)

    
    # node_connect = gen_conn_matrix(trussCord, meshsize, conn_node_tol)
    

    print('Creation Complete')

def gen_conn_matrix(trussCord, meshsize, conn_node_tol):
    # generation of micro-structure connectivity matrix
    node_connect = np.zeros([0, 2])
    for i in range(len(trussCord)):
        start_point = trussCord[i]
        for j in range(len(trussCord)):
            if i < j:
                end_point = trussCord[j]
                dis = dis_btw_points(start_point, end_point)
                if dis < (1 + conn_node_tol) * meshsize:
                    # [[i+1,j+1]] means python indexing is taken care of
                    node_connect = np.append(
                        node_connect, [[i + 1, j + 1]], axis=0)
    node_connect = node_connect.astype(int)
    return node_connect


if __name__ == '__main__':
    main()

    # params
    layers = 1
    meshsize = 1
    conn_node_tol = 0.5

    joint_nodes_file = np.load('data/dynamic/new_square/l%d_m%.3f_joint_nodes_final.npz' %(layers, meshsize))
    joint_nodes_extra_file = np.load('data/dynamic/new_square/l%d_m%.3f_joint_nodes_extra.npz' %(layers, meshsize))
    # trussCord_file = np.load('data/dynamic/l%d_m%.3f_trussCord_final.npz' %(layers, meshsize))


    trussCord = np.zeros([0,2])
    for i in range(1):
        trussCord_file = np.load('data/dynamic/new_square/l%d_m%.3f_trussCord%d.npz' %(layers, meshsize, i+1))
        trussCord = np.append(trussCord, trussCord_file['x'], axis=0)

    trussCord = np.append(trussCord, joint_nodes_file['x'], axis=0)

    # # deleting points in trussCord; layers = 2 --------------------------------------------------------
    # trussCord = np.delete(trussCord, np.array([8139, 2873]), 0)
# '''    
    node_connect = gen_conn_matrix(trussCord, meshsize, conn_node_tol)

    # # adding lines; layers = 3 -----------------------
    # node_connect = np.append(node_connect, np.array([[5421, 10833], [8593, 10833], [3822, 5424], [3822, 14034], [10826, 13995], [8591, 13999], [1584, 14014], [1592, 14017], [1, 3], [3, 14009]]), 0)

    # saving trussCord:
    file_name = os.path.join('data/dynamic/new_square/', 'l%d_m%.3f_final_trussCord.npz' %(layers, meshsize))
    np.savez(file_name, x=trussCord)

    # saving node_connect:
    file_name = os.path.join('data/dynamic/new_square/', 'l%d_m%.3f_final_node_connect.npz' %(layers, meshsize))
    np.savez(file_name, x=node_connect)

    # plotting:
    meshplot(trussCord, node_connect)

# '''