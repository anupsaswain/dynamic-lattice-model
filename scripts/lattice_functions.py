import matplotlib.pyplot as plt
import numpy as np


def calc_slope(i,EnConn,parent_trussCord):
    xi, xf = parent_trussCord[EnConn[i, 0] - 1, 0], \
            parent_trussCord[EnConn[i, 1] - 1, 0]
    yi, yf = parent_trussCord[EnConn[i, 0] - 1, 1], \
            parent_trussCord[EnConn[i, 1] - 1, 1]
    if xi-xf == 0:
        slope = 'infinity'
    else:
        slope = (yf - yi)/(xf - xi)
    return slope, xi, xf, yi, yf


def centroid_correction(transform_temp,first_point,last_point):
    error_centroid_x, error_centroid_y = get_centroid_error(transform_temp,first_point,last_point)
    while abs(error_centroid_x) >= 1e-8 or abs(error_centroid_y) >= 1e-8:
        transform_temp = translate(transform_temp,error_centroid_x, error_centroid_y)
        error_centroid_x, error_centroid_y = get_centroid_error(transform_temp,first_point,last_point)
    return transform_temp


def create_vertex(modf_NeConn, node_num, counts, layers_line_const, centre_line_const, layers):
    joint_nodes = np.zeros([0,2])     
    row_min = layers + 1

# Neglecting elements connected to Vertical and Horizontal elements
    for i in range(len(node_num)):
        store_element_dir = []
        list_slope = np.array([])
        store_element_sim = []
        for j in range(counts[i]):
            element_num = modf_NeConn[j,1]
            element_slope = centre_line_const[element_num-1, 0]
            list_slope = np.append(list_slope, element_slope)
            store_element_dir.append(element_num)
# Chategorizing nodes to solve using Direct Method and Simultaneous Equation
        infinity_slope_column = np.where(list_slope == 'infinity')
        inf_check = np.size(infinity_slope_column[0])
        zero_slope_column = np.where(list_slope == '0.0')
        zero_check = np.size(zero_slope_column[0])

        if inf_check > 0 and zero_check > 0:
            first_element = store_element_dir[infinity_slope_column[0][0]]
            second_element = store_element_dir[zero_slope_column[0][0]]
            if first_element > second_element:
                temp = second_element
                second_element = first_element
                first_element = temp

            first_element_temp = int(first_element - 1)
            first_slope = centre_line_const[first_element_temp, 0]
            start_point_1 = int(first_element_temp*row_min)
            end_point_1 = start_point_1 + int(row_min)

            second_element_temp = int(second_element - 1)
            second_slope = centre_line_const[second_element_temp, 0]
            start_point_2 = int(second_element_temp*row_min)
            end_point_2 = start_point_2 + int(row_min)
            row_num_1, column_num = np.where(layers_line_const[start_point_1:end_point_1] == first_slope)
            row_num_2, column_num = np.where(layers_line_const[start_point_2:end_point_2] == second_slope)
            row_num_1 = [x + start_point_1 for x in row_num_1]
            row_num_2 = [x + start_point_2 for x in row_num_2]

            c_y_intercepts = np.zeros([0, 2])
            for k in range(row_min):
                c_y_intercept_1 = float(layers_line_const[row_num_1[k], 1])
                c_y_intercept_2 = float(layers_line_const[row_num_2[k], 1])
                c_y_intercepts = np.append(c_y_intercepts, [[c_y_intercept_1, c_y_intercept_2]], axis=0)
            C_y_intercepts = np.zeros([2, 0])
            for p in range(layers + 1):
                c_1 = c_y_intercepts[p, 0]
                for q in range(layers + 1):
                    c_2 = c_y_intercepts[q, 1]
                    C_y_intercepts = np.append(C_y_intercepts, [[c_1], [c_2]], axis=1)
            if second_slope == 'infinity':
                C_y_intercepts[[1, 0], :] = C_y_intercepts[[0, 1], :]
            for r in range(np.shape(C_y_intercepts)[1]):
                sol = np.array(C_y_intercepts[:, r],dtype=float)
                joint_nodes = np.append(joint_nodes, [sol], axis=0)
        else:
# arranging arrays in order to solve simultaneous eq
            list_slope_temp = []
            store_element_dir_temp = []
            for k in range(len(list_slope)):
                first_element = store_element_dir[k]
                least_element = np.amin(store_element_dir[k:len(store_element_dir)])
                if first_element > least_element:
                    store_element_dir_temp.append(least_element)
                    list_slope_temp.append(list_slope[k])
                else:
                    store_element_dir_temp.append(str(first_element))
                    list_slope_temp.append(list_slope[k])
            store_element_sim = np.array(store_element_dir_temp,dtype=int)
            list_slope = np.array(list_slope_temp,dtype=float)
            
            for k in range(len(store_element_dir)-1):
                for p in range(len(store_element_dir)):
                    if k < p:
                        first_element = store_element_dir[k]
                        first_element_temp = int(first_element - 1)
                        first_slope = centre_line_const[first_element_temp, 0]
                        start_point_1 = int(first_element_temp*row_min)
                        end_point_1 = start_point_1 + int(row_min)

                        second_element = store_element_dir[p]
                        second_element_temp = int(second_element - 1)
                        second_slope = centre_line_const[second_element_temp, 0]
                        start_point_2 = int(second_element_temp*row_min)
                        end_point_2 = start_point_2 + int(row_min)

                        Matrix = np.array([[first_slope, -1],
                                        [second_slope, -1]])
                        Matrix = Matrix.astype(float)
                        row_num_1, column_num = np.where(layers_line_const[start_point_1:end_point_1] == first_slope)
                        row_num_2, column_num = np.where(layers_line_const[start_point_2:end_point_2] == second_slope)
                        row_num_1 = [x + start_point_1 for x in row_num_1]
                        row_num_2 = [x + start_point_2 for x in row_num_2]

                        c_y_intercepts = np.zeros([0, 2])
                        for x in range(row_min):
                            c_y_intercept_1 = layers_line_const[row_num_1[x], 1]
                            c_y_intercept_2 = layers_line_const[row_num_2[x], 1]
                            c_y_intercepts = np.append(
                                c_y_intercepts, [[c_y_intercept_1, c_y_intercept_2]], axis=0)
                    
                    
                        C_y_intercepts = np.zeros([2, 0])
                        for y in range(row_min):
                            c_1 = c_y_intercepts[y, 0]
                            for q in range(layers + 1):
                                c_2 = c_y_intercepts[q, 1]
                                C_y_intercepts = np.append(
                                    C_y_intercepts, [[c_1], [c_2]], axis=1)
                        C_y_intercepts = C_y_intercepts.astype(float)
                        
                        for z in range(np.shape(C_y_intercepts)[1]):
                            sol = np.matmul(np.linalg.inv(Matrix), -C_y_intercepts[:, z]).reshape(1, 2)
                            joint_nodes = np.append(joint_nodes, sol, axis=0)
                

        modf_NeConn = delete_duplicate_row(np.arange(counts[i]), modf_NeConn)
  
    joint_nodes = np.array(joint_nodes,dtype=float)
    return joint_nodes


def delete_corner_nodes(trussCord_temp,c1,c2):
    duplicate_node = np.zeros([0,2])
# f and b (forward and backward) implies the direction in which nodes are generated
# along longitudinal axis
    duplicate_node_f = np.zeros([0,2])
    duplicate_node_b = np.zeros([0,2])
    if len(trussCord_temp) < 100:
        iter = int(len(trussCord_temp)*0.5)
    else:
        iter = int(len(trussCord_temp)*0.3)
    for m in range(2):              # one element is connected to only 2 nodes
        if m == 0 & np.size(c1) != 0:
            for n in range(iter):
                micro_Cord = trussCord_temp[n]
# c1[0] is the cordinate point farther from centroid of truss element than c1[1].
                dis1 = dis_btw_points(micro_Cord, c1[0])
                dis2 = dis_btw_points(micro_Cord, c1[1])
                if dis1 > dis2:
                    duplicate_node_f = np.append(duplicate_node_f, [micro_Cord],axis=0)
                else:
                    break
        elif m == 1 & np.size(c2) != 0:
            for n in range(len(trussCord_temp)-1, iter+1, -1):
                micro_Cord = trussCord_temp[n]
# c2[0] is the cordinate point farther from centroid of truss element than c2[1].
                dis1 = dis_btw_points(micro_Cord, c2[0])
                dis2 = dis_btw_points(micro_Cord, c2[1])
                if dis1 > dis2:
                    duplicate_node_f = np.append(duplicate_node_f, [micro_Cord],axis=0)
                else:
                    break
    duplicate_node = np.append(duplicate_node, duplicate_node_f,axis=0)
    duplicate_node = np.append(duplicate_node, duplicate_node_b,axis=0)
    print('duplicate_node',len(duplicate_node))
    duplicate_row = list_duplicate_node(duplicate_node, trussCord_temp)
    trussCord_temp = delete_duplicate_row(duplicate_row, trussCord_temp)
    return trussCord_temp


def delete_corner_points(
        up_bound_corr_factor,
        layers_line_const,
        layers,
        meshsize,
        centre_line_const,
        trussCord):
    up_bound_corr = up_bound_corr_factor * meshsize
    corner_node = np.zeros([0, 2])
    c_initz_max = 1e9
    c_initz_min = -1e9
    out_line_const = np.zeros([0, 2])
    lines_per_strut = int(len(layers_line_const) / len(centre_line_const))
    for i in range(len(centre_line_const)):
        row_jump = i * lines_per_strut
        slope = centre_line_const[i, 0]

        if slope > 0.0:
            c_max = c_initz_min
            for j in range(layers + 1):
                c = layers_line_const[row_jump + j, 1]
                if c > c_max:
                    c_max = c
            out_line_const = np.append(
                out_line_const, [[slope, c_max]], axis=0)
            print('slope>0', c_max)
            for k in range(len(trussCord)):
                x, y = trussCord[k, 0], trussCord[k, 1]
                if y - slope * x - c_max - up_bound_corr > 0:
                    corner_node = np.append(corner_node, [[x, y]], axis=0)
        elif slope < 0.0:
            c_max = c_initz_min
            for j in range(layers + 1):
                c = layers_line_const[row_jump + j, 1]
                if c > c_max:
                    c_max = c
            out_line_const = np.append(
                out_line_const, [[slope, c_max]], axis=0)
            print('slope<0', c_max)
            for k in range(len(trussCord)):
                x, y = trussCord[k, 0], trussCord[k, 1]
                if y - slope * x - c_max - up_bound_corr > 0:
                    corner_node = np.append(corner_node, [[x, y]], axis=0)
        else:
            c_min = c_initz_max
            for j in range(layers + 1):
                c = layers_line_const[row_jump + j, 1]
                if c < c_min:
                    c_min = c
            out_line_const = np.append(
                out_line_const, [[slope, c_min]], axis=0)
            print('slope=0', c_min)
            for k in range(len(trussCord)):
                x, y = trussCord[k, 0], trussCord[k, 1]
                if y - slope * x - c_min < 0:
                    corner_node = np.append(corner_node, [[x, y]], axis=0)
    print("len(corner_node)", len(corner_node))
    return corner_node, out_line_const


def delete_duplicate_row(duplicate_row, trussCord):
    trussCord = np.delete(trussCord, duplicate_row, 0)
    return trussCord


def dis_btw_points(iniPoint, finPoint):
    dis = np.sqrt( (finPoint[0]-iniPoint[0])**2 + (finPoint[1]-iniPoint[1])**2 )
    return dis


def far_near(Cord_array,x,y):
    dis_max = 0.0
    dis_min = 1e9
    for i in range(len(Cord_array)):
        dis = dis_btw_points(Cord_array[i], [x, y])
        if dis > dis_max:
            dis_max = dis
            far_node = Cord_array[i]
        if dis < dis_min:
            dis_min = dis
            near_node = Cord_array[i]
    return far_node, near_node


def gen_conn_matrix(parent_trussCord, trussCord, meshsize, conn_node_tol, w, tol_app_radius):
    # generation of micro-structure connectivity matrix
    node_connect = np.zeros([0, 2])
    for i in range(len(trussCord)):
        start_point = trussCord[i]
        for j in range(len(trussCord)):
            if i < j:
                end_point = trussCord[j]
                dis = dis_btw_points(start_point, end_point)
                if dis < (1 + conn_node_tol[0]) * meshsize:
                    # [[i+1,j+1]] means python indexing is taken care of
                    node_connect = np.append(
                        node_connect, [[i + 1, j + 1]], axis=0)
                else:
                    for k in range(len(parent_trussCord)):
                        joint = parent_trussCord[k]
                        rad1 = dis_btw_points(joint, start_point)
                        rad2 = dis_btw_points(joint, end_point)
                        if ((rad1 < tol_app_radius[k, 1]) & (rad2 < tol_app_radius[k, 1])):
                            if ((rad1 >= tol_app_radius[k, 0]) & (rad2 >= tol_app_radius[k, 0])):
                                if dis < (1 + conn_node_tol[k+1])*meshsize: # k+1 since [0] is tolerance for main body
                                    node_connect = np.append(
                                        node_connect, [[i+1, j+1]], axis=0)
    node_connect = node_connect.astype(int)
    return node_connect


def gen_HzStrut_cord(dis,meshsize,w,layers):
    # iteration for layer generation
    itr = layers + 1

    INCR_dis = meshsize
    numPoints = dis/meshsize + 2 # *(layers)       # No. of points to generate in longitudinal direction
    HzStrut_cord = np.zeros([0,2])

    # generation of mesh coordinates
    for i in range(int(numPoints)):
        for j in range(int(itr)):
            x = -2*w + i*INCR_dis   
            y = -w*0.5 + j*meshsize
            HzStrut_cord = np.append(HzStrut_cord, [[x,y]],axis=0)

    return HzStrut_cord, numPoints 


def gen_layers_line(centre_line_const, trussCord, layers, points_per_strut):
    layers_line_const = np.zeros([0, 2])
    print('centre_line_const',centre_line_const)
    tot_centre_line = len(centre_line_const)  # number of strut
    row_jump = 0
    for i in range(tot_centre_line):
        row_jump += points_per_strut[i]
        a_slope = centre_line_const[i, 0]
        if a_slope == 'infinity':
            for j in range(layers +1):  # j represents y intercept for layers
                if i == 0:
                    row = j
                    a_point = trussCord[row]
                    a_x_intercept_C = a_point[0]
                    layers_line_const = np.append(layers_line_const, [['infinity', a_x_intercept_C]], axis=0)
                else:
                    row = row_jump + j #+ 1
                    a_point = trussCord[row]
                    a_x_intercept_C = a_point[0]
                    layers_line_const = np.append(layers_line_const, [['infinity', a_x_intercept_C]], axis=0)
        else:
            for j in range(layers +1):  # j represents y intercept for layers
                if i == 0:
                    row = j
                    a_point = trussCord[row]
                    a_y_intercept_C = a_point[1] - float(a_slope) * a_point[0]
                    layers_line_const = np.append(layers_line_const, [[a_slope, a_y_intercept_C]], axis=0)
                else:
                    row = row_jump + j #+ 1
                    a_point = trussCord[row]
                    a_y_intercept_C = a_point[1] - float(a_slope) * a_point[0]
                    layers_line_const = np.append(layers_line_const, [[a_slope, a_y_intercept_C]], axis=0)
    return layers_line_const


def gen_NeConn(EnConn):
    NeConn = np.zeros([0, 2])
    Searched_rows = np.array([])
    for i in range(len(EnConn)):
        for j in range(2):
            if EnConn[i, j] not in Searched_rows:
                row_num, column_num = np.where(EnConn == EnConn[i, j])
                for k in range(len(row_num)):
                    NeConn = np.append(NeConn, [[EnConn[i, j], row_num[k]+1]], axis=0)                
                Searched_rows = np.append(Searched_rows, [EnConn[i, j]], axis=0)
    return NeConn.astype(int)


def gen_parent_line_eq(EnConn,parent_trussCord):
    centre_line_const = np.zeros([0, 2])
    for i in range(len(EnConn)):
        start_node = parent_trussCord[EnConn[i,0]-1]
        end_node = parent_trussCord[EnConn[i,1]-1]
        if end_node[0] == start_node[0]:
            slope_m = 'infinity'
            y_intercept_C = 'no'
        else:
            slope_m = (end_node[1] - start_node[1]) / (end_node[0] - start_node[0])
            y_intercept_C = start_node[1] - slope_m * start_node[0]
        centre_line_const = np.append(centre_line_const, [[slope_m, y_intercept_C]], axis=0)
    return centre_line_const


def get_centroid_error(transform_temp,first_point,last_point):
    # Calculate the average centroid of the coordinates
    sum_x, sum_y = 0.0, 0.0
    for i in range(len(transform_temp)):
        sum_x += transform_temp[i,0]
        sum_y += transform_temp[i,1]
    avg_x = sum_x/(i+1)
    avg_y = sum_y/(i+1)

    true_centroid_x, true_centroid_y = solid_strut_centroid(first_point, last_point)

    error_centroid_x = true_centroid_x - avg_x
    error_centroid_y = true_centroid_y - avg_y
    return error_centroid_x, error_centroid_y


def get_modf_NeConn(NeConn):
# Neglecting idle node connected to only one element
    node_num, counts = np.unique(NeConn[:,0], return_counts=True)
    counts_column = np.where(counts == 1)
    NeConn_row = np.where(NeConn[:,0] == node_num[counts_column[0]])
    modf_NeConn = delete_duplicate_row(NeConn_row[0], NeConn)
    node_num = delete_duplicate_row(counts_column[0], node_num)
    counts = delete_duplicate_row(counts_column[0], counts)
    return modf_NeConn, node_num, counts


def gen_vertex(
        HzStrut_cord,
        parent_trussCord,
        meshsize,
        centre_line_const,
        trussCord,
        layers,
        NeConn,
        delete_mesh_tol,
        origin_tol,
        up_bound_corr_factor,
        layers_line_const):
    say = 0

    if say == 1:
        if layers != 1:
            duplicate_node = id_duplicate_nodes(
                meshsize, trussCord, delete_mesh_tol)
            duplicate_row = list_duplicate_node(duplicate_node, trussCord)
            trussCord = delete_duplicate_row(duplicate_row, trussCord)

            corner_node, out_line_const = delete_corner_points(
                up_bound_corr_factor, layers_line_const, layers, meshsize, centre_line_const, trussCord)
            duplicate_row = list_duplicate_node(corner_node, trussCord)
            trussCord = delete_duplicate_row(duplicate_row, trussCord)
    else:
        trussCord = create_vertex(
            trussCord,
            NeConn,
            layers_line_const,
            centre_line_const, layers)
    return trussCord


def id_duplicate_nodes(meshsize,trussCord,delete_mesh_tol):
    duplicate_node = np.zeros([0,2])
    for i in reversed(range(len(trussCord))):
        start_point = trussCord[i]
        for j in reversed(range(len(trussCord))):
            if i > j:
                end_point = trussCord[j]
                dis = dis_btw_points(start_point, end_point)
                if dis <= delete_mesh_tol*meshsize:
                    duplicate_node = np.append(duplicate_node, [start_point, end_point],axis=0)
    return duplicate_node


def list_duplicate_node(duplicate_node, trussCord):
    arr = trussCord
    ant = duplicate_node
    duplicate_row = np.array([])
    for i in range(len(duplicate_node)):
        result = np.where((arr[:, 0] == ant[i, 0]) & (arr[:, 1] == ant[i, 1]))
        if len(result[0]) > 1:
            result = result[0]
            print("result_1",result)
            result = result[:-1]
            print("result",result)
            duplicate_row = np.append(duplicate_row, result, axis=0)
        else:
            duplicate_row = np.append(duplicate_row, result[0], axis=0)
    duplicate_row = duplicate_row.astype(int)
    return duplicate_row


def id_dup_node_Cord(duplicate_node, start_point,end_point,meshsize,delete_mesh_tol):
    dis = dis_btw_points(start_point, end_point)
    if dis <= delete_mesh_tol*meshsize:
        duplicate_node = np.append(duplicate_node, [end_point],axis=0)
    return duplicate_node


def delete_dup_vertex(trussCord_temp, joint_nodes, meshsize, joint_mesh_tol):
    # for intersection nodes
    duplicate_node = np.zeros([0,2])
    for j in range(len(joint_nodes)):
        start_point = joint_nodes[j]
        for k in range(len(trussCord_temp)):
            end_point = trussCord_temp[k]
            dis = dis_btw_points(start_point, end_point)
            if dis <= joint_mesh_tol*meshsize:
                duplicate_node = np.append(duplicate_node, [end_point],axis=0)
    duplicate_row = list_duplicate_node(duplicate_node, trussCord_temp)
    trussCord_temp = delete_duplicate_row(duplicate_row, trussCord_temp)
    return trussCord_temp


def find_short_dis(start_point,parent_trussCord):
    dis_temp = 1e9
    node_numb = 0
    for i in range(len(parent_trussCord)):
        dis = dis_btw_points(start_point, parent_trussCord[i])
        if dis <= dis_temp:
            dis_temp = dis
            node_numb = i
    return dis_temp, node_numb


def get_trussCord(start_points,end_points,meshsize,delete_INmesh_tol,delete_OUTmesh_tol,parent_trussCord, w, upper_threshold, lower_threshold):
    duplicate_node = np.zeros([0,2])
    for j in range(len(start_points)):
        start_point = start_points[j]                       #Cordinates of j th element
        start_point_dis, node_1 = find_short_dis(start_point,parent_trussCord)      # find the shortest distance between 'start_point' and vertex nodes

        for k in range(len(end_points)):
            end_point = end_points[k]                       #Cordinates of k th element
            end_point_dis, node_2 = find_short_dis(end_point,parent_trussCord)

            if start_point_dis**2 < 1e18 and end_point_dis**2 < 1e18 and end_point_dis < upper_threshold*w and node_1 == node_2: #  and start_point_dis**2 > end_point_dis**2
# Delete nodes from k th element closer to parent_trussCord
                dis = dis_btw_points(start_point, end_point)
                if dis**2 <= (delete_INmesh_tol*meshsize)**2:
                    duplicate_node = np.append(duplicate_node, [end_point],axis=0)
            if start_point_dis**2 < 1e18 and end_point_dis**2 < 1e18 and end_point_dis > lower_threshold*w and node_1 == node_2:
# Delete nodes from k th element farther to parent_trussCord
                dis = dis_btw_points(start_point, end_point)
                if dis**2 <= (delete_OUTmesh_tol*meshsize)**2:
                    duplicate_node = np.append(duplicate_node, [end_point],axis=0)
 
    print('len(duplicate_node)', len(duplicate_node))

    if len(duplicate_node) != 0:
        duplicate_row_2 = list_duplicate_node(duplicate_node, end_points)
        end_points = delete_duplicate_row(duplicate_row_2, end_points)
        print('yes')
    return start_points, end_points


def mid_point(first_pnt,second_pnt):
    return [0.5*(first_pnt[0]+second_pnt[0]), 0.5*(first_pnt[1]+second_pnt[1])]


def my_plot(trussCord):
    x_scatter = trussCord[:,0]
    y_scatter = trussCord[:,1]
    x_scatter = np.vstack([x_scatter]).flatten()
    y_scatter = np.vstack([y_scatter]).flatten()
    plt.figure(figsize=(50,50))
    plt.scatter(x_scatter, y_scatter, s = 1)

    plt.grid()
    plt.show()
    

def reflect_HzStrut_cord(HzStrut_cord,ref_vector):
    # Look Wikipedia for the reflection formula
    reflect_M = 1/(ref_vector[1]**2 + ref_vector[0]**2) * np.array([[ ref_vector[0]**2 -ref_vector[1]**2, 2*ref_vector[0]*ref_vector[1]],
                                                                    [ 2*ref_vector[0]*ref_vector[1], ref_vector[1]**2- ref_vector[0]**2] ])
    reflect_temp = np.zeros([0,2])
    for i in range(len(HzStrut_cord)):
        matrix = HzStrut_cord[i]
        cord_dash = np.dot(reflect_M,matrix)
        reflect_temp = np.append(reflect_temp, [[cord_dash[0],cord_dash[1]]], axis=0)
    return reflect_temp


def solid_strut_centroid(first_point, last_point):
    true_centroid_x = 0.5*(first_point[0]+last_point[0])
    true_centroid_y = 0.5*(first_point[1]+last_point[1])
    return true_centroid_x, true_centroid_y


def transform_HzStrut_cord(HzStrut_cord,slope,layers):
    if slope == 'infinity':
        angle = np.pi/2
    else:
        angle = np.arctan(-slope) # Bcz of matrix multiplication issue, I am getting negative slope. Hence this modification

    if (layers % 2) != 0:
        HzStrut_cord = reflect_HzStrut_cord(HzStrut_cord,[1,0])# reflect about ref_vector
    transform_temp = np.zeros([0,2])
    # Rotate coordinates
    cord_trans = np.array([[np.cos(angle), np.sin(angle)],
                    [-np.sin(angle), np.cos(angle)]])
    for i in range(len(HzStrut_cord)):
        matrix = HzStrut_cord[i]
        cord_dash = np.dot(cord_trans,matrix)
        transform_temp = np.append(transform_temp, [[cord_dash[0],cord_dash[1]]], axis=0)
    return transform_temp


def translate(matrix,x, y):
    matrix = matrix + np.array([x,y])
    return matrix