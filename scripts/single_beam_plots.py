# Date: 06 June 2024
# Author: Anupsa Swain

import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb 
import single_beam_continuous as tbar
from tf_unit import u_t

if __name__ == "__main__": 

    L = 10 
    E = 7e10 # Aluminium standard
    A = 0.25 # m^2
    rho = 2710 # kg/m^3, Aluminium

    P = 1e9

    u_static = P*L/E/A

    t_array = np.linspace(0.01, .1, 100)

    freq = t_array[-1]/len(t_array)

    u_L_array = np.zeros_like(t_array)

    for i in range(len(u_L_array)):
        t = t_array[i]
        u_L_array[i] = u_t(t)

    colors = ['r', 'b']

    # for 2 layers:
    for freq in [50, 100, 150, 200]:
        layers = 1
        meshsize = 1
        conn_node_tol = 0.5
        H_meshsize = meshsize #*(3**.5)/2
        w = H_meshsize*layers

        u1_time_array = np.load('data/dynamic/new_square/l%d_m%.3f_tf_data_%gHz.npz' %(layers, meshsize, freq))['x']
        # nodes = np.load('data/dynamic/new_square/l%d_m%.3f_final_trussCord.npz' %(layers, meshsize))['x']

        # u1L_idx = np.zeros([0])

        # for i in range(len(nodes)):
        #     x0 = nodes[i, 0]
        #     y0 = nodes[i, 1]
        #     if (x0 - L)**2 + y0**2 < w**2 :
        #         u1L_idx = np.append(u1L_idx, [2*i], axis=0) # need x_axis displacement only

        # u1L_idx = u1L_idx.astype('int32')    

        # u1L_array = np.average(u1_time_array[:, u1L_idx], axis=1)

        # u1L_array = u1_time_array[0]

        # for i in range(int(len(u1_time_array)/2)-1):
        #     u1L_array += u1_time_array[2*(i+1)]

        # for i in range(len(u1_time_array)):
        #     print(np.sum(u1_time_array[i]))

        u1L_array = np.average(u1_time_array, axis=0)

        ax = sb.lineplot(x=t_array, y=u1L_array)

    # # for 4 layers:
    # layers = 4

    # u4_time_array = np.load('data/dynamic/l%d_m%.3f_displacement_data.npz' %(layers, meshsize))['x']
    # nodes = np.load('data/dynamic/l%d_m%.3f_final_trussCord.npz' %(layers, meshsize))['x']

    # u4L_idx = np.zeros([0])

    # for i in range(len(nodes)):
    #     x0 = nodes[i, 0]
    #     y0 = nodes[i, 1]
    #     if (x0 - L)**2 + y0**2 < w**2 :
    #         u4L_idx = np.append(u4L_idx, [2*i], axis=0) # need x_axis displacement only

    # u4L_idx = u4L_idx.astype('int32')    

    # u4L_array = np.average(u4_time_array[:, u4L_idx], axis=1)

    # file_name = os.path.join('data/dynamic/', 'c_2_4_end_displacement_data.npz')
    # np.savez(file_name, t=t_array, x=u_L_array, y=u2L_array, z=u4L_array)
 
    # ax1 = sb.lineplot(x=t_array, y=u_L_array)
    # ax1 = sb.lineplot(x=t_array, y=u1L_array)
    # ax4 = sb.lineplot(x=t_array, y=u4L_array)

    ax.set(xlabel = "Time, t (s)", ylabel = "Displacement, u (m)")

    # xlims = (2*np.array((ax.get_xlim())))
    # ylims = (1.5*np.array((ax.get_ylim())))

    # ax.set_xlim(xlims)
    # ax.set_ylim(ylims)

    plt.legend(("50 Hz", "100 Hz", "150 Hz", "200 Hz"), loc='best')


    # plt.title("Input Frequency = 20 Hz")

    plt.grid()
    plt.show()
    