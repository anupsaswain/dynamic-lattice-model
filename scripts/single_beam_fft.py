# Date: 11 June 2024
# Author: Anupsa Swain

import numpy as np
from scipy.fft import fft, ifft, fftfreq
import matplotlib.pyplot as plt
import single_beam_continuous as tbar

if __name__ == "__main__": 

    L = 10 
    E = 7e10 # Aluminium standard
    A = 0.25 # m^2
    rho = 2710 # kg/m^3, Aluminium

    P = 1e9

    u_static = P*L/E/A

    t_array = np.linspace(0, .1, 100)

    # freq = t_array[-1]/len(t_array)

    u_L_array = np.zeros_like(t_array)

    for i in range(len(u_L_array)):
        t = t_array[i]
        u_L_array[i] = tbar.u(P, L, E, A, rho, L, t) #*u_static

    # for 2 layers:
    for freq in np.linspace(10, 800, 80):
        layers = 1
        meshsize = .1
        conn_node_tol = 0.5
        H_meshsize = meshsize #*(3**.5)/2
        w = H_meshsize*layers

# (50/freq)**1.5*
        u2_time_array = np.load('data/dynamic/new_square/laplace/l%d_m%.3f_fft_data_f%gHz.npz' %(layers, meshsize, freq))['x']
        # nodes = np.load('data/dynamic/new_square/l%d_m%.3f_final_trussCord.npz' %(layers, meshsize))['x']
        t_array = np.load('data/dynamic/new_square/laplace/l%d_m%.3f_fft_data_f%gHz.npz' %(layers, meshsize, freq))['t']
        
        # print(sum(sum(u2_time_array)))

        # u2L_idx = np.zeros([0])

        # for i in range(len(nodes)):
        #     x0 = nodes[i, 0]
        #     y0 = nodes[i, 1]
        #     if (x0 - L)**2 + y0**2 < w**2 :
        #         u2L_idx = np.append(u2L_idx, [2*i], axis=0) # need x_axis displacement only

        # u2L_idx = u2L_idx.astype('int32')    

        # u2L_array = np.average(u2_time_array[:, u2L_idx], axis=1)

        # fs = fft(u2L_array)

        # u2L_array = u2_time_array[0]

        # for i in range(int(len(u2_time_array)/2)-1):
        #     u2L_array += u2_time_array[2*(i+1)]

        # u2L_array = 2*np.average(u2_time_array, axis=0)
        u2L_array = (u2_time_array[0] + u2_time_array[2])/2

        x = t_array
        y = u2L_array

        N = len(t_array)
        T = t_array[1] - t_array[0]

        yf = fft(y)
        xf = fftfreq(N, T)[:N//2]

        # plt.plot(xf[:], 2.0/N* np.abs(yf[:N//2]))
        plt.plot(xf[:], 2.0/N*np.abs(yf[:N//2]))

    # plt.legend(("25 Hz", "50 Hz", "75 Hz", "100 Hz", "125 Hz", "150 Hz", "175 Hz", "200 Hz"), loc='best')
    # plt.title("Fourier")
    # plt.legend(["Layers 1", "Layers 2", "Layers 3"])

    plt.grid()
    plt.show()