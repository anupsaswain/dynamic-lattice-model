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

    freqs_array = np.linspace(1, 50, 5000)
    amps_array = np.zeros_like(freqs_array)

    i = 0
    for freq in freqs_array:
        layers = 1
        meshsize = .1
        conn_node_tol = 0.5
        H_meshsize = meshsize #*(3**.5)/2
        w = H_meshsize*layers

# (50/freq)**1.5*
        u2_time_array = np.load('data/dynamic/new_square/laplace/trial1_f%g.npz' %(freq))['x']
        t_array = np.load('data/dynamic/new_square/laplace/trial1_f%g.npz' %(freq))['t']
        
        u2L_array = (u2_time_array[0] + u2_time_array[2])/2

        y = u2L_array

        N = len(t_array)
        yf = fft(y)

        # amps = (np.sort(2.0/N*np.abs(yf[:N//2]))/0.5)

        # amps_array[i] = np.average(amps[-1:])
        amps_array[i] = np.max(2.0/N*np.abs(yf[:N//2]))/0.5

        i +=1

    plt.plot(freqs_array, amps_array)
    plt.plot([0, 50], [1, 1], "--", 0.5)
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Amplitude of output/input")

    # plt.legend(("25 Hz", "50 Hz", "75 Hz", "100 Hz", "125 Hz", "150 Hz", "175 Hz", "200 Hz"), loc='best')
    # plt.title("Fourier")
    # plt.legend(["Layers 1", "Layers 2", "Layers 3"])

    plt.grid()
    plt.show()