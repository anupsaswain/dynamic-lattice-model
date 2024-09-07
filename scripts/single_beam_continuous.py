# Author: Anupsa Swain
# Date: 29 May 2024

import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

def u(p0, L, E, A, rho, x, t, k=1000):
    
    m = rho*A
    c = np.sqrt(E*A/m)
    elems = 8*p0*L/((np.pi**2)*E*A)*np.array([((-1)**(s-1))/((2*s-1)**2)
                      *np.sin((2*s-1)*np.pi*x/2/L)
                      *(1-np.cos((2*s-1)*np.pi*c*t/2/L)) for s in range(k)])

    return np.sum(elems)

def u_array(p0, L, E, A, rho, x_array, t, k=1000):
    us = np.zeros_like(x_array)

    for i in range(len(x_array)):
        us[i] = u(p0, L, E, A, rho, x_array[i], t, k)

    return us

def update(val):
    i = val/0.05 # h value set
    i = int(i) 
    line.set_xdata(qs[i])
    fig.canvas.draw_idle()

if __name__ == "__main__":

    n = 10 # number of elements
    L = 10 
    E = 7e10 # Aluminium standard
    A = 0.25 # m^2
    rho = 2710 # kg/m^3, Aluminium

    xs0 = np.linspace(0, L, n)

    P = 1e9

    us0 = u_array(P, L, E, A, rho, xs0, t=0)

    # print(us0)

    # '''

    h = 0.05
    n_iterations = 200

    T = h*n_iterations # total time

    qs = np.zeros([n_iterations+1, n])

    for t_count in range(n_iterations):
        t = t_count*h
        qs[t_count] = xs0 + u_array(P, L, E, A, rho, xs0, t)

    # Create the figure and the line that we will manipulate
    fig, ax = plt.subplots()
    line, = ax.plot(qs[0], np.ones_like(qs[0]), lw=2, marker='o')
    # annotations = 
    # for i in range(N+1):
    #     ax.annotate(i, xy=[qs[0, i], 1],c='magenta')
    ax.set_xlabel("Position (m)")
    ax.set_xlim([0, 3*L/2])
    ax.set_ylim([0, 2])
    ax.set_aspect(True)
    ax.set_yticks([])

    # Make a horizontal slider to control the frequency.
    axfreq = fig.add_axes([.27, 0.2, 0.5, 0.05])
    freq_slider = Slider(
        ax=axfreq,
        label='Time (s)',
        valmin=0,
        valmax=10,
        valinit=0,
        valstep=np.linspace(0, 10, 201),
    )

    freq_slider.on_changed(update)

    plt.show()

    # '''


