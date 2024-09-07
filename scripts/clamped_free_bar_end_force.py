# from dynamic_solve import newmark_integration
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

def bar_newmark_integration(nodes, M, C, K, q0_dot, p, h, n_iterations, gamma=0.5, beta=0.25):
    t = 0
    q0 = nodes - nodes
    p0 = p(0) # p is force, function of time
    q0_doubledot = np.linalg.inv(M)@(p0 - C@q0_dot - K@q0)
    dim = int(len(nodes))
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
        qn = np.round(qn, 2)
        # print(qn)
        qs = np.append(qs, [qn], axis=0)
        q_dots = np.append(q_dots, [qn_dot], axis=0)
        q_doubledots = np.append(q_doubledots, [qn_doubledot], axis=0)

    qs = qs + np.array([nodes for i in range(n_iterations+1)])

    return qs, q_dots, q_doubledots

# An example: bar with multiple elements

# initial parameters

L = 10 # m
N = 10 # number of elements
le = L/N # length of each element

E = 7e10 # Aluminium standard
A = 0.25 # m^2

rho = 2710 # kg/m^3, Aluminium

m = rho*A # mass per unit length (kg/m)

# frequency of force applied
w = np.pi/6

# stiffness matrix
K = np.zeros([N+1, N+1])

K[0, 0:2] = E*A/le*np.array([2, -1])

K[1, 0:3] = E*A/le*np.array([-1, 2, -1])

for i in range(N-2):
    K[i+2, i+1:(i+4)] = E*A/le*np.array([-1, 2, -1])

K[N, (N-1):(N+1)] = E*A/le*np.array([-1, 1])

# mass matrix
M = np.zeros([N+1, N+1])

M[0, 0:2] = m*le/6*np.array([4, 1])

for i in range(N-1):
    M[i+1, i:(i+3)] = m*le/6*np.array([1, 4, 1])

M[N, (N-1):(N+1)] = m*le/6*np.array([1, 2])

# damping matrix 
# C = Dalpha*K + Dbeta*M
C = np.zeros([N+1, N+1])

# external force on end only. 
# force as a function of time

# print(M)
# print(K)

def force(time):
    p = np.zeros(N+1)
    p[N] = 1e9 #*np.sin(w*time)
    return p

# 10 elements => 11 nodes

nodes = np.linspace(0, L, N+1)

q0_dot = np.zeros([N+1])

qs, q_dots, q_doubledots = bar_newmark_integration(nodes, M, C, K, q0_dot, force, h=0.05, n_iterations=200)

# print("Displacements:")
# qs = np.round(qs, 2)
# print(qs)
# print("Velocities:")
# print(q_dots)
# print("Accelerations:")
# print(q_doubledots)

# plotting displacement response

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

def update(val):
    i = val/0.05 # h value set
    i = int(i) 
    line.set_xdata(qs[i])
    fig.canvas.draw_idle()

freq_slider.on_changed(update)

plt.show()