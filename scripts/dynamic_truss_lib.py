"""
Plane Truss Analysis (Dynamic)
Version 4
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix # ,lil_array
from scipy.sparse.linalg import spsolve
# from clamped_free_bar_end_force import bar_newmark_integration
from matplotlib.widgets import Slider

class PlaneTruss:
    def __init__(self, joints, members, E, A, rho):
        self.dim = 2

        self.joints = joints
        self.members = members

        self.n_joints = len(joints)
        self.n_members = len(members)

        if type(E) == np.ndarray:
            self.Es = np.array(E)
        else:
            self.Es = E*np.ones(self.n_members)

        if type(A) == np.ndarray:
            self.As = np.array(A)
        else:
            self.As = A*np.ones(self.n_members)

        self.Ls = self.compute_lengths()
        self.thetas = self.compute_angles()

        self.n_dofs = self.dim*self.n_joints
        self.dofs = np.array([np.nan for _ in range(self.n_dofs)])

        self.M = np.zeros((self.n_dofs, self.n_dofs))
        self.C = np.zeros((self.n_dofs, self.n_dofs))
        self.K = np.zeros((self.n_dofs, self.n_dofs))
        self.F = np.zeros(self.n_dofs)

        self.q0 = np.zeros(self.n_dofs)
        self.q0_dot = np.zeros(self.n_dofs)
        
        if type(rho) == np.ndarray:
            self.ms = np.array(rho)*self.Ls
        else: 
            self.ms = rho*self.Ls

    def length(self, idx):
        id1, id2 = self.members[idx]
        coords1 = self.joints[id1]
        coords2 = self.joints[id2]
        return np.sqrt(np.sum((coords2 - coords1)**2))

    def angle(self, idx):
        L = self.length(idx)
        id1, id2 = self.members[idx]
        # print(id1)
        # print(id2)
        coords1 = self.joints[id1]
        coords2 = self.joints[id2]
        dy = coords2[1] - coords1[1]
        dx = coords2[0] - coords1[0]
        if dx >=0: return np.arcsin(dy/L)
        else: return np.arcsin(-dy/L)

    def compute_lengths(self):
        Ls = np.zeros(self.n_members)
        for i in range(self.n_members):
            Ls[i] = self.length(i)
        return Ls

    def compute_angles(self):
        thetas = np.zeros(self.n_members)
        for i in range(self.n_members):
            thetas[i] = self.angle(i)
        return thetas

    def apply_constraints(self, constraints):
        for c in constraints:
            self.dofs[int(self.dim*c[0] + c[1])] = c[2]

    def apply_loads(self, loads):
        for l in loads:
            self.F[int(self.dim*l[0] + l[1])] = l[2]

    def member_stiffness(self, idx):
        K = np.zeros((4, 4))
        c = np.cos(self.thetas[idx])
        s = np.sin(self.thetas[idx])

        K[0, 0] = c*c
        K[0, 1] = c*s
        K[0, 2] = -c*c
        K[0, 3] = -c*s

        K[1, 0] = c*s
        K[1, 1] = s*s
        K[1, 2] = -c*s
        K[1, 3] = -s*s

        K[2, 0] = -c*c
        K[2, 1] = -c*s
        K[2, 2] = c*c
        K[2, 3] = c*s

        K[3, 0] = -c*s
        K[3, 1] = -s*s
        K[3, 2] = c*s
        K[3, 3] = s*s

        return self.Es[idx]*self.As[idx]*K/self.Ls[idx]

    def member_mass(self, idx):
        M = np.zeros((4, 4))

        M[0, 0] = 2
        M[0, 2] = 1

        M[1, 1] = 2
        M[1, 3] = 1

        M[2, 0] = 1
        M[2, 2] = 2

        M[3, 1] = 1
        M[3, 3] = 2

        return self.ms[idx]*self.Ls[idx]*M/6

    def compute_stiffness(self):
        for i in range(self.n_members):
            K_member = self.member_stiffness(i)

            id1, id2 = self.members[i]
            ids = []
            ids.append(self.dim*id1)
            ids.append(self.dim*id1 + 1)
            ids.append(self.dim*id2)
            ids.append(self.dim*id2 + 1)

            for j in range(2*self.dim):
                for k in range(2*self.dim):
                    self.K[ids[j], ids[k]] += K_member[j, k]

    def compute_mass(self):
        for i in range(self.n_members):
            M_member = self.member_mass(i)

            id1, id2 = self.members[i]
            ids = []
            ids.append(self.dim*id1)
            ids.append(self.dim*id1 + 1)
            ids.append(self.dim*id2)
            ids.append(self.dim*id2 + 1)

            for j in range(2*self.dim):
                for k in range(2*self.dim):
                    self.M[ids[j], ids[k]] += M_member[j, k]

    def enforce_constraints(self):
        n_support = 0
        K_support = []
        M_support = []
        C_support = []

        for i in range(self.n_dofs):
            if not np.isnan(self.dofs[i]):
                n_support += 1
                K_support.append(np.array(self.K[i]))
                M_support.append(np.array(self.M[i]))
                C_support.append(np.array(self.C[i]))
                for j in range(self.n_dofs):
                    if j == i:
                        self.K[i, j] = 1.0
                        self.C[i, j] = 1.0
                        self.M[i, j] = 1.0
                    else:
                        self.K[i, j] = 0.0
                        self.C[i, j] = 0.0
                        self.M[i, j] = 0.0

        self.n_support = n_support
        self.K_support = np.array(K_support)
        self.M_support = np.array(M_support)
        self.C_support = np.array(C_support)

    def compute_reactions(self):
        reactions = np.zeros(self.n_support)
        for i in range(self.n_support):
            reactions[i] = self.K_support[i] @ self.dofs
        self.reactions = reactions

    def compute_member_forces(self):
        member_forces = np.zeros(self.n_members)
        for i in range(self.n_members):
            id1, id2 = self.members[i]
            u1 = self.dofs[self.dim*id1 + 0]
            v1 = self.dofs[self.dim*id1 + 1]
            u2 = self.dofs[self.dim*id2 + 0]
            v2 = self.dofs[self.dim*id2 + 1]
            d1 = u1*np.cos(self.thetas[i]) + v1*np.sin(self.thetas[i])
            d2 = u2*np.cos(self.thetas[i]) + v2*np.sin(self.thetas[i])
            member_forces[i] = self.Es[i]*self.As[i]*(d2 - d1)/self.Ls[i]
        self.member_forces = member_forces

    def dynamic_constant_force(self, t):
        return self.F
    
    def dynamic_periodic_force(self, t):
        f = 345
        return self.F*np.cos(2*np.pi*f*t)


    def solve(self, h, n_iterations, q0=[0], q0_dot=[0], gamma=0.5, beta=0.25):
        p = self.dynamic_constant_force
        self.compute_stiffness()
        self.compute_mass()
        self.enforce_constraints()
        # self.C = Dalpha*self.K + Dbeta*self.M # find what Dalpha and Dbeta are supposed to be
        
        # M, C, K defined

        if len(q0) == 1: q0 = self.q0
        if len(q0_dot)==1: q0_dot = self.q0_dot

        # initial setup
        t = 0
    
        p0 = p(0)
        q0_ddot = np.linalg.inv(self.M)@(p0 - self.C@q0_dot - self.K@q0)

        qs = np.zeros([0, self.n_dofs])
        q_dots = np.zeros([0, self.n_dofs])
        q_ddots = np.zeros([0, self.n_dofs])

        qs = np.append(qs, [q0], axis=0)
        q_dots = np.append(q_dots, [q0_dot], axis=0)
        q_ddots = np.append(q_ddots, [q0_ddot], axis=0)

        for i in range(n_iterations):
            t = t + h
            q_dot_star = q_dots[i] + (1-gamma)*h*q_ddots[i]
            q_star = qs[i] + h*q_dots[i] + (0.5 - beta)*h*h*q_ddots[i]

            S = self.M + h*gamma*self.C + h*h*beta*self.K

            qn_ddot = np.linalg.inv(S)@(p(t) - self.C@q_dot_star - self.K@q_star)
            qn_dot = q_dot_star + h*gamma*qn_ddot
            qn = q_star + h*h*beta*qn_ddot

            qs = np.append(qs, [qn], axis=0)
            q_dots = np.append(q_dots, [qn_dot], axis=0)
            q_ddots = np.append(q_ddots, [qn_ddot], axis=0)

        self.dofs = qs
        return qs, q_dots, q_ddots

    def plot(self, deformed=True, mag=1):
        TOL = 1e-6
        self.fig, self.ax = plt.subplots()

        for i in range(self.n_members):
            id1, id2 = self.members[i]
            self.ax.plot(
                [self.joints[id1, 0], self.joints[id2, 0]],
                [self.joints[id1, 1], self.joints[id2, 1]],
                '-', color='gray', linewidth=5)

        self.ax.scatter(self.joints[:, 0], self.joints[:, 1], c='b', s=20)
        

        if deformed:
            for i in range(self.n_members):
                id1, id2 = self.members[i]
                x1 = self.joints[id1, 0] + mag*self.dofs[0][2*id1 + 0]
                y1 = self.joints[id1, 1] + mag*self.dofs[0][2*id1 + 1]
                x2 = self.joints[id2, 0] + mag*self.dofs[0][2*id2 + 0]
                y2 = self.joints[id2, 1] + mag*self.dofs[0][2*id2 + 1]

                color = 'black'
                # if self.member_forces[i] > TOL:
                #     color = 'blue'
                # elif self.member_forces[i] < -TOL:
                #     color = 'red'

                self.ax.plot(
                    [x1, x2], [y1, y2],
                    '-', color=color, linewidth=1)
                
        x_ends = np.array(self.ax.get_xlim())
        y_ends = np.array(self.ax.get_ylim())

        self.ax.set_xlim(2*x_ends)
        self.ax.set_ylim(35*y_ends)
        self.ax.set_aspect(True)

        # Make a horizontal slider to control the frequency.
        axfreq = self.fig.add_axes([.27, 0.25, 0.5, 0.05]) # left, bottom, width, height
        freq_slider = Slider(
            ax=axfreq,
            label='Time (s)',
            valmin=0,
            valmax=0.1,
            valinit=0,
            valstep=np.linspace(0, 0.1, 101),
        )

        freq_slider.on_changed(self.update)

        plt.show()

    def update(self, val, deformed=True, mag=1):
        TOL = 1e-6
        # fig, ax = plt.subplots()

        val = val/10*200
        val = int(val)

        x_ends = np.array(self.ax.get_xlim())
        y_ends = np.array(self.ax.get_ylim())

        self.ax.cla()

        self.ax.set_xlim(x_ends)
        self.ax.set_ylim(y_ends)

        for i in range(self.n_members):
            id1, id2 = self.members[i]
            self.ax.plot(
                [self.joints[id1, 0], self.joints[id2, 0]],
                [self.joints[id1, 1], self.joints[id2, 1]],
                '-', color='gray', linewidth=3)

        self.ax.scatter(self.joints[:, 0], self.joints[:, 1], c='b', s=20)
        

        if deformed:
            for i in range(self.n_members):
                id1, id2 = self.members[i]
                x1 = self.joints[id1, 0] + mag*self.dofs[val][2*id1 + 0]
                y1 = self.joints[id1, 1] + mag*self.dofs[val][2*id1 + 1]
                x2 = self.joints[id2, 0] + mag*self.dofs[val][2*id2 + 0]
                y2 = self.joints[id2, 1] + mag*self.dofs[val][2*id2 + 1]

                color = 'black'
                # if self.member_forces[i] > TOL:
                #     color = 'blue'
                # elif self.member_forces[i] < -TOL:
                #     color = 'red'

                self.ax.plot(
                    [x1, x2], [y1, y2],
                    '-', color=color, linewidth=1)

def p(t):
    p = np.zeros(12)
    p[-1] = 1e4
    return p

if __name__ == '__main__':
    L = 1.0
    rho = 2710

    joints = np.array([
        [0.0, 0.0],
        [0.0, L], 
        [L, 0.0],
        [L, L], 
        [2*L, 0.0], 
        [2*L, L]
    ])

    members = np.array([
        [0, 2], 
        [0, 3], 
        [1, 3], 
        [2, 3], 
        [2, 4],
        [3, 4], 
        [3, 5], 
        [4, 5],
    ], dtype=int)

    E = 2.1e11
    A = 1e-4

    truss = PlaneTruss(joints, members, E, A, rho)

    constraints = [
        [0, 0, 0.0],
        [0, 1, 0.0],
        [1, 0, 0.0],
        [1, 1, 0.0]
    ]
    truss.apply_constraints(constraints)

    P = 1e4
    Q = 0
    loads = [
        [4, 1, -P]
    ]
    truss.apply_loads(loads)

    truss.solve(h=0.05, n_iterations=200)

    # q0dot = np.zeros([2, truss.n_joints])
    # qs, qdots, qddots = bar_newmark_integration(joints, truss.M, truss.C, truss.K, q0dot, p, h=0.05, n_iterations=200) # need to redefine same thing within truss class definition.

    # print('Reactions: ')
    # print(truss.reactions)

    # print('Internal forces: ')
    # print(truss.member_forces)

    truss.plot()
