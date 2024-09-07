
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_array, lil_matrix
from scipy.sparse.linalg import spsolve

class PlaneTruss:
    def __init__(self, joints, members, E, A):
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

        self.K = np.zeros((self.n_dofs, self.n_dofs))
        self.F = np.zeros(self.n_dofs)

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

    def enforce_constraints(self):
        n_support = 0
        K_support = []

        for i in range(self.n_dofs):
            if not np.isnan(self.dofs[i]):
                n_support += 1
                K_support.append(np.array(self.K[i]))
                for j in range(self.n_dofs):
                    if j == i:
                        self.K[i, j] = 1.0
                    else:
                        self.K[i, j] = 0.0
                        self.K[j, i] = 0.0
 
        self.n_support = n_support
        self.K_support = np.array(K_support)

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

    def solve(self):
        self.compute_stiffness()
        # print(self.K)
        self.enforce_constraints()
        # print(self.K)

        # breakpoint()
        A = self.K 
        A = lil_matrix(A)
        A = A.tocsc()
        # A = A.tolist()
        # A = A.tocsr()
        self.dofs = spsolve(A, self.F)

        self.compute_reactions()
        self.compute_member_forces()


    def plot(self, deformed=True, mag=25):
        TOL = 1e-6
        fig = plt.figure(figsize=(6, 6))

        for i in range(self.n_members):
            id1, id2 = self.members[i]
            plt.plot(
                [self.joints[id1, 0], self.joints[id2, 0]],
                [self.joints[id1, 1], self.joints[id2, 1]],
                '-', color='gray', linewidth=5)

        plt.scatter(self.joints[:, 0], self.joints[:, 1], c='b', s=20)
        

        if deformed:
            for i in range(self.n_members):
                id1, id2 = self.members[i]
                x1 = self.joints[id1, 0] + mag*self.dofs[2*id1 + 0]
                y1 = self.joints[id1, 1] + mag*self.dofs[2*id1 + 1]
                x2 = self.joints[id2, 0] + mag*self.dofs[2*id2 + 0]
                y2 = self.joints[id2, 1] + mag*self.dofs[2*id2 + 1]

                color = 'black'
                if self.member_forces[i] > TOL:
                    color = 'blue'
                elif self.member_forces[i] < -TOL:
                    color = 'red'

                plt.plot(
                    [x1, x2], [y1, y2],
                    '-', color=color, linewidth=3)

        plt.show()

if __name__ == '__main__':
    L = 1.0

    joints = np.array([
        [0.0, 0.0],
        [0.0, L], 
        [L, 0.0],
        [L, L], 
        [2*L, 0.0], 
        [2*L, L]
    ])

    members = np.array([
        [0, 1], 
        [0, 2], 
        [1, 3], 
        [2, 3], 
        [2, 4], 
        [3, 5], 
        [4, 5],
    ], dtype=int)

    E = 2.1e11
    A = 1e-4

    truss = PlaneTruss(joints, members, E, A)

    constraints = [
        [0, 0, 0.0],
        [0, 1, 0.0],
        [1, 0, 0.0],
        [1, 1, 0.0]
    ]
    truss.apply_constraints(constraints)

    P = 1e-10
    Q = 0
    loads = [
        [4, 1, -P]
    ]
    truss.apply_loads(loads)

    truss.solve()

    print('Reactions: ')
    print(truss.reactions)

    print('Internal forces: ')
    print(truss.member_forces)

    truss.plot()
