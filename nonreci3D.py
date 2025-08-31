import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manim
import math


L = 30
sweep = L*L*L
N = 600
J = 1
K = 2
T = 1

ising_A =  np.ones((L,L,L)) #np.random.choice([1,-1], size = (L,L,L)) #np.ones((L,L,L))
ising_B =  (-1)*np.ones((L,L,L)) #np.random.choice([1,-1], size = (L,L,L)) #(-1)*np.ones((L,L,L))
lattice = [ising_A, ising_B]


def self_eng_open(spin, i, j, k, lc): #Open boundaries
    dE = 0
    s = spin[lc][i][j][k]
    if i>0:
        dE += 2*J*s*spin[lc][i-1][j][k]
    if j>0:
        dE += 2*J*s*spin[lc][i][j-1][k]
    if i<L-1:
        dE += 2*J*s*spin[lc][i+1][j][k]
    if j<L-1:
        dE += 2*J*s*spin[lc][i][j+1][k]
    if k>0:
        dE += 2*J*s*spin[lc][i][j][k-1]
    if k<L-1:
        dE += 2*J*s*spin[lc][i][j][k+1]

    dE += (2-4*lc)*K*s*spin[lc-1][i][j][k]
    return dE

def self_eng_per(spin, i, j, k, lc): #Periodic boundaries
    dE = 0
    s = spin[lc][i][j][k]
    dE += 2*J*s*spin[lc][i-1][j][k]
    dE += 2*J*s*spin[lc][i][j-1][k]
    dE += 2*J*s*spin[lc][(i+1)%L][j][k]
    dE += 2*J*s*spin[lc][i][(j+1)%L][k]
    dE += (2-4*lc)*K*s*spin[lc-1][i][j][k]
    dE += 2*J*s*spin[lc][i][j][k-1]
    dE += 2*J*s*spin[lc][i][j][(k+1)%L]

    return dE
    
m_A = np.zeros(N)
m_B = np.zeros(N)

m_A[0] = np.array(ising_A).sum()/(L*L*L)
m_B[0] = np.array(ising_B).sum()/(L*L*L)
avgm_A = m_A[0]
avgm_B = m_B[0]
print(avgm_A, avgm_B)


def metro():
    for _ in range(0, sweep):
        lc = np.random.randint(0,2)
        i = np.random.randint(0,L)
        j = np.random.randint(0,L)
        k = np.random.randint(0,L)
        r = np.random.random()
        w = 1/(1+np.exp((self_eng_per(lattice, i, j, k, lc))/T))
        if r <= w:
            lattice[lc][i][j][k] = (-1)*lattice[lc][i][j][k]
        else:
            pass
    



#Magnetisation section

# '''
l = np.zeros(N-1)

for s in range(1,N):
    metro()
    m_A[s] = np.array(ising_A).sum()/(L*L*L)
    m_B[s] = np.array(ising_B).sum()/(L*L*L)
    l[s-1] = (m_A[s-1]*m_B[s] - m_B[s-1]*m_A[s])
    
avgm_A = np.sum(m_A)/N
avgm_B = np.sum(m_B)/N
R = math.sqrt((avgm_A)**2 + (avgm_B)**2)
avg_L = np.sum(l)/(N-1)


print(avgm_A, avgm_B, R, avg_L)
plt.plot(m_A, color = 'b')
plt.plot(m_B, color = 'r')
plt.plot(l, color = 'g')
plt.title('L = 60, J = 1, K = 3')
plt.show()
# '''

#Plane animation
'''
mid = int(L/2)
fig, ax = plt.subplots(1,2)
img1 = ax[0].imshow(lattice[0][:][:][mid], cmap='gray', interpolation='none', vmin=-1, vmax=1)
img2 = ax[1].imshow(lattice[1][:][:][mid], cmap='gray', interpolation='none', vmin=-1, vmax=1)


for axe in ax:
    axe.axis('off') 


def update(frame):
    metro()

    img1.set_data(lattice[0][:][:][mid])
    img2.set_data(lattice[1][:][:][mid])
    return [img1, img2]

ani = manim.FuncAnimation(fig, update, frames=N, blit=True)
plt.tight_layout()
plt.show()

'''