import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manim
import math


L = 100
m = int(L/2)
sweep = L*L
N = 700
J = 2
K = 1
T = 1
H = 0

ising_A =  np.ones((L,L)) #np.random.choice([1,-1], size = (L,L)) 
ising_B =  (-1)*np.ones((L,L)) #np.random.choice([1,-1], size = (L,L)) 
lattice = [ising_A, ising_B]

'''
for i in range (0, L):
    for j in range (m, L):
        ising_A[i][j] *= -1
        ising_B[i][j] *= -1
'''


def self_eng_open(spin, i, j, lc): #Open boundaries
    dE = 0
    s = spin[lc][i][j]
    if i>0:
        dE += 2*J*s*spin[lc][i-1][j]
    if j>0:
        dE += 2*J*s*spin[lc][i][j-1]
    if i<L-1:
        dE += 2*J*s*spin[lc][i+1][j]
    if j<L-1:
        dE += 2*J*s*spin[lc][i][j+1]
    dE += (2-4*lc)*K*s*spin[lc-1][i][j]
    return dE

def self_eng_per(spin, i, j, lc): #Periodic boundaries
    dE = 0
    s = spin[lc][i][j]
    dE += 2*J*s*spin[lc][i-1][j]
    dE += 2*J*s*spin[lc][i][j-1]
    dE += 2*J*s*spin[lc][(i+1)%L][j]
    dE += 2*J*s*spin[lc][i][(j+1)%L]
    dE += ((2-4*lc)*K*s*spin[lc-1][i][j] + 2*H*(1-lc)*s)
    return dE
    
m_A = np.zeros(N)
m_B = np.zeros(N)

m_A[0] = np.array(ising_A).sum()/(L*L)
m_B[0] = np.array(ising_B).sum()/(L*L)
avgm_A = m_A[0]
avgm_B = m_B[0]
print(avgm_A, avgm_B)


def metro():
    for x in range(0, sweep):
        lc = np.random.randint(0,2)
        i = np.random.randint(0,L)
        j = np.random.randint(0,L)
        r = np.random.random()
        w = 1/(1+np.exp((self_eng_per(lattice, i, j, lc))/T))
        if r <= w:
            lattice[lc][i][j] = (-1)*lattice[lc][i][j]
        else:
            pass
    

R = 20

for i in range(m-2*R, m+2*R):
    for j in range(m-2*R, m+2*R):
        if ((i-m)**2 + (j-m)**2) < R**2:
            lattice[0][i][j] *= -1
            # lattice[1][i][j] *= -1

'''
fig, ax = plt.subplots(1,2)
img1 = ax[0].imshow(lattice[0], cmap='gray', interpolation='none', vmin=-1, vmax=1)
img2 = ax[1].imshow(lattice[1], cmap='gray', interpolation='none', vmin=-1, vmax=1)


for axe in ax:
    axe.axis('off') 


def update(frame):
    metro()

    img1.set_data(lattice[0])
    img2.set_data(lattice[1])
    return [img1, img2]

ani = manim.FuncAnimation(fig, update, frames=N, blit=True)
plt.tight_layout()
plt.show()

'''
l = np.zeros(N-1)

# '''
for s in range(1,N):
    metro()
    m_A[s] = np.array(ising_A).sum()/(L*L)
    m_B[s] = np.array(ising_B).sum()/(L*L)
    l[s-1] = 10*(m_A[s-1]*m_B[s] - m_B[s-1]*m_A[s])
    
avgm_A = np.sum(m_A)/N
avgm_B = np.sum(m_B)/N
R = math.sqrt((avgm_A)**2 + (avgm_B)**2)
avg_L = np.sum(l)/(N-1)


print(avgm_A, avgm_B, R, avg_L)
plt.plot(m_A, color = 'b')
plt.plot(m_B, color = 'r')
# plt.plot(l, color = 'g')
plt.title('L = 50, J = 0.75, K = 1')
plt.show()
# '''









        



    