import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manim

L = 100
t_step = L*L
N = 200
J = 2
H = -0.6
T = 1
m = int(L/2)
ising = np.random.choice([1,-1], size = (L,L))

# plt.imshow(ising, cmap='gray', interpolation='none')
# plt.show()

def self_eng_per(spin, i, j): #Periodic boundaries
    s = spin[i][j]
    dE = 2*J*s*spin[i-1][j]
    dE += 2*J*s*spin[i][j-1]
    dE += 2*J*s*spin[(i+1)%L][j]
    dE += 2*J*s*spin[i][(j+1)%L]
    dE += 2*H*s
    return dE

def self_eng_hel(spin, i, j): #Helical boundaries
    s = np.zeros(N)
    s[i + L*j] = spin[i][j]
    dE = 2*J*s*spin[i-1][j]
    dE += 2*J*s*spin[i][j-1]
    dE += 2*J*s*spin[(i+1)%L][j]
    dE += 2*J*s*spin[i][(j+1)%L]
    dE += 2*H*s
    return dE

def metro():
    for _ in range(t_step):
        i = np.random.randint(0,L)
        j = np.random.randint(0,L)
        r = np.random.random()
        dE = self_eng_per(ising, i, j)
        if dE <= 0:
            ising[i][j] = -ising[i][j]
        else:
            if r < np.exp(-dE/T):
                ising[i][j] = -ising[i][j]
            else:
                ising[i][j] = ising[i][j]

def update_1():
    for _ in range(0,N):
        metro()

update_1()

# plt.imshow(ising, cmap='gray', interpolation='none')
# plt.show()

#Switch magnetic field:
H = -H

R = 15

for i in range(m-2*R, m+2*R):
    for j in range(m-2*R, m+2*R):
        if ((i-m)**2 + (j-m)**2) < R**2:
            ising[i][j] *= -1

'''
mag = np.empty(N)

for i in range(0,N):
    mag[i] = np.sum(ising)/(t_step)
    metro()

plt.plot(mag)
plt.show()
'''

fig, ax = plt.subplots()
img = ax.imshow(ising, cmap='gray', interpolation='none', vmin=-1, vmax=1)


def update(frame):
    metro()
    img.set_data(ising)
    return [img]

ani = manim.FuncAnimation(fig, update, frames=50)
plt.tight_layout()
plt.show()


