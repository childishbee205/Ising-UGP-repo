import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manim
from numba import njit, prange

L = 100
pi = 3.14159

t_step = L*L
steps = 1500
n_snaps = 50
snap = int(steps/n_snaps)
copies = 5

J = 1
J_1 = 1 #J_1 = T/J
J_2 = -0.7 #J_2 = H/J
X = 10 #Initial bubble radius

T = J_1*J
H = J_2*J

m = int(L/2)

ising = np.ones((L,L)) #Initialised stable 
#For random initialisation: np.random.choice([1,-1], size = (L,L))

@njit
def self_eng_per(spin, i, j, J_1, J_2, L): #Periodic boundaries
    s = spin[i][j]
    dE = 2*s*spin[i-1][j]/J_1
    dE += 2*s*spin[i][j-1]/J_1
    dE += 2*s*spin[(i+1)%L][j]/J_1
    dE += 2*s*spin[i][(j+1)%L]/J_1
    dE += 2*J_2*s/J_1
    return dE

@njit
def metro(ising, L, t_step, J_1, J_2, mode = 'g'): #Monte carlo update, L*L sweeps
    for _ in range(t_step):
        i = np.random.randint(0,L)
        j = np.random.randint(0,L)
        r = np.random.random()
        dE = self_eng_per(ising, i, j, J_1, J_2, L)
        if mode == 'g':
            if r < 1/(1+np.exp(dE)):
                ising[i][j] = -ising[i][j]
            else:
                ising[i][j] = ising[i][j]
        if mode == 'k':
            if dE <= 0 or r < np.exp(-dE):
                temp = ising[i][j]
                d = np.random.randint(0,2)
                r_2 = np.random.random()
                if d == 0:
                    z = int(np.sign(r_2-0.5))
                    ising[i][j] = ising[i][(j+z)%L]
                    ising[i][(j+z)%L] = temp
                if d == 1:
                    z = int(np.sign(r_2-0.5))
                    ising[i][j] = ising[(i+z)%L][j]
                    ising[(i+z)%L][j] = temp

@njit
def bubble_gen(y, ising, m): #Creates bubble of radius y
    R = int(y)+1
    for i in range(m-2*R, m+2*R):
        for j in range(m-2*R, m+2*R):
            if ((i-m)**2 + (j-m)**2) < y**2:
                ising[i][j] *= -1


@njit
def radius_t(mag,L): #Finds "radius" given magnetisation, assuming all flipped spins are in a disk
    R_t = float((1/np.sqrt(2*pi))*L*np.sqrt(1-mag))
    return R_t

bubble_gen(X, ising, m)
mag = float(np.sum(ising)/(t_step))
# plt.imshow(ising, cmap = 'grey', vmax=1, vmin=-1)
# plt.show()
print(mag)

avg_ising = np.zeros((n_snaps,L,L))
avg_R = np.zeros(n_snaps)

@njit 
def sweep(copies,X,m,steps,t_step,J_1,J_2,snap,L): #Returns averag
    R = np.zeros(n_snaps)
    memory = np.zeros((n_snaps,L,L))
    for k in range(0, copies):
        ising = np.ones((L,L))
        bubble_gen(X, ising, m)
        for i in range (0, steps):
            metro(ising, L, t_step, J_1, J_2)
            mag = float(np.sum(ising)/(t_step))
            if i%snap == 0:
                memory[i//snap] += np.copy(ising)
                R[i//snap] += radius_t(mag,L)

    return(memory/copies, R/copies)


#To save data
avg_ising, avg_R = sweep(copies,X,m,steps,t_step,J_1,J_2,snap,L)

np.savez_compressed('Nucleation_data_13.npz',ising=avg_ising,R_t=avg_R,J_1=J_1,J_2=J_2,X=X,snap=snap)

#Data visualisation
'''
print(J, H, T, np.sum(avg_ising[n_snaps-1])/(L*L))
plt.plot(avg_R-X)
plt.xlabel("Time")
plt.ylabel("Radius of bubble")
plt.title(f"L = {L}, T/J = {J_1}, H/J = {J_2}, T = {T}, R_0 = {X}")
plt.show()
plt.imshow(avg_ising[n_snaps//2], cmap = 'grey', vmax=1, vmin=-1)
plt.title(f"L = {L}, T/J = {J_1}, H/J = {J_2}")
plt.show()
# '''

#Obtains growth exponent of radius
'''
time = np.arange(len(avg_R))*snap
log_R_A = np.log(avg_R[3:10]-X)
log_time = np.log(time[3:10])
a, b = np.polyfit(log_time, log_R_A, 1)
print(a, b)
plt.plot(log_time, log_R_A)
plt.xlabel("log(Time)")
plt.ylabel("log(Radius)")
plt.title(f"L = {L}, J_1 = {J_1}, J_2 = {J_2}, T = {T}, R_0 = {X}")
plt.show()
'''
#Animation sector
# '''
ising = np.ones((L,L))
bubble_gen(X, ising, m)

fig, ax = plt.subplots()
ax.set_title(f'Nucleation in a magnetic field; T/J = {J_1}, H/J = {J_2}')
img1 = ax.imshow(ising, cmap='gray', origin = 'lower', interpolation='none', vmin=-1, vmax=1)

def update(frame):
    metro(ising, L, t_step, J_1, J_2)
    img1.set_data(ising)
    return img1

ani = manim.FuncAnimation(fig, update, frames=50)
# ani.save('H_nuc.gif')
plt.tight_layout()
plt.show()

# '''
