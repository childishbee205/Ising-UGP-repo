import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

L = 700
N = 500
J = 1
K = -3
T = 1
ising =  np.ones(L) #np.random.choice(([-1,1]), size = L)

def energy(i, ising):
    dE = 0
    dE = ising[i]*(2*J*ising[(i+1)%L] + 2*K*ising[i-1])
    return dE

def montecarlo():
    for _ in range(0, L):
        i = np.random.randint(0, L)
        dE = energy(i, ising)
        p = 1/(1+np.exp(-dE/T))
        r = np.random.random()
        if r<p:
            ising[i] = -ising[i]
        else:
            pass

steps = N
history = np.zeros((steps, L), dtype=int)

for t in range(steps):
    montecarlo()
    history[t] = ising.copy()

# Prepare scatter plot data
x_coords, y_coords, colors = [], [], []
for t in range(steps):
    for i in range(L):
        x_coords.append(i)
        y_coords.append(t)
        colors.append('black' if history[t, i] == 1 else 'white')

m = sum(ising)
print(m/700)

plt.figure(figsize=(10, 6))
plt.scatter(x_coords, y_coords, c=colors, s=10, edgecolors='k', linewidths=0.2)
plt.xlabel("Spin index")
plt.ylabel("Time step")
plt.gca().invert_yaxis()  # time flows downward
plt.show()


