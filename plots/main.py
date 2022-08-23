import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os


fig, axs = plt.subplots(1, 2, figsize=(15, 10))
data_dir = 'data/'
files = os.listdir(data_dir)
x, u_init, b_init = np.loadtxt(data_dir + 'output00000.dat', unpack=True)
num_data = len(files)
print('Found', num_data, 'data files')

minx = np.min(x)
maxx = np.max(x)
min1 = np.min(u_init)
max1 = np.max(u_init)
min_tot = min(min1, min1)
min2 = np.min(b_init)
max2 = np.max(b_init)
max_tot = max(max1, max2)
line11, = axs[0].plot(x, u_init, alpha=0.2)
line12, = axs[0].plot([], [])
line21, = axs[1].plot(x, b_init, alpha=0.2)
line22, = axs[1].plot([], [])
line = [line11, line12, line21, line22]


def fix_size():
    for ax in axs:
        ax.set_ylim(min_tot - 0.1 * np.abs(min_tot), max_tot + 0.1 * max_tot)
        ax.set_xlim(minx, maxx)


fix_size()


def animate(i):
    leni = len(str(i))
    ID = '0'*(5 - leni) + str(i)
    _, u, b = np.loadtxt(data_dir + 'output' + ID + '.dat', unpack=True)
    fix_size()
    line[1].set_data(x, u)
    line[3].set_data(x, b)
    return line


ani = animation.FuncAnimation(
    fig, animate, frames=num_data, interval=100, repeat=False)

plt.show()
