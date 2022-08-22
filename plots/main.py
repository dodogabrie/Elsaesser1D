import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


fig, axs = plt.subplots(1, 2, figsize=(15, 10))
data_dir = 'data/'
x, u_init, b_init = np.loadtxt(data_dir + 'output00000.dat', unpack=True)
num_data = 500

minx = np.min(x)
maxx = np.max(x)
min1 = np.min(u_init)
max1 = np.max(u_init)
min2 = np.min(b_init)
max2 = np.max(b_init)
line11, = axs[0].plot(x, u_init, alpha=0.2)
line12, = axs[0].plot([], [])
line21, = axs[1].plot(x, b_init, alpha=0.2)
line22, = axs[1].plot([], [])
line = [line11, line12, line21, line22]


def fix_size():
    axs[0].set_ylim(min1 - 0.5 * np.abs(min1), max1 + 0.5 * max1)
    axs[0].set_xlim(minx, maxx)
    axs[1].set_ylim(min2 - 0.5 * np.abs(min2), max2 + 0.5 * max2)
    axs[1].set_xlim(minx, maxx)


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
    fig, animate, frames=num_data, interval=40, repeat=False)

plt.show()
