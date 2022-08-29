import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack

def plot_template(row = 1, col = 2, figsize = (12, 3)):
    fig, axs = plt.subplots(row, col, figsize = figsize)
    ax = axs[0]
    ax.set_xlabel('x', fontsize=15)
    ax.set_ylabel('u', fontsize=15)
    ax.grid(alpha = 0.3)
    ax.minorticks_on()
    ax.tick_params('x', which='major', direction='in', length=5)
    ax.tick_params('y', which='major', direction='in', length=5)
    ax.tick_params('y', which='minor', direction='in', length=3, left=True)
    ax.tick_params('x', which='minor', direction='in', length=3, bottom=True)
    ax = axs[1]
    ax.grid(alpha = 0.3)
    ax.minorticks_on()
    ax.tick_params('x', which='major', direction='in', length=5)
    ax.tick_params('y', which='major', direction='in', length=5)
    ax.tick_params('y', which='minor', direction='in', length=3, left=True)
    ax.tick_params('x', which='minor', direction='in', length=3, bottom=True)
    ax.set_xlabel('k', fontsize=15)
    ax.set_ylabel(r'$\left|u_k\right|^2$', fontsize=15)
    ax.set_yscale('log')
    ax.set_xscale('log')
    return fig, axs

def my_real_fft(u, dx):
    """
    Real Fast Fourier Transform.
    """
    N = len(u)
    # x : 2 * pi = y : 1 ---> unitary rate
    # => dy = dx/(2*pi)
    dy =  dx / (2 * np.pi)
    # fft(j) = (u * exp(-2*pi*i*j*np.arange(n)/n)).sum()
    fft = fftpack.fft(u) # Discret fourier transform 
    k = fftpack.fftfreq(N, dy) 
    return fft, k


def evaluate_energy_density_spectrum(u_t, dx, zero_energy = False):
    # final fft
    fft_u, k_u = my_real_fft(u_t, dx) # FFT of final u 
    mod_u_k = np.abs(fft_u) # |u_k|
    mod_fft2_u = (mod_u_k)**2 # |u_k|^2
    if zero_energy:
        mask_pos_k_u = k_u >= 0 # mask for positive k final
    else:
        mask_pos_k_u = k_u > 0 # mask for positive k final
    mod_fft2_u = mod_fft2_u[mask_pos_k_u]
    k_u = k_u[mask_pos_k_u]
    return mod_fft2_u, k_u

def plot_results(u, x, dx, label, axs = [None], ref = None, **kwargs):
    if (axs == None).all():
        fig, axs = plot_template()
    u_k2   , k = evaluate_energy_density_spectrum(u, dx)
    axs[0].plot(x, u, label = label, **kwargs)
    if ref != None:
        u_k2ref, k = evaluate_energy_density_spectrum(ref, dx)
        axs[1].scatter(k, u_k2/u_k2ref, s = 10, label = label)
    else:
        axs[1].scatter(k, u_k2, s = 10, label = label)
    return axs


if __name__ == '__main__':
    fig, axs = plot_template()
    x, y, _, _, _ = np.loadtxt('./data/output00000.dat', unpack=True)
    _, yf, _, _, _ = np.loadtxt('./data/output00001.dat', unpack=True)
    dx = x[1]-x[0]
    plot_results(y, x, dx, label = "init", axs = axs)
    plot_results(yf, x, dx, label = "final", axs = axs)
    plt.show()

