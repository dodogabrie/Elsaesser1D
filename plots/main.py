import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import chart_studio.tools as tls
import datapane as dp 
import os

def matplotlib_animation():
    fig, axs = plt.subplots(1, 2, figsize=(15, 10))
    data_dir = 'data/'
    files = os.listdir(data_dir)
    x, zp, zm, u_init, b_init = np.loadtxt(data_dir + 'output00000.dat', unpack=True)
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
        _, zp, zm, u, b = np.loadtxt(data_dir + 'output' + ID + '.dat', unpack=True)
        fix_size()
        line[1].set_data(x, u)
        line[3].set_data(x, b)
        return line
    
    
    ani = animation.FuncAnimation(
        fig, animate, frames=num_data, interval=100, repeat=False)
    
    plt.show()

def ID_files(f): # Number of files
    return int(f.split('.')[0].split('output')[1])

def sort_files(ff): # Sort list of files
    ff.sort(key=lambda f: ID_files(f))
    return ff
def get_data(filename, Elsasser, data_dir):
    x, zp, zm, u, b = np.loadtxt(data_dir + filename, unpack=True)
    if Elsasser:
        return x, zp, zm
    else: 
        return x, u, b


def plotly_animation(get_embed = False, Elsasser = True):

    _, _, dt, T_steps, _ , _ = np.loadtxt('input.dat', unpack=True, comments = '!')
    def scatter_from_file(filename, Elsasser, data_dir):
        x, y, z = get_data(filename, Elsasser, data_dir)
        if Elsasser: 
            return [go.Scatter(x=x, y=y, name='zm'), go.Scatter(x=x, y=z, name='zp')]
        else:
            return [go.Scatter(x=x, y=y, name='u'), go.Scatter(x=x, y=z, name='b')]

    # Extract list of data sorted by time
    data_dir = 'data/nlin/'
    files = os.listdir(data_dir)
    files = sort_files(files)
    num_files = ID_files(files[-1]) + 1
    real_step = round(T_steps/num_files)

    tt = np.arange(0, T_steps*dt, real_step*dt)

    x, y, z = get_data(files[0], Elsasser, data_dir)
##    x, zp, zm, u, b = np.loadtxt(data_dir + files[0], unpack=True)
    fig = make_subplots(rows=1, cols=2, subplot_titles=('Z+', 'Z-'),
                        horizontal_spacing=0.051)
    fig.add_trace(go.Scatter(x=x, y=y, line_width=1, marker_size = 0), row=1, col=1)
    fig.add_trace(go.Scatter(x=x, y=z, line_width=1, marker_size = 0), row=1, col=2)

    frames = [dict(name = t,
              data = scatter_from_file(f, Elsasser, data_dir),
              traces=[0, 1] # the elements of the list [0,1,2] give info on the traces in fig.data
                                      # that are updated by the above three go.Scatter instances
              ) for f, t in zip(files, tt)]

    ymin = np.min(y) 
    ymax = np.max(y) 
    yA = ymax-ymin
    ymin -= 10/100 * yA
    ymax += 10/100 * yA
    zmin = np.min(z) 
    zmax = np.max(z) 
    zA = zmax-zmin
    zmin -= 10/100 * zA
    zmax += 10/100 * zA

    fig.layout.yaxis1.range = [ymin, ymax]
    fig.layout.yaxis2.range = [zmin, zmax]
    fig.layout.xaxis1.range = [np.min(x), np.max(x)]
    fig.layout.xaxis2.range = [np.min(x), np.max(x)]
    figa = go.Figure(data=fig.data, frames=frames, layout=fig.layout)
    # add slider
    sliders=[{"active": 0,
              "currentvalue": {"prefix": "Tempo Fisico="},
              "len": 0.9,
              "steps": [{"args": [[fr.name],{
                                    "frame": {"duration": 0, "redraw": True},
                                    "mode": "immediate",
                                    "fromcurrent": True,
                                    },
                                  ],
                        "label": fr.name,
                        "method": "animate",
                        }
                    for fr in figa.frames
                ]}]


    figa.update_layout(sliders=sliders)
    if get_embed:
        report = dp.Report(dp.Plot(figa))
        report.upload(name="copy Simulazione Z+-", open=True)
    else:
        figa.show()
    return

def plot_instant(Elsasser = True):
    _, _, dt, T_steps, _ , _ = np.loadtxt('input.dat', unpack=True, comments = '!')
    data_dir_nlin = 'data/nlin/'
    data_dir_lin = 'data/line/'
    files = os.listdir(data_dir_nlin)
    files = sort_files(files)
    num_files = ID_files(files[-1]) + 1
    real_step = round(T_steps/num_files)
    tt = np.arange(0, T_steps*dt, real_step*dt)
    I_list = [0., 2.2, 3.5, 5.75, 7.8]

    plt.rc('font', **{'size'   : 19})
    fig, axs = plt.subplots(2, len(I_list), figsize=(18, 8), sharey= 'row', sharex=True)
    for i, I in enumerate(I_list):
        idx = np.argmin(np.abs(tt - I))
        x, zp, zm = get_data(files[idx], Elsasser, data_dir_nlin)
        x, zp_l, zm_l = get_data(files[idx], Elsasser, data_dir_lin)
        axs[0][i].plot(x, zp, label='z+', c = 'Tab:blue')
        axs[0][i].plot(x, zp_l, lw = 0.9, ls = '--', alpha = 0.6, c = 'k')
        axs[1][i].plot(x, zm, label='z-', c = 'Tab:orange')
        axs[1][i].plot(x, zm_l, lw = 0.9, ls = '--', alpha = 0.6, c = 'k')
        axs[0][i].set_title('t = ' + str(round(tt[idx], 2)))
        axs[1][i].set_xlabel('x')
        if i == 0:
            axs[0][i].set_ylabel(r'$Z^{+}$')
            axs[1][i].set_ylabel(r'$Z^{-}$')
        beautify(axs[0][i])
        beautify(axs[1][i])
    plt.suptitle(r"Evoluzione della soluzione nelle variabili $Z^{\pm}$")
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.04)
    plt.savefig("figure/Zpm_evo_full.png", dpi = 200)
    plt.show()
    return

def beautify(ax):
    ax.grid(alpha = 0.3)
    ax.minorticks_on()
    ax.tick_params('x', which='major', direction='in', length=5)
    ax.tick_params('y', which='major', direction='in', length=5)
    ax.tick_params('x', which='minor', direction='in', length=3, bottom=True)


if __name__ == '__main__':
#    matplotlib_animation()
    plotly_animation(get_embed=False, Elsasser=True)
#    plot_instant()
