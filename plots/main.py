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

def plotly_animation(get_embed = False, Elsasser = True):

    _, _, dt, T_steps, _ , _ = np.loadtxt('input.dat', unpack=True, comments = '!')
    def ID_files(f): # Number of files
        return int(f.split('.')[0].split('output')[1])

    def sort_files(ff): # Sort list of files
        ff.sort(key=lambda f: ID_files(f))
        return files

    def get_data(filename, Elsasser = True):
        x, zp, zm, u, b = np.loadtxt(data_dir + filename, unpack=True)
        if Elsasser:
            return x, zp, zm
        else: 
            return x, u, b

    def scatter_from_file(filename, Elsasser):
        x, y, z = get_data(filename, Elsasser)
        if Elsasser: 
            return [go.Scatter(x=x, y=y, name='zm'), go.Scatter(x=x, y=z, name='zp')]
        else:
            return [go.Scatter(x=x, y=y, name='u'), go.Scatter(x=x, y=z, name='b')]

    # Extract list of data sorted by time
    data_dir = 'data/'
    files = os.listdir(data_dir)
    files = sort_files(files)
    num_files = ID_files(files[-1]) + 1
    real_step = round(T_steps/num_files)

    tt = np.arange(0, T_steps*dt, real_step*dt)

    x, y, z = get_data(files[0], Elsasser = Elsasser)
##    x, zp, zm, u, b = np.loadtxt(data_dir + files[0], unpack=True)
    fig = make_subplots(rows=1, cols=2, subplot_titles=('Z+', 'Z-'),
                        horizontal_spacing=0.051)
    fig.add_trace(go.Scatter(x=x, y=y, line_width=1, marker_size = 0), row=1, col=1)
    fig.add_trace(go.Scatter(x=x, y=z, line_width=1, marker_size = 0), row=1, col=2)

    frames = [dict(name = t,
              data = scatter_from_file(f, Elsasser),
              traces=[0, 1] # the elements of the list [0,1,2] give info on the traces in fig.data
                                      # that are updated by the above three go.Scatter instances
              ) for f, t in zip(files, tt)]

    fig.layout.yaxis1.range = [-1.1, 1.1]
    fig.layout.yaxis2.range = [-1.1, 1.1]
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
        report.upload(name="figura di esempio", open=True)
    else:
        figa.show()
    return

if __name__ == '__main__':
#    matplotlib_animation()
    plotly_animation(get_embed=False, Elsasser=True)
