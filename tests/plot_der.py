import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    data = np.loadtxt('./data/output00001.dat')
    x = data[:, 0]
    f = data[:, 1]
    df = data[:, 2]
    plt.plot(x, f)
    plt.plot(x, df)
    plt.show()

