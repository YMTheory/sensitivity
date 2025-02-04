import matplotlib.pyplot as plt
import numpy as np

def get_contour(delta_chisquare, level):
    dm_square_cents         = np.logspace(-2, 1, 100)
    sin2theta_square_cents  = np.logspace(-2, 0, 100)

    fig, ax = plt.subplots()
    X, Y = np.meshgrid(sin2theta_square_cents, dm_square_cents)
    CS = ax.contour(X, Y, delta_chisquare, levels=[level], cmap='Spectral', linewidths=3, linestyles='-')
    plt.close()
    
    p = CS.get_paths()[0]
    v = p.vertices
    x, y = v[:,0], v[:,1]

    return x, y

