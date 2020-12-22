# plot_functions.py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def new_cmap():
    ''' Creates custom colormap for HEPD telescope and MEPD detector plot. '''
    jet = plt.cm.get_cmap('jet', 256)
    newcolors = jet(np.linspace(0, 1, 256))
    white = np.array([256/256, 256/256, 256/256, 1])
    newcolors[0, :] = white
    new_cmap = matplotlib.colors.ListedColormap(newcolors)
    return new_cmap

