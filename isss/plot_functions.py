# plot_functions.py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from data_functions import *

def new_cmap():
    ''' Creates custom colormap for HEPD telescope and MEPD detector plot. '''
    jet = plt.cm.get_cmap('jet', 256)
    newcolors = jet(np.linspace(0, 1, 256))
    white = np.array([256/256, 256/256, 256/256, 1])
    newcolors[0, :] = white
    new_cmap = matplotlib.colors.ListedColormap(newcolors)
    return new_cmap

def plot_pc1(fig, outer, HEPD, MEPD):
    inner = gridspec.GridSpecFromSubplotSpec(1, 2,
                    subplot_spec=outer[0], wspace=0.1, hspace=0.1)
    # Plot HEPD PC1 data.
    ax = plt.Subplot(fig, inner[0])
    ax.plot(HEPD.pc1, HEPD.time - HEPD.time[0], '-k')
    fig.add_subplot(ax)
    ax.set_title('HEPD: Time vs PC1')
    ax.set_xlabel('PC1')
    ax.set_ylabel('Time (sec)')
    # Divide PC1 data into MEPD-A and MEPD-B
    MEPD_pc1_A, MEPD_pc1_B = sliceAB(MEPD.pc1, MEPD.dt, 3)
    MEPD_time_A, MEPD_time_B = sliceAB(MEPD.time, MEPD.dt, 3)
    # Plot MEPD-A and MEPD-B PC1 data.
    ax = plt.Subplot(fig, inner[1])
    ax.plot(MEPD_pc1_A, MEPD_time_A - MEPD_time_A[0], '-k', label='MEPD-A')
    ax.plot(MEPD_pc1_B, MEPD_time_B - MEPD_time_B[0], '-r', label='MEPD-B')
    fig.add_subplot(ax)
    ax.set_title('MEPD: Time vs PC1')
    ax.set_xlabel('PC1')
    plt.legend()

def plot_mag(fig, outer, HEPD, MEPD):
    inner = gridspec.GridSpecFromSubplotSpec(1, 1,
                    subplot_spec=outer[2], wspace=0.1, hspace=0.1)
    mag_avg = B_avg(HEPD.mag)
    # Plot magnetic field data.
    ax = plt.Subplot(fig, inner[0])
    ax.plot(HEPD.time, HEPD.mag[:, 0], 'k', label='Bx')
    ax.plot(HEPD.time, HEPD.mag[:, 1], 'b', label='By')
    ax.plot(HEPD.time, HEPD.mag[:, 2], 'r', label='Bz')
    ax.plot(HEPD.time, HEPD.mag[:, 4], '--k', label='IGRF Bx')
    ax.plot(HEPD.time, HEPD.mag[:, 5], '--b', label='IGRF By')
    ax.plot(HEPD.time, HEPD.mag[:, 6], '--r', label='IGRF Bz')
    ax.plot(HEPD.time, mag_avg, '--y', label='IGRF|B|')
    fig.add_subplot(ax)
    plt.ylabel('Magnetic Field (nT)')
    plt.ylim(-60000, 60000)
    plt.legend(loc='upper center', ncol=7, prop={'size': 8})

## ADD FUNCTIONS


def plot_title(fig, start_time, end_time, HEPD, MEPD, plot_file_type):
    a = start_time.strftime('%Y/%m/%d %H:%M:%S')
    b = end_time.strftime('%H:%M:%S')
    c = MEPD.time[-1] - MEPD.time[0]

    orbit_no = MEPD.orbit_no
    if orbit_no < 10000 :
        title = 'Orbit: 0' + str(orbit_no) + '   Date: ' + str(a) + ' - ' \
            + str(b) + 'UT (' + str(c) + ' sec)'
        savename = 'isss/plots/ORB_0' + str(orbit_no)
    else:
        title = 'Orbit: ' + str(orbit_no) + '   Date: ' + str(a) + ' - ' \
            + str(b) + 'UT (' + str(c) + ' sec)'
        savename = 'isss/plots/ORB_' + str(orbit_no)
    fig.suptitle(title, fontsize=20)
    # Filetype
    if plot_file_type == 'PDF':
        plt.savefig(savename + '.pdf')
    elif plot_file_type == 'PNG':
        plt.savefig(savename + '.png')
    else:
        print('Error: Currently, only PDF and PNG files are supported.')
        exit()

