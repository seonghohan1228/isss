# plot_functions.py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
from datetime import datetime
from mpl_toolkits.basemap import Basemap
from data_functions import sliceAB, sliceAB2, B_avg, closest, geomag_lat, \
    slice_tel, start_end_time
from plot_constants import *

def new_cmap():
    ''' Creates custom colormap for HEPD telescope and MEPD detector plot. '''
    jet = plt.cm.get_cmap('jet', 256)
    newcolors = jet(np.linspace(0, 1, 256))
    white = np.array([256/256, 256/256, 256/256, 1])
    newcolors[0, :] = white
    new_cmap = matplotlib.colors.ListedColormap(newcolors)
    return new_cmap

def plot_msg(name, state):
    ''' Prints out message about the state of plot. Either 'start' or 'end' is 
        given to state. '''
    if state == 'start':
        print(name, ' plotting...')
    elif state == 'end':
        print(name, ' plot completed')
    else:
        print('Error: plot_msg() state is invalid')
        exit()

def make_ticks(time, lat, lon, data_num):
    ''' Creates ticks for UT, GLAT, GLONG. The number of ticks are given by 
        constant TICK_NUM in plot_constants.py. '''
    ticks = []
    labels = []
    time_data_num = len(time)
    pos_data_num = len(lat)
    time_interval = round(time_data_num / TICK_NUM)
    pos_interval = round(pos_data_num / TICK_NUM)
    data_interval = round(data_num / TICK_NUM)
    for i in range(TICK_NUM):
        time_index = time_interval * i
        pos_index = pos_interval * i
        data_index = data_interval * i
        t = datetime.fromtimestamp(time[time_index])
        if i == 0:
            ut = 'UT'
            glat = 'GLAT'
            glong = 'GLONG'
        else:
            ut = t.strftime('%H:%M')
            glat = str('{:.1f}'.format(lat[pos_index]))
            glong = str('{:.1f}'.format(lat[pos_index]))
        ticks.append(data_index)
        # '\n' for multiline ticks.
        labels.append(ut + '\n' + glat + '\n' + glong)
    return ticks, labels

def plot_pc1(fig, outer_grid, HEPD, MEPD):
    ''' Plots PC1 data of HEPD and MEPD data onto figure. '''
    plot_msg('PC1', 'start')
    plot_msg('| HEPD PC1', 'start')
    inner_grid = outer_grid[0, 0].subgridspec(1, 2, wspace=0.2, hspace=0)
    axes = inner_grid.subplots()
    # Plot HEPD PC1 data.
    axes[0].plot(HEPD.pc1, HEPD.time - HEPD.time[0], '-k')
    axes[0].set_title('HEPD: Time vs PC1', fontsize=AXES_FONT)
    axes[0].set_xlabel('PC1', fontsize=AXES_FONT)
    axes[0].set_ylabel('Time (sec)', fontsize=AXES_FONT)
    axes[0].tick_params (axis = 'x', direction='in', labelsize =TICK_FONT)
    axes[0].tick_params (axis = 'y', direction='in', labelsize =TICK_FONT)
    plot_msg('| HEPD PC1', 'end')
    # Divide PC1 data into MEPD-A and MEPD-B
    MEPD_pc1_A, MEPD_pc1_B = sliceAB(MEPD.pc1, MEPD.dt, 3)
    MEPD_time_A, MEPD_time_B = sliceAB(MEPD.time, MEPD.dt, 3)
    # Plot MEPD-A and MEPD-B PC1 data.
    plot_msg('| MEPD PC1', 'start')
    axes[1].plot(MEPD_pc1_A, MEPD_time_A - MEPD_time_A[0], '-k', label='MEPD-A')
    axes[1].plot(MEPD_pc1_B, MEPD_time_B - MEPD_time_B[0], '-r', label='MEPD-B')
    axes[1].set_title('MEPD: Time vs PC1', fontsize=AXES_FONT)
    axes[1].set_xlabel('PC1', fontsize=AXES_FONT)
    axes[1].tick_params (axis = 'x', direction='in', labelsize =TICK_FONT)
    axes[1].tick_params (axis = 'y', direction='in', labelsize =TICK_FONT)
    plt.legend(fontsize=TICK_FONT)
    plot_msg('| HEPD PC1', 'end')
    plot_msg('PC1', 'end')

def plot_mag(fig, outer_grid, data):
    ''' Plots magnetic field data of the given data (HEPD or MEPD). '''
    plot_msg('Magnetic field', 'start')
    inner_grid = outer_grid[1, 0].subgridspec(1, 1, wspace=0, hspace=0)
    ax = inner_grid.subplots()
    time = data.time
    mag = data.mag
    mag_avg = B_avg(data.mag)
    lat = data.pos[:, 0]
    lon = data.pos[:, 1]
    # Plot magnetic field data.
    ax.plot(time, mag[:, 0], 'k', label='Bx')
    ax.plot(time, mag[:, 1], 'b', label='By')
    ax.plot(time, mag[:, 2], 'r', label='Bz')
    ax.plot(time, mag[:, 4], '--k', label='IGRF Bx')
    ax.plot(time, mag[:, 5], '--b', label='IGRF By')
    ax.plot(time, mag[:, 6], '--r', label='IGRF Bz')
    ax.plot(time, mag_avg, '--y', label='IGRF|B|')
    ax.set_ylabel('Magnetic Field (nT)', fontsize=AXES_FONT)
    ax.tick_params (axis = 'x', direction='in', labelsize =TICK_FONT)
    ax.tick_params (axis = 'y', direction='in', labelsize =TICK_FONT)
    plt.ylim(-60000, 60000)
    plt.xlim(time[-1], time[0])
    ticks, label = make_ticks(time, lat, lon, time[-1] - time[0])
    ax.set_xticks(ticks + time[0])
    ax.set_xticklabels(label)
    plt.setp(ax.get_xticklabels(), visible=False, fontsize=TICK_FONT)    
    plt.legend(loc='upper center', ncol=7, prop={'size': 4})
    plot_msg('Magnetic field', 'end')

def plot_geomag(ax, m, alt, start_time, conv_module):
    ''' Plots geomagnetic latitudes onto map using geomag_lat function. 
        Conversion module can be selected. Either Spacepy or AACGMV2 is used. 
        Also plots terminator onto map. '''
    # Plot geomagnetic latitude.
    mat = geomag_lat(alt[0], start_time, conv_module)
    for i in range(5):
        x, y = m(np.arange(360) - 180, mat[i, :])
        m.plot(x, y, 'b', linewidth=1)
        if i < 2:
            ax.annotate(str(30 * i - 60) + 'S', (x[0], y[0]), color='b', 
                    fontsize=TICK_FONT)
        elif i == 2:
            ax.annotate('0', (x[0], y[0]), color='b', fontsize=TICK_FONT)
        elif i > 2:
            ax.annotate(str(30 * i - 60) + 'N', (x[0], y[0]), color='b', 
                    fontsize=TICK_FONT)
    # Plot terminator
    m.nightshade(start_time)

def plot_pos(fig, outer_grid, data, pole, conv_module, mag=True):
    ''' Plots satellite position onto the map using basemap. The plots are in 
        Mercador and orthographic projection. The pole of interest is given as 
        and argument by the variable 'pole'. 'pole' is 90 for north pole 
        projection, and -90 for south pole projection in the orthographic 
        projection plot. Either MEPD or HEPD data can be used. To save time, 
        argument 'mag' can be set to False in order to plot the position 
        without geomagnetic latitudes. It is set to True by default. '''
    plot_msg('Satellite position', 'start')
    # Position data
    lat = data.pos[:, 0]
    lon = data.pos[:, 1]
    alt = data.pos[:, 2]
    start_time, end_time = start_end_time(data.time)
    # Mercador projection
    plot_msg('| Mercador projection', 'start')
    inner_grid = outer_grid[0:2, 1].subgridspec(2, 1, wspace=0, hspace=0.2)
    axes = inner_grid.subplots()
    m = Basemap(projection='merc',llcrnrlat=-85,urcrnrlat=85, llcrnrlon=-180,\
            urcrnrlon=180, ax=axes[0])
    m.drawcoastlines()
    m.drawmeridians(np.arange(0,360,45), labels=[False, False, False, True], 
                fontsize=TICK_FONT)
    m.drawparallels(np.arange(-90,90,30), labels=[False, True, False, False], 
                fontsize=TICK_FONT)
    # Plot satellite position.
    X, Y = m(lon, lat)
    m.scatter(X, Y, marker='.', c=data.time, cmap=plt.get_cmap('jet'))
    axes[0].annotate(start_time.strftime('%H:%M'), (X[0], Y[0]), 
                fontsize=AXES_FONT)
    axes[0].annotate(end_time.strftime('%H:%M'), (X[-1], Y[-1]), 
                fontsize=AXES_FONT)
    #Plot geomagnetic latitude.
    if mag == True:
        plot_msg('| | Geomagnetic latitude', 'start')
        plot_geomag(axes[0], m, alt, start_time, conv_module)
        plot_msg('| | Geomagnetic latitude', 'end')
    axes[0].set_title('Orbit (Mercador projection)', fontsize=AXES_FONT)
    plot_msg('| Mercador projection', 'end')
    # Orthographic projection.
    plot_msg('| Orthographic projection', 'start')
    if pole == 'north':
        p = 90
    elif pole == 'south':
        p = -90
    else:
        print('Error: pole should be either north or south')
        exit()
    m = Basemap(projection='ortho', lat_0=p, lon_0=0, ax=axes[1])
    m.drawcoastlines()
    m.drawmeridians(np.arange(0,360,45), labels=[False, False, False, True], 
                fontsize=TICK_FONT)
    # Cannot label parallels on Orthographic basemap.
    m.drawparallels(np.arange(-90,90,30))
    # Plot satellite position.
    X, Y = m(lon, lat)
    m.scatter(X, Y, marker='.', c=data.time, cmap=plt.get_cmap('jet'))
    axes[1].annotate(start_time.strftime('%H:%M'), (X[0], Y[0]), 
                fontsize=AXES_FONT)
    axes[1].annotate(end_time.strftime('%H:%M'), (X[-1], Y[-1]), 
                fontsize=AXES_FONT)
    #cbar = fig.colorbar(s1, ax=ax2, label='Time (UNIX)')
    # Plot geomagnetic latitude.
    if mag == True:
        plot_msg('| | Geomagnetic latitude', 'start')
        plot_geomag(axes[1], m, alt, start_time, conv_module)
        plot_msg('| | Geomagnetic latitude', 'end')
    axes[1].set_title('Orbit (Orthographic projection)', fontsize=AXES_FONT)
    plot_msg('| Orthographic projection', 'end')
    plot_msg('Satellite position', 'end')

def yticks(axes, data):
    ''' Plots appropriate yticks for HEPD telescope and MEPD detector plot. 
        Set 'data' to 'p', 'e', or 'det' for proton, electron, and detector 
        plot, respectively. '''
    # Telescope
    if data == 'p':
        labels = [3, 7, 11, 14, 20, 24, 38]
        ticks = list(range(len(labels)))
    elif data == 'e':
        labels = [0.35, 0.45, 0.56, 0.66, 0.77, 0.87, 0.98, 1.08, 1.19, 1.32, 1.52, 
                2, 2.5, 3, 3.5, 4.5]
        ticks = list(range(len(labels)))
    # Detector
    elif data == 'det':
        labels = [0, 100, 200, 300, 400]
        ticks = np.linspace(1, 64, 5)
    else:
        print('Error: invalid data to create yticks')
        exit()
    for ax in axes:
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)

def nplot(fig, axes, data, time, lat, lon, label, xlabel, ylabel1, ylabel2, 
        plot_data):
    ''' Plots n plots vertically (stacked). lable and ylabel can be entered by 
        setting variables. Depending on the number of subplots, ylabel1 and 
        ylabel2 can either be combined or separated. The number of plots are 
        automatically set to the length of data, which should be an array. 
        xlabel can be set to True to show values on the x axis, or False to 
        hide values. '''
    n = len(data)
    ticks, labels = make_ticks(time, lat, lon, len(data[0]))
    for i in range(n):
        axes[i].imshow(X=np.transpose(data[i]), origin='lower', 
            cmap=new_cmap(), aspect='auto', interpolation='none')
        axes[i].set_xticks(ticks)
        axes[i].tick_params(axis = 'x', direction='in', labelsize =TICK_FONT)
        axes[i].tick_params(axis = 'y', direction='in', labelsize =TICK_FONT)
        axes[i].text(0.02, 0.85, (label + ' ' + str(i)), 
            horizontalalignment='left', transform=axes[i].transAxes, 
            fontsize=TICK_FONT)
        if i != n-1:
            plt.setp(axes[i].get_xticklabels(), visible=False)
    # ylabel center allignment.
    if n == 3:
        axes[1].set_ylabel(ylabel1 + ylabel2, fontsize=AXES_FONT)
    elif n == 4:
        axes[1].set_ylabel(ylabel2, fontsize=AXES_FONT)
        axes[2].set_ylabel(ylabel1, fontsize=AXES_FONT)
    # Show xtick labels on the bottom-most plot.
    axes[n-1].set_xticklabels(labels)
    plt.setp(axes[n-1].get_xticklabels(), visible=xlabel, fontsize=TICK_FONT)
    # Set yticks.
    yticks(axes, plot_data)
    #cb_ax = fig.add_axes([0.92, 0.05, 0.02, 0.9])
    #cbar = fig.colorbar(im, cax=cb_ax)

def plot_tel(fig, outer_grid, data):
    ''' Plots HEPD telescope data. '''
    plot_msg('HEPD telescope', 'start')
    # Proton and electron data.
    tel = data.tel
    p, e = slice_tel(tel)
    time = data.time
    lat = data.pos[:, 0]
    lon = data.pos[:, 1]
    # Plot HEPD proton data
    plot_msg('| HEPD proton', 'start')
    inner_grid_p = outer_grid[2, 0].subgridspec(3, 1, wspace=0, hspace=0)
    axes_p = inner_grid_p.subplots()
    nplot(fig, axes_p, p, time, lat, lon, 'Telescope', False, 
        'HEPD Proton ', 'Energy [MeV]', 'p')
    plot_msg('| HEPD proton', 'end')
    # Plot HEPD electron data
    plot_msg('| HEPD electron', 'start')
    inner_grid_e = outer_grid[3, 0].subgridspec(3, 1, wspace=0, hspace=0)
    axes_e = inner_grid_e.subplots()
    nplot(fig, axes_e, e, time, lat, lon, 'Telescope', True, 
        'HEPD Electron ', 'Energy [MeV]', 'e')
    plot_msg('| HEPD electron', 'end')
    plot_msg('HEPD telescope', 'end')

def plot_det(fig, outer_grid, data):
    ''' Plots detector data. '''
    plot_msg('MEPD detector', 'start')
    # Divide detector data into MEPD-A and MEPD-B based on dt.
    det = data.det
    det_A = []
    det_B = []
    for i in range(len(det)):
        A, B = sliceAB2(det[i], data.dt, 3)
        det_A.append(A)
        det_B.append(B)
    lat = data.pos[:, 0]
    lon = data.pos[:, 1]
    time_A, time_B = sliceAB(data.time, data.dt, 3)
    # MEPD-A
    plot_msg('| MEPD-A detector', 'start')
    inner_grid_A = outer_grid[2, 1].subgridspec(4, 1, wspace=0, hspace=0)
    axes_A = inner_grid_A.subplots()
    nplot(fig, axes_A, det_A, time_A, lat, lon, 'Detector', False, 
        'MEPD-A ', 'Energy [keV]', 'det')
    plot_msg('| MEPD-A detector', 'end')
    # MEPD-B
    plot_msg('| MEPD-B detector', 'start')
    inner_grid_B = outer_grid[3, 1].subgridspec(4, 1, wspace=0, hspace=0)
    axes_B = inner_grid_B.subplots()
    nplot(fig, axes_B, det_B, time_B, lat, lon, 'Detector', True, 
        'MEPD-B ', 'Energy [keV]', 'det')
    plot_msg('| MEPD-B detector', 'end')
    plot_msg('MEPD detector', 'end')

def plot_title(fig, data, plot_file_type):
    ''' Creates and plots title and saves the file in the desired plot file 
        type. The file type can be selected from PDF and PNG. '''
    start_time, end_time = start_end_time(data.time)
    a = start_time.strftime('%Y/%m/%d %H:%M:%S')
    b = end_time.strftime('%H:%M:%S')
    c = data.time[-1] - data.time[0]
    # Creating title
    orbit_no = data.orbit_no
    if orbit_no < 10000 :
        title = 'Orbit: 0' + str(orbit_no) + '   Date: ' + str(a) + ' - ' \
            + str(b) + 'UT (' + str(c) + ' sec)'
        savename = 'isss/plots/ORB_0' + str(orbit_no)
    else:
        title = 'Orbit: ' + str(orbit_no) + '   Date: ' + str(a) + ' - ' \
            + str(b) + 'UT (' + str(c) + ' sec)'
        savename = 'isss/plots/ORB_' + str(orbit_no)
    fig.suptitle(title, fontsize=TITLE_FONT)
    # Filetype
    if plot_file_type == 'pdf':
        plt.savefig(savename + '.pdf')
    elif plot_file_type == 'png':
        plt.savefig(savename + '.png')
    else:
        print('Error: Currently, only PDF and PNG files are supported.')
        exit()

