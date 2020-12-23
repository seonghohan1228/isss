# plot_functions.py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
from datetime import datetime
from mpl_toolkits.basemap import Basemap
from data_functions import slice_AB, sliceAB2, B_avg, closest, geomag_lat, \
    slice_tel, start_end_time

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

def plot_pc1(fig, outer, HEPD, MEPD):
    ''' Plots PC1 data of HEPD and MEPD data onto figure. '''
    plot_msg('PC1', 'start')
    plot_msg('| HEPD PC1', 'start')
    inner = gridspec.GridSpecFromSubplotSpec(1, 2,
                    subplot_spec=outer[0], wspace=0.1, hspace=0.1)
    # Plot HEPD PC1 data.
    ax = plt.Subplot(fig, inner[0])
    ax.plot(HEPD.pc1, HEPD.time - HEPD.time[0], '-k')
    fig.add_subplot(ax)
    ax.set_title('HEPD: Time vs PC1')
    ax.set_xlabel('PC1')
    ax.set_ylabel('Time (sec)')
    plot_msg('| HEPD PC1', 'end')
    # Divide PC1 data into MEPD-A and MEPD-B
    MEPD_pc1_A, MEPD_pc1_B = sliceAB(MEPD.pc1, MEPD.dt, 3)
    MEPD_time_A, MEPD_time_B = sliceAB(MEPD.time, MEPD.dt, 3)
    # Plot MEPD-A and MEPD-B PC1 data.
    plot_msg('| MEPD PC1', 'start')
    ax = plt.Subplot(fig, inner[1])
    ax.plot(MEPD_pc1_A, MEPD_time_A - MEPD_time_A[0], '-k', label='MEPD-A')
    ax.plot(MEPD_pc1_B, MEPD_time_B - MEPD_time_B[0], '-r', label='MEPD-B')
    fig.add_subplot(ax)
    ax.set_title('MEPD: Time vs PC1')
    ax.set_xlabel('PC1')
    plt.legend()
    plot_msg('| HEPD PC1', 'end')
    plot_msg('PC1', 'end')

def plot_mag(fig, outer, data):
    ''' Plots magnetic field data of the given data (HEPD or MEPD). '''
    plot_msg('Magnetic field', 'start')
    inner = gridspec.GridSpecFromSubplotSpec(1, 1,
                    subplot_spec=outer[2], wspace=0.1, hspace=0.1)
    mag_avg = B_avg(data.mag)
    # Plot magnetic field data.
    ax = plt.Subplot(fig, inner[0])
    ax.plot(data.time, data.mag[:, 0], 'k', label='Bx')
    ax.plot(data.time, data.mag[:, 1], 'b', label='By')
    ax.plot(data.time, data.mag[:, 2], 'r', label='Bz')
    ax.plot(data.time, data.mag[:, 4], '--k', label='IGRF Bx')
    ax.plot(data.time, data.mag[:, 5], '--b', label='IGRF By')
    ax.plot(data.time, data.mag[:, 6], '--r', label='IGRF Bz')
    ax.plot(data.time, mag_avg, '--y', label='IGRF|B|')
    fig.add_subplot(ax)
    plt.ylabel('Magnetic Field (nT)')
    plt.ylim(-60000, 60000)
    plt.legend(loc='upper center', ncol=7, prop={'size': 8})
    plot_msg('Magnetic field', 'end')

def plot_geomag(ax, m, alt, start_time, conv_module):
    ''' Plots geomagnetic latitudes onto map using geomag_lat function. Conversion
        module can be selected. Either Spacepy or AACGMV2 is used. Also plots 
        terminator onto map. '''
    # Plot geomagnetic latitude.
    mat = geomag_lat(alt[0], start_time, conv_module)
    for i in range(5):
        x, y = m(np.arange(360) - 180, mat[i, :])
        m.plot(x, y, 'b', linewidth=1)
        if i < 2:
            ax.annotate(str(30 * i - 60) + 'S', (x[0], y[0]), color='b')
        elif i == 2:
            ax.annotate('0', (x[0], y[0]), color='b')
        elif i > 2:
            ax.annotate(str(30 * i - 60) + 'N', (x[0], y[0]), color='b')
    # Plot terminator
    m.nightshade(start_time)

def plot_pos(fig, outer, data, pole, conv_module):
    ''' Plots satellite position onto the map using basemap. The plots are in 
        Mercador and orthographic projection. The pole of interest is given as 
        and argument by the variable 'pole'. 'pole' is 90 for north pole 
        projection, and -90 for south pole projection in the orthographic 
        projection plot. Either MEPD or HEPD data can be used. '''
    plot_msg('Satellite position', 'start')
    # Position data
    lat = data.pos[:, 0]
    lon = data.pos[:, 1]
    alt = data.pos[:, 2]
    start_time, end_time = start_end_time(data.time)

    # Mercador projection
    plot_msg('| Mercador projection', 'start')
    inner = gridspec.GridSpecFromSubplotSpec(1, 1,
                    subplot_spec=outer[1], wspace=0.1, hspace=0.1)
    ax = plt.Subplot(fig, inner[0])
    m = Basemap(projection='merc',llcrnrlat=-85,urcrnrlat=85, llcrnrlon=-180,\
            urcrnrlon=180, ax=ax)
    m.drawcoastlines()
    m.drawmeridians(np.arange(0,360,45), labels=[False, False, False, True])
    m.drawparallels(np.arange(-90,90,30), labels=[False, True, False, False])
    # Plot satellite position.
    X, Y = m(lon, lat)
    m.scatter(X, Y, marker='.', c=data.time, cmap=plt.get_cmap('jet'))
    ax.annotate(start_time.strftime('%H:%M'), (X[0], Y[0]))
    ax.annotate(end_time.strftime('%H:%M'), (X[-1], Y[-1]))
    #Plot geomagnetic latitude.
    plot_msg('| | Geomagnetic latitude', 'start')
    plot_geomag(ax, m, alt, start_time, conv_module)
    plot_msg('| | Geomagnetic latitude', 'end')
    ax.set_title('Orbit (Mercador projection)')
    fig.add_subplot(ax)
    plot_msg('| Mercador projection', 'end')

    # Orthographic projection.
    plot_msg('| Orthographic projection', 'start')
    inner = gridspec.GridSpecFromSubplotSpec(1, 1,
                    subplot_spec=outer[3], wspace=0.1, hspace=0.1)
    ax = plt.Subplot(fig, inner[0])
    if pole == 'north':
        p = 90
    elif pole == 'south':
        p = -90
    else:
        print('Error: pole should be either north or south')
        exit()
    m = Basemap(projection='ortho', lat_0=p, lon_0=0, ax=ax)
    m.drawcoastlines()
    m.drawmeridians(np.arange(0,360,45), labels=[False, False, False, True])
    # Cannot label parallels on Orthographic basemap.
    m.drawparallels(np.arange(-90,90,30))
    # Plot satellite position.
    X, Y = m(lon, lat)
    m.scatter(X, Y, marker='.', c=data.time, cmap=plt.get_cmap('jet'))
    ax.annotate(start_time.strftime('%H:%M'), (X[0], Y[0]))
    ax.annotate(end_time.strftime('%H:%M'), (X[-1], Y[-1]))
    #cbar = fig.colorbar(s1, ax=ax2, label='Time (UNIX)')
    # Plot geomagnetic latitude.
    plot_msg('| | Geomagnetic latitude', 'start')
    plot_geomag(ax, m, alt, start_time, conv_module)
    plot_msg('| | Geomagnetic latitude', 'end')
    ax.set_title('Orbit (Orthographic projection)')
    fig.add_subplot(ax)
    plot_msg('| Orthographic projection', 'end')
    plot_msg('Satellite position', 'end')

def nplot(fig, inner, data, xmin, xmax, label, xlabel, ylabel1, ylabel2):
    ''' Plots n plots vertically (stacked). lable and ylabel can be entered by 
        setting variables. Depending on the number of subplots, ylabel1 and 
        ylabel2 can either be combined or separated. The number of plots are 
        automatically set to the length of data, which should be an array. 
        xlabel can be set to True to show values on the x axis, or False to 
        hide values. '''
    n = len(data)
    for i in range(n):
        ax = plt.Subplot(fig, inner[i])
        im = ax.imshow(X=np.transpose(data[i]), origin='lower', cmap=new_cmap(), 
            aspect='auto', interpolation='none', extent = [xmin, xmax, 0, 400])
        ax.text(0.02, 0.85, (label + ' ' + str(i)), 
            horizontalalignment='left', transform=ax.transAxes)
        if (i != n-1) or (xlabel == False):
            plt.setp(ax.get_xticklabels(), visible=False)
        # ylabel center allignment.
        if (n == 3) and (i == 1):
            ax.set_ylabel(ylabel1 + ylabel2)
        fig.add_subplot(ax)
    fig.subplots_adjust(left=0.08, right=0.9, bottom=0.01, top=0.99, 
                    wspace=0.0, hspace=0.0)
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    #cb_ax = fig.add_axes([0.92, 0.05, 0.02, 0.9])
    #cbar = fig.colorbar(im, cax=cb_ax)

def plot_tel(fig, outer, data):
    ''' Plots HEPD telescope data. '''
    plot_msg('HEPD telescope', 'start')
    # Proton and electron data.
    tel = data.tel
    p, e = slice_tel(tel)
    start_time, end_time = start_end_time(data.time)
    xmin = mdates.date2num(start_time)
    xmax = mdates.date2num(end_time)
    # Plot HEPD proton data
    plot_msg('| HEPD proton', 'start')
    inner = gridspec.GridSpecFromSubplotSpec(3, 1,
                    subplot_spec=outer[4], wspace=0.05, hspace=0.01)
    nplot(fig, inner, p, xmin, xmax, 'Telescope', False, 'HEPD Proton ', 
        'Energy [MeV]')
    plot_msg('| HEPD proton', 'end')
    # Plot HEPD electron data
    plot_msg('| HEPD electron', 'start')
    inner = gridspec.GridSpecFromSubplotSpec(3, 1,
                    subplot_spec=outer[6], wspace=0.05, hspace=0.01)
    nplot(fig, inner, e, xmin, xmax, 'Telescope', True, 'HEPD Electron ', 
        'Energy [MeV]')
    plot_msg('| HEPD electron', 'end')
    plot_msg('HEPD telescope', 'end')

def plot_det(fig, outer, data):
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
    time_A, time_B = sliceAB(data.time, data.dt, 3)
    # MEPD-A
    plot_msg('| MEPD-A detector', 'start')
    start_time, end_time = start_end_time(time_A)
    xmin = mdates.date2num(start_time)
    xmax = mdates.date2num(end_time)
    inner = gridspec.GridSpecFromSubplotSpec(4, 1,
                    subplot_spec=outer[5], wspace=0.05, hspace=0.01)
    nplot(fig, inner, det_A, xmin, xmax, 'Detector', False, 'MEPD-A ', 
        'Energy [keV]')
    plot_msg('| MEPD-A detector', 'end')
    # MEPD-B
    plot_msg('| MEPD-B detector', 'start')
    start_time, end_time = start_end_time(time_B)
    xmin = mdates.date2num(start_time)
    xmax = mdates.date2num(end_time)
    inner = gridspec.GridSpecFromSubplotSpec(4, 1,
                    subplot_spec=outer[7], wspace=0.05, hspace=0.01)
    nplot(fig, inner, det_B, xmin, xmax, 'Detector', True, 'MEPD-B ', 
        'Energy [keV]')
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
    fig.suptitle(title, fontsize=20)
    # Filetype
    if plot_file_type == 'pdf':
        plt.savefig(savename + '.pdf')
    elif plot_file_type == 'png':
        plt.savefig(savename + '.png')
    else:
        print('Error: Currently, only PDF and PNG files are supported.')
        exit()

