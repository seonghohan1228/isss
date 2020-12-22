# main.py

#### MODULES ############################################################

import matplotlib.dates as mdates
from mpl_toolkits.basemap import Basemap
from datetime import datetime

from functions import *
from classes import *
from constants import *
from data_functions import *
'''
# Plot graphs
def plot_graph(orbit_no, HEPD_time,HEPD_pc1, HEPD_pos, HEPD_mag, tel0, tel1, tel2,
            MEPD_time, dt, MEPD_pc1, MEPD_pos, MEPD_mag, det0, det1, det2, det3, pole, filetype):

    ## Satellite position
    # Mercador projection
    inner = gridspec.GridSpecFromSubplotSpec(1, 1,
                    subplot_spec=outer[1], wspace=0.1, hspace=0.1)

    lat = MEPD_pos[:, 0]
    lon = MEPD_pos[:, 1]
    alt = MEPD_pos[:, 2]
    
    start_time = datetime.fromtimestamp(MEPD_time[0]) # Converts UNIX datetime to UTC datetime
    end_time = datetime.fromtimestamp(MEPD_time[-1])

    ax = plt.Subplot(fig, inner[0])
    
    m = Basemap(projection='merc',llcrnrlat=-85,urcrnrlat=85, llcrnrlon=-180,urcrnrlon=180, ax=ax)
    m.drawcoastlines()
    m.drawmeridians(np.arange(0,360,45), labels=[False, False, False, True])
    m.drawparallels(np.arange(-90,90,30), labels=[False, True, False, False])

    # Plot satellite position.
    X, Y = m(lon, lat)
    m.scatter(X, Y, marker='.', c=MEPD_time, cmap=plt.get_cmap('jet'))
    ax.annotate(start_time.strftime('%H:%M'), (X[0], Y[0]))
    ax.annotate(end_time.strftime('%H:%M'), (X[-1], Y[-1]))

    # Plot geomagnetic latitude.
    mat = geo_lat(alt[0], start_time)
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

    ax.set_title('Orbit (Mercador projection)')

    fig.add_subplot(ax)

    # Orthographic projection.
    inner = gridspec.GridSpecFromSubplotSpec(1, 1,
                    subplot_spec=outer[3], wspace=0.1, hspace=0.1)
    
    ax = plt.Subplot(fig, inner[0])

    m = Basemap(projection='ortho', lat_0=pole, lon_0=0, ax=ax)
    m.drawcoastlines()
    m.drawmeridians(np.arange(0,360,45), labels=[False, False, False, True])
    m.drawparallels(np.arange(-90,90,30)) # Cannot label parallels on Orthographic basemap.
    
    # Plot satellite position.
    X, Y = m(lon, lat)
    m.scatter(X, Y, marker='.', c=MEPD_time, cmap=plt.get_cmap('jet'))
    ax.annotate(start_time.strftime('%H:%M'), (X[0], Y[0]))
    ax.annotate(end_time.strftime('%H:%M'), (X[-1], Y[-1]))
    
    #cbar = fig.colorbar(s1, ax=ax2, label='Time (UNIX)')

    # Plot geomagnetic latitude.
    mat = geo_lat(alt[0], start_time)
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

    ax.set_title('Orbit (Orthographic projection)')

    fig.add_subplot(ax)

    ## Telescope
    p0, p1, p2, e0, e1, e2 = div_tel(tel0, tel1, tel2)
    
    # Plot HEPD proton data
    inner = gridspec.GridSpecFromSubplotSpec(3, 1,
                    subplot_spec=outer[4], wspace=0.05, hspace=0.01)

    xmin = mdates.date2num(datetime.fromtimestamp(HEPD_time[0]))
    xmax = mdates.date2num(datetime.fromtimestamp(HEPD_time[-1]))

    p = [p0, p1, p2]
    for i in range(3):
        ax = plt.Subplot(fig, inner[i])
        ax.imshow(X=np.transpose(p[i]), origin='lower', cmap=new_cmap(), aspect='auto',
            interpolation='none', extent = [xmin, xmax, 0, 400])
        ax.text(0.02, 0.85, ('Telescope ' + str(i)), horizontalalignment='left', transform=ax.transAxes)
        if i == 1:
            ax.set_ylabel('HEPD Proton Energy [MeV]')
        fig.add_subplot(ax)

    fig.subplots_adjust(left=0.08, right=0.9, bottom=0.01, top=0.99, wspace=0.0, hspace=0.0)
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    #cb_ax_p = fig_p.add_axes([0.92, 0.05, 0.02, 0.9])
    #cbar_p = fig_p.colorbar(im_p, cax=cb_ax_p)

    # Plot HEPD electron data
    inner = gridspec.GridSpecFromSubplotSpec(3, 1,
                    subplot_spec=outer[6], wspace=0.05, hspace=0.01)

    e = [e0, e1, e2]
    for i in range(3):
        ax = plt.Subplot(fig, inner[i])
        ax.imshow(X=np.transpose(e[i]), origin='lower', cmap=new_cmap(), aspect='auto', 
            interpolation='none', extent = [xmin, xmax, 0, 400])
        ax.text(0.02, 0.85, ('Telescope ' + str(i)), horizontalalignment='left', transform=ax.transAxes)
        if i == 1:
            ax.set_ylabel('HEPD Electron Energy [MeV]')
        fig.add_subplot(ax)
        
    fig.subplots_adjust(left=0.08, right=0.9, bottom=0.01, top=0.99, wspace=0.0, hspace=0.0)
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    #cb_ax_e = fig_e.add_axes([0.92, 0.05, 0.02, 0.9])
    #cbar_e = fig_e.colorbar(im_e, cax=cb_ax_e)


    ## Detector
    # Divide PC1 data into MEPD-A and MEPD-B
    det0_A, det0_B = sliceAB2(det0, dt, 3)
    det1_A, det1_B = sliceAB2(det1, dt, 3)
    det2_A, det2_B = sliceAB2(det2, dt, 3)
    det3_A, det3_B = sliceAB2(det3, dt, 3)
    MEPD_time_A, MEPD_time_B = sliceAB(MEPD_time, dt, 3)

    # MEPD-A
    inner = gridspec.GridSpecFromSubplotSpec(4, 1,
                    subplot_spec=outer[5], wspace=0.05, hspace=0.01)
    
    # Plot MEPD-A
    xmin = mdates.date2num(datetime.fromtimestamp(MEPD_time_A[0]))
    xmax = mdates.date2num(datetime.fromtimestamp(MEPD_time_A[-1]))
    
    det = [det0_A, det1_A, det2_A, det3_A]
    for i in range(4):
        ax = plt.Subplot(fig, inner[i])
        ax.imshow(X=np.transpose(det[i]), origin='lower',cmap=new_cmap(), 
            aspect='auto', interpolation='none', extent = [xmin, xmax, 0, 400])
        ax.text(0.02, 0.85, ('Detector ' + str(i)), horizontalalignment='left', transform=ax.transAxes)
        if i == 1:
            ax.set_ylabel('Energy [keV]')
        if i == 2:
            ax.set_ylabel('MEPD-A')
        fig.add_subplot(ax)
    
    fig.subplots_adjust(left=0.08, right=0.9, bottom=0.01, top=0.99, wspace=0.0, hspace=0.0)
    #cb_ax = fig.add_axes([0.92, 0.05, 0.02, 0.9])
    #cbar = fig.colorbar(im, cax=cb_ax)
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))


    # MEPD-B
    inner = gridspec.GridSpecFromSubplotSpec(4, 1,
                    subplot_spec=outer[7], wspace=0.05, hspace=0.01)
    
    # Plot MEPD-B
    det = [det0_B, det1_B, det2_B, det3_B]
    for i in range(4):
        ax = plt.Subplot(fig, inner[i])
        ax.imshow(X=np.transpose(det[i]), origin='lower',cmap=new_cmap(), 
            aspect='auto', interpolation='none', extent = [xmin, xmax, 0, 400])
        ax.text(0.02, 0.85, ('Detector ' + str(i)), horizontalalignment='left', transform=ax.transAxes)
        if i == 1:
            ax.set_ylabel('Energy [keV]')
        if i == 2:
            ax.set_ylabel('MEPD-B')
        fig.add_subplot(ax)

    fig.subplots_adjust(left=0.08, right=0.9, bottom=0.01, top=0.99, wspace=0.0, hspace=0.0)
    #cb_ax = fig.add_axes([0.92, 0.05, 0.02, 0.9])
    #cbar = fig.colorbar(im, cax=cb_ax)
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))


    
'''
####################################################################################

def main():
    env_setup(CONDA_PATH)
    HEPD, MEPD = data(ORBIT_NO, DATA_PATH, True, True)
    plot(HEPD, MEPD, SOUTH_POLE, 'plot_combined', PDF)

if __name__ == '__main__':
    main()
