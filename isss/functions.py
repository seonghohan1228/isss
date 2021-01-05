# functions.py
import os
import warnings
import spacepy.datamodel as dm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from file_functions import *
from plot_functions import plot_msg, plot_pc1, plot_mag, plot_pos, plot_tel, \
    plot_det, plot_title
from classes import isssData
from constants import *

def env_setup():
    ''' Creates an environment suitable for isss program. '''
    os.environ["PROJ_LIB"] = CONDA_PATH
    warnings.filterwarnings('ignore')

def data(orbit_no, show_avail, show_tree):
    ''' Gets the HEPD and MEPD data of the corresponding orbit number in the 
        data path in SpacePy datamodel format. Also capable of showing 
        available files and the selected file's file structure tree if 
        show_avail or show_tree is set to True, respectively. '''
    files = get_files(DATA_PATH, show_avail)
    HEPD_file_path, MEPD_file_path = get_file_paths(orbit_no, DATA_PATH, files)
    print('Selected orbit number: ', str(orbit_no), '\n')
    # Create SpacePy datamodel.
    HEPD_data = dm.fromHDF5(HEPD_file_path)
    MEPD_data = dm.fromHDF5(MEPD_file_path)
    if show_tree == True:
        file_tree(HEPD_data, MEPD_data)
    HEPD = isssData(orbit_no, HEPD_data, GROUP1, DATASET1, DATASET2)
    MEPD = isssData(orbit_no, MEPD_data, GROUP2, DATASET1, DATASET2)
    return HEPD, MEPD

def plot(HEPD, MEPD, pole, conv_module, plot_name, plot_file_type):
    ''' Plots ISSS data into one figure and saves it into a desired file. 
        Currently, the file name can be chosen by user, but only PDF and PNG
        files are supported. '''
    plot_msg('Combined graph', 'start')
    print('')
    fig = plt.figure(dpi=200)
    fig.set_size_inches(WIDTH, HEIGHT) 
    outer_grid = fig.add_gridspec(4, 2, wspace=0.2, hspace=0.2)
    # Plot each data.
    plot_pc1(fig, outer_grid, HEPD, MEPD)
    plot_mag(fig, outer_grid, HEPD)
    plot_pos(fig, outer_grid, MEPD, pole, conv_module, mag=False)
    plot_tel(fig, outer_grid, HEPD)
    plot_det(fig, outer_grid, MEPD)

    plot_title(fig, MEPD, plot_file_type)
    print('')
    plot_msg('Combined graph', 'end')
    #plt.show()

