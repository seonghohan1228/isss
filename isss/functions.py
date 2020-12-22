# functions.py
import os
import warnings
import spacepy.datamodel as dm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from datetime import datetime
from file_functions import *
from plot_functions import *
from classes import *

def env_setup(conda_path):
    ''' Creates an environment suitable for isss program. '''
    os.environ["PROJ_LIB"] = conda_path;
    warnings.filterwarnings('ignore')

def data(orbit_no, data_path, show_avail, show_tree):
    ''' Gets the HEPD and MEPD data of the corresponding orbit number in the 
        data path in SpacePy datamodel format. Also capable of showing 
        available files and the selected file's file structure tree if 
        show_avail or show_tree is set to True, respectively. '''
    files = get_files(DATA_PATH, show_avail)
    HEPD_file_path, MEPD_file_path = get_file_paths(ORBIT_NO, DATA_PATH, files)
    print('Selected orbit number: ', str(ORBIT_NO), '\n\n')
    # Create SpacePy datamodel.
    HEPD_data = dm.fromHDF5(HEPD_file_path)
    MEPD_data = dm.fromHDF5(MEPD_file_path)
    if show_tree == True:
        file_tree(HEPD_data, MEPD_data)
    HEPD = isssData(orbit_no, HEPD_data, GROUP1, DATASET1, DATASET2)
    MEPD = isssData(orbit_no, MEPD_data, GROUP2, DATASET1, DATASET2)
    return HEPD, MEPD

def plot(HEPD, MEPD, pole, plot_name, plot_file_type):
    ''' Plots ISSS data into one figure and saves it into a desired file. 
        Currently, the file name can be chosen by user, but only PDF and PNG
        files are supported. '''
    fig = plt.figure(figsize=(20, 30))
    outer = gridspec.GridSpec(4, 2, wspace=0.1, hspace=0.3)
    # Convert UNIX datetime to UTC datetime
    start_time = datetime.fromtimestamp(MEPD.time[0])
    end_time = datetime.fromtimestamp(MEPD.time[-1])
    # Plot each data.
    plot_pc1(fig, outer, HEPD, MEPD)
    plot_mag(fig, outer, HEPD, MEPD)

    plot_title(fig, start_time, end_time, HEPD, MEPD, plot_file_type)
    plt.show()

    print('Pole: ', pole)

