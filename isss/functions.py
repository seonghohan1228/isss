# functions.py
import os
import warnings
import spacepy.datamodel as dm
from file_functions import *
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



def plot(HEPD, MEPD, pole, plot_file_type):
    print('Orbit Number: ', HEPD.orbit_no)
    print(HEPD.time)
    print(MEPD.pc1)
    print('Pole: ', pole)
    print('Plot file type: ', plot_file_type)



