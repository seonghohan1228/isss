# file_functions.py

import os
import warnings

def env_setup(conda_path):
    os.environ["PROJ_LIB"] = conda_path;
    warnings.filterwarnings('ignore')

def show_avail_file():
    """ This function shows available files in the data directory. """
    files = os.listdir(os.getcwd() + '../data/')
    print('\nAvailable files:')
    for filename in files:
        print(filename)
    print('\n')
    