# file_functions.py

import os
import warnings

def env_setup(conda_path):
    ''' This function sets the environment for proper results. '''
    os.environ["PROJ_LIB"] = conda_path;
    warnings.filterwarnings('ignore')

def show_avail_file(data_path):
    ''' This function shows available files in the data directory. '''
    files = os.listdir(data_path)
    print('\nAvailable files:')
    for filename in files:
        print('- ', filename)
    print('\n')
    return files

def get_file_paths(orbit_no, data_path, files):
    ''' This function returns the HEPD and MEPD filenames of the corresponding
        orbit number. '''
    count = 0
    # Check for invalid orbit number.
    if orbit_no < 0:
        print('Error: Invalid orbit number.')
        exit()
    # Convert orbit number string (for filename).
    if orbit_no < 10000:
        orbit_no = '0' + str(orbit_no)
    else:
        orbit_no = str(orbit_no)
    # Get file names for HEPD and MEPD.
    for filename in files:
        if filename[27:32] == orbit_no:
            if filename[0:4] == 'HEPD':
                hepd_file_path = data_path + filename
                count += 1
                continue
            elif filename[0:4] == 'MEPD':
                mepd = data_path + filename
                count += 1
                continue
        # No match found.
        if filename == files[-1]:
            if count != 2:
                print('Error: No such file. Check orbit number or file.')
                exit()
    return hepd_file_path, mepd_file_path

def file_tree(HEPD_data, MEPD_data):
    ''' This function prints the file structure tree. '''
    print('File tree:')
    HEPD_data.tree(attrs=False)
    print('\n')
    MEPD_data.tree(attrs=False)
    print('\n')


