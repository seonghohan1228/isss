# file_functions.py

import os

def show_avail_file():
    """ This function shows available files in the data directory. """
    files = os.listdir(os.getcwd() + '../data/')
    print('\nAvailable files:')
    for filename in files:
        print(filename)
    print('\n')
    