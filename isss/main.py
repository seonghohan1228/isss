# main.py
from functions import env_setup, data, plot
orbit_no = 8795

def main():
    env_setup()
    HEPD, MEPD = data(orbit_no, True, True)
    plot(HEPD, MEPD, 'south', 'spacepy', 'plot_combined', 'pdf')

if __name__ == '__main__':
    main()

