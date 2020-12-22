# constants.py

CONDA_PATH = '/home/seonghohan/anaconda3/share'
DATA_PATH = './isss/data/'

ORBIT_NO = 8795

NORTH_POLE = 90
SOUTH_POLE = -90

GROUP1 = 'HEPD_DIV'
GROUP2 = 'MEPD_SCI'
DATASET1 = 'block1_values'
DATASET2 = 'block2_values'
PDF = 'PDF'
PNG = 'PNG'

# Data indices
HEPD_TIME = 7       # Time (UNIX)
MEPD_TIME = 10      # Time (UNIX)
DT = 4              # Subunit ID
PC1 = 5             # Packet count 1
POS = 16            # Position (deg)
POS_LEN = 3         # Position data length
MAG = 0             # Magnetic field
MAG_LEN = 8         # Magnetic field data length
DET0 = 13           # Detector 0
DET1 = 81           # Detector 1
DET2 = 149          # Detector 2
DET3 = 217          # Detector 3
DET_LEN = 64        # Detector data length
TEL0 = 9            # Telescope 0
TEL1 = 50           # Telescope 1
TEL2 = 91           # Telescope 2
TEL_LEN = 40        # Telescope data length