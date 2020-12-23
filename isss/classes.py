# classes.py
import numpy as np
from data_constants import *

class isssData():
    def __init__(self, orbit_no, data, group, ds1, ds2):
        self.orbit_no = orbit_no
        self.dataset1 = data[group][ds1]
        self.dataset2 = data[group][ds2]
        self.pc1 = np.array(self.dataset1[:, PC1])
        self.pos = np.array(self.dataset2[:, POS:POS + POS_LEN])
        self.mag = np.array(self.dataset2[:, MAG:MAG + MAG_LEN])
        if group == 'HEPD_DIV':
            self.time = np.array(self.dataset1[:, HEPD_TIME])
            tel0 = np.array(self.dataset1[:, TEL0:TEL0 + TEL_LEN])
            tel1 = np.array(self.dataset1[:, TEL1:TEL1 + TEL_LEN])
            tel2 = np.array(self.dataset1[:, TEL2:TEL2 + TEL_LEN])
            self.tel = np.array([tel0, tel1, tel2])
        else:
            self.time = np.array(self.dataset1[:, MEPD_TIME])
            self.dt = np.array(self.dataset1[:, DT])
            det0 = np.array(self.dataset1[:, DET0:DET0 + DET_LEN])
            det1 = np.array(self.dataset1[:, DET1:DET1 + DET_LEN])
            det2 = np.array(self.dataset1[:, DET2:DET2 + DET_LEN])
            det3 = np.array(self.dataset1[:, DET3:DET3 + DET_LEN])
            self.det = np.array([det0, det1, det2, det3])
    
