# Author: Lucas Zenichi Terada
# Institution: University of Campinas
# Description: This program use the config global variables

from data import config
import pandas as pd
import numpy as np

def lines_admittance():
    config.loadfiles(config)  #Retirar
    size = config.bus.__len__()
    Ynet = np.zeros([size,size],dtype=complex)
    for bar in config.bus.itertuples():
        Ynet[bar.num-1,bar.num-1] = bar.Gsh+1j*bar.Bsh
    for line in config.branches.itertuples():
        k = line.start - 1
        m = line.end - 1
        Ynet[k, m] = -(1/(line.r + 1j*line.x))/line.tap
        Ynet[m, k] = -(1/(line.r + 1j*line.x))/line.tap
        Ynet[k, k] = Ynet[k, k] + 1j*line.bsh/2 + (1/(line.r + 1j*line.x))/(line.tap**2)
        Ynet[m, m] = Ynet[m, m] + 1/(line.r + 1j*line.x) + 1j*line.bsh/2
    return Ynet

def switches_addmitance(pos):
    config.loadfiles(config)  #Retirar
    size = config.bus.__len__()
    Ysw = np.zeros([size, size], dtype=complex)
    for line in config.switches.itertuples():
        if pos[line.Index]:
            k = line.start - 1
            m = line.end - 1
            Ysw[k, m] = -(1/(line.r + 1j * line.x))/line.tap
            Ysw[m, k] = -(1/(line.r + 1j * line.x))/line.tap
            Ysw[k, k] = Ysw[k, k] + 1j*line.bsh/2 + (1/(line.r + 1j*line.x))/(line.tap**2)
            Ysw[m, m] = Ysw[m, m] + 1/(line.r + 1j * line.x) + 1j*line.bsh/2
    return Ysw

def near_bars():
    K = []
    for bar in config.bus.itertuples():
        Kaux = [bar.num]
        for line in config.branches.itertuples():
            if bar.num == line.start:
                Kaux.append(line.end)
            if bar.num == line.end:
                Kaux.append(line.start)
        K.append(Kaux)
    return K

def near_switches(pos):
    K = []
    for bar in config.bus.itertuples():
        Kaux = []
        for line in config.switches.itertuples():
            if bar.num == line.start and pos[line.Index] == 1:
                Kaux.append(line.end)
            if bar.num == line.end and pos[line.Index] == 1:
                Kaux.append(line.start)
        K.append(Kaux)
    return K


