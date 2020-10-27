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

def power_calculation(v, ang, size, K, G, B):
    Pcalc = np.zeros(size)
    Qcalc = np.zeros(size)
    for k in range(0, size, 1):
        for m in K[k]:
            m = m - 1
            km = ang[k] - ang[m]
            Pcalc[k] = Pcalc[k] + v[k] * v[m] * (G[k, m] * np.cos(km) + B[k, m] * np.sin(km))
            Qcalc[k] = Qcalc[k] + v[k] * v[m] * (G[k, m] * np.sin(km) - B[k, m] * np.cos(km))
    return Pcalc, Qcalc


def H_calculation(v, ang, size, G, B, P, Q, Pcalc, Qcalc, pos):
    H = np.zeros([size, size])
    for k in range(0, size):
        if config.bus.type[k] == 'SL':
            H[k, k] = np.infty
        elif config.bus.type[k] == 'PV':
            H[k, k] = -B[k, k] * v[k] ** 2 - Qcalc[k]
        elif config.bus.type[k] == 'PQ':
            H[k, k] = -B[k, k] * v[k] ** 2 - Qcalc[k]
    for line in config.branches.itertuples():
        k = line.start - 1
        m = line.end - 1
        km = ang[k] - ang[m]
        H[k, m] = v[k] * v[m] * (G[k, m] * np.sin(km) - B[k, m] * np.cos(km))
        H[m, k] = -v[k] * v[m] * (G[k, m] * np.sin(km) + B[k, m] * np.cos(km))
    if config.switches is not None:
        for line in config.switches.itertuples():
            if pos[line.Index]:
                k = line.start - 1
                m = line.end - 1
                km = ang[k] - ang[m]
                H[k, m] = v[k] * v[m] * (G[k, m] * np.sin(km) - B[k, m] * np.cos(km))
                H[m, k] = -v[k] * v[m] * (G[k, m] * np.sin(km) + B[k, m] * np.cos(km))
    return H


def L_calculation(v, ang, size, G, B, P, Q, Pcalc, Qcalc, pos):
    L = np.zeros([size, size])
    for k in range(0, size):
        if config.bus.type[k] == 'SL':
            L[k, k] = np.infty
        elif config.bus.type[k] == 'PV':
            L[k, k] = np.infty
        elif config.bus.type[k] == 'PQ':
            L[k, k] = -Q[k] + (Qcalc[k] - B[k, k] * v[k] ** 2) / v[k]
    for line in config.branches.itertuples():
        k = line.start - 1
        m = line.end - 1
        km = ang[k] - ang[m]
        L[k, m] = v[k] * (G[k, m] * np.sin(km) - B[k, m] * np.cos(km))
        L[m, k] = -v[m] * (G[k, m] * np.sin(km) + B[k, m] * np.cos(km))
    if config.switches is not None:
        for line in config.switches.itertuples():
            if pos[line.Index]:
                k = line.start - 1
                m = line.end - 1
                km = ang[k] - ang[m]
                L[k, m] = v[k] * (G[k, m] * np.sin(km) - B[k, m] * np.cos(km))
                L[m, k] = -v[m] * (G[k, m] * np.sin(km) + B[k, m] * np.cos(km))
    return L


def M_calculation(v, ang, size, G, B, P, Q, Pcalc, Qcalc, pos):
    M = np.zeros([size, size])
    for k in range(0, size):
        M[k, k] = -G[k, k] * v[k] ** 2 + Pcalc[k]
    for line in config.branches.itertuples():
        k = line.start - 1
        m = line.end - 1
        km = ang[k] - ang[m]
        M[k, m] = -v[k] * v[m] * (G[k, m] * np.cos(km) + B[k, m] * np.sin(km))
        M[m, k] = -v[k] * v[m] * (G[k, m] * np.cos(km) - B[k, m] * np.sin(km))
    if config.switches is not None:
        for line in config.switches.itertuples():
            if pos[line.Index]:
                k = line.start - 1
                m = line.end - 1
                km = ang[k] - ang[m]
                M[k, m] = -v[k] * v[m] * (G[k, m] * np.cos(km) + B[k, m] * np.sin(km))
                M[m, k] = -v[k] * v[m] * (G[k, m] * np.cos(km) - B[k, m] * np.sin(km))
    return M


def N_calculation(v, ang, size, G, B, P, Q, Pcalc, Qcalc, pos):
    N = np.zeros([size, size])
    for k in range(0, size):
        N[k, k] = -P[k] + (Pcalc[k] + G[k, k] * v[k] ** 2) / v[k]
    for line in config.branches.itertuples():
        k = line.start - 1
        m = line.end - 1
        km = ang[k] - ang[m]
        N[k, m] = v[k] * (G[k, m] * np.cos(km) + B[k, m] * np.sin(km))
        N[m, k] = v[m] * (G[k, m] * np.cos(km) - B[k, m] * np.sin(km))
    if config.switches is not None:
        for line in config.switches.itertuples():
            if pos[line.Index]:
                k = line.start - 1
                m = line.end - 1
                km = ang[k] - ang[m]
                N[k, m] = v[k] * (G[k, m] * np.cos(km) + B[k, m] * np.sin(km))
                N[m, k] = v[m] * (G[k, m] * np.cos(km) - B[k, m] * np.sin(km))
    return N


def mismatche(DP, DQ):
    mism = 0
    for bar in config.bus.itertuples():
        if bar.type == 'PV':
            mism = np.maximum(mism, np.abs(DP[bar.Index]))
        elif bar.type == 'PQ':
            mism = np.maximum(mism, np.abs(DP[bar.Index]))
            mism = np.maximum(mism, np.abs(DQ[bar.Index]))
    return mism


