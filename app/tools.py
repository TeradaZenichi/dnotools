# Author: Lucas Zenichi Terada
# Insitution: Univeristy of Campinas

from data import config
import numpy as np


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


def active_power_calculation(v, ang, size, K, G, B):
    Pcalc = np.zeros(size)
    for k in range(0, size, 1):
        for m in K[k]:
            m = m - 1
            km = ang[k] - ang[m]
            Pcalc[k] = Pcalc[k] + v[k] * v[m] * (G[k, m] * np.cos(km) + B[k, m] * np.sin(km))
    return Pcalc


def reactive_power_calculation(v, ang, size, K, G, B):
    Qcalc = np.zeros(size)
    for k in range(0, size, 1):
        for m in K[k]:
            m = m - 1
            km = ang[k] - ang[m]
            Qcalc[k] = Qcalc[k] + v[k] * v[m] * (G[k, m] * np.sin(km) - B[k, m] * np.cos(km))
    return Qcalc


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


def B_calculation(size, pos, B):
    B1 = np.zeros([size, size])
    B2 = np.zeros([size, size])
    for k in range(0, size):
        if config.bus.type[k] == 'SL':
            B1[k, k] = np.infty
            B2[k, k] = np.infty
        elif config.bus.type[k] == 'PV':
            B2[k, k] = np.infty
            aux = 0;
            if config.branches is not None:
                for line in config.branches.itertuples():
                    if config.bus.num[k] == line.start:
                        aux = aux + 1/line.x
                    if config.bus.num[k] == line.end:
                        aux = aux+1/line.x
            if config.switches is not None:
                for line in config.switches.itertuples():
                    if pos[line.Index]:
                        if config.bus.num[k] == line.start:
                            aux = aux + 1/line.x
                        if config.bus.num[k] == line.end:
                            aux = aux+1/line.x
            B1[k, k] = aux
        elif config.bus.type[k] == 'PQ':
            aux = 0;
            if config.branches is not None:
                for line in config.branches.itertuples():
                    if config.bus.num[k] == line.start:
                        aux = aux + 1 / line.x
                    if config.bus.num[k] == line.end:
                        aux = aux + 1 / line.x
            if config.switches is not None:
                for line in config.switches.itertuples():
                    if pos[line.Index]:
                        if config.bus.num[k] == line.start:
                            aux = aux + 1 / line.x
                        if config.bus.num[k] == line.end:
                            aux = aux + 1 / line.x
            B1[k, k] = aux
            B2[k, k] = -B[k, k]
    for line in config.branches.itertuples():
        k = line.start - 1
        m = line.end - 1
        B1[k, m] = -1/line.x
        B1[m, k] = -1/line.x
        B2[k, m] = -B[k, m]
        B2[m, k] = -B[k, m]
    if config.switches is not None:
        for line in config.switches.itertuples():
            if pos[line.Index]:
                k = line.start - 1
                m = line.end - 1
                B1[k, m] = -1/line.x
                B1[m, k] = -1/line.x
                B2[k, m] = -B[k, m]
                B2[m, k] = -B[k, m]
    return B1, B2


def mismatche(DP, DQ):
    mism = 0
    for bar in config.bus.itertuples():
        if bar.type == 'PV':
            mism = np.maximum(mism, np.abs(DP[bar.Index]))
        elif bar.type == 'PQ':
            mism = np.maximum(mism, np.abs(DP[bar.Index]))
            mism = np.maximum(mism, np.abs(DQ[bar.Index]))
    return mism


def active_mismatche(DP):
    mism = 0
    for bar in config.bus.itertuples():
        if bar.type == 'PV':
            mism = np.maximum(mism, np.abs(DP[bar.Index]))
        elif bar.type == 'PQ':
            mism = np.maximum(mism, np.abs(DP[bar.Index]))
    return mism


def reactive_mismatche(DQ):
    mism = 0
    for bar in config.bus.itertuples():
        if bar.type == 'PQ':
            mism = np.maximum(mism, np.abs(DQ[bar.Index]))
    return mism
