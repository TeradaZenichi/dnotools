# Author: Lucas Zenichi Terada
# Institution: University of Campinas

import numpy as np
from numpy.distutils.command.config import config

from data import config


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


def L_calculation(v, ang, size, G, B, P, Q, Pcalc, Qcalc,pos):
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


def M_calculation(v, ang, size, G, B, P, Q, Pcalc, Qcalc,pos):
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


def N_calculation(v, ang, size, G, B, P, Q, Pcalc, Qcalc,pos):
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

def mismatche(DP,DQ):
    mism = 0
    for bar in config.bus.itertuples():
        if bar.type == 'PV':
            mism = np.maximum(mism,np.abs(DP[bar.Index]))
        elif bar.type == 'PQ':
            mism = np.maximum(mism,np.abs(DP[bar.Index]))
            mism = np.maximum(mism,np.abs(DQ[bar.Index]))
    return mism

def newton_raphson(Y, K, maxiter, maxerror, pos):
    P = (config.bus.pg - config.bus.pd)/config.param['snom']
    Q = (config.bus.qg - config.bus.qd)/config.param['snom']
    size = len(config.bus)
    ang = np.asarray(config.bus.ang)
    v = np.asarray(config.bus.v)
    G = Y.real
    B = Y.imag
    error = np.infty
    iter = 0
    (Pcalc, Qcalc) = power_calculation(v, ang, size, K, G, B)
    while error > maxerror and iter < maxiter:
        iter = iter + 1
        H = H_calculation(v, ang, size, G, B, P, Q, Pcalc, Qcalc,pos)
        L = L_calculation(v, ang, size, G, B, P, Q, Pcalc, Qcalc,pos)
        M = M_calculation(v, ang, size, G, B, P, Q, Pcalc, Qcalc,pos)
        N = N_calculation(v, ang, size, G, B, P, Q, Pcalc, Qcalc,pos)
        J1 = np.hstack((H, N))
        J2 = np.hstack((M, L))
        J = np.vstack((J1, J2))
        if np.linalg.det(J) == 0:
            error = np.infty
            break
        DP = P - Pcalc
        DQ = Q - Qcalc
        DS = np.concatenate((DP, DQ))
        delta = np.linalg.solve(J, DS)
        ang = ang + delta[0:size]
        v = v + delta[size:2 * size]
        (Pcalc, Qcalc) = power_calculation(v, ang, size, K, G, B)
        DP = P - Pcalc
        DQ = Q - Qcalc
        error = mismatche(DP,DQ)
    return v, ang, Pcalc, Qcalc, error, iter
