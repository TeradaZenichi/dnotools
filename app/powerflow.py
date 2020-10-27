# Author: Lucas Zenichi Terada
# Institution: University of Campinas

from app import circuit
from data import config
import numpy as np


def newton_raphson(Y, K, maxiter, maxerror,pos):
    P = config.bus.pg - config.bus.pd
    Q = config.bus.qg - config.bus.qd
    size = len(config.bus)
    Pcalc = np.zeros(size)
    Qcalc = np.zeros(size)
    ang = np.asarray(config.bus.ang)
    v = np.asarray(config.bus.v)
    G = Y.real
    B = Y.imag
    error = np.infty
    H = np.zeros([size, size])
    N = np.zeros([size, size])
    M = np.zeros([size, size])
    L = np.zeros([size, size])
    iter = 0
    for k in range(0, size, 1):
        for m in K[k]:
            m = m-1
            km = ang[k] - ang[m]
            Pcalc[k] = Pcalc[k] + v[k]*v[m]*(G[k, m]*np.cos(km) + B[k, m]*np.sin(km))
            Qcalc[k] = Qcalc[k] + v[k]*v[m]*(G[k, m]*np.sin(km)-B[k, m]*np.cos(km))
    while error > maxerror and iter < maxiter:
        iter = iter + 1
        for k in range(0,size):
            if config.bus.type[k] == 'SL':
                H[k,k] = np.infty
                L[k,k] = np.infty
            elif config.bus.type[k] == 'PV':
                L[k,k] = np.infty
                H[k,k] = -B[k,k]*v[k]**2 - Qcalc[k]
            elif config.bus.type[k] == 'PQ':
                H[k, k] = -B[k, k] * v[k]**2 - Qcalc[k]
                L[k, k] = -Q[k]+(Qcalc[k]-B[k, k]*v[k]**2)/v[k]
            N[k, k] = -P[k]+(Pcalc[k]+G[k, k]*v[k]** 2)/v[k]
            M[k, k] = -G[k, k]*v[k]**2 + Pcalc[k]
        for line in config.branches.itertuples():
            k = line.start-1
            m = line.end-1
            km = ang[k] - ang[m]
            H[k, m] = v[k] * v[m] * (G[k, m] * np.sin(km) - B[k, m] * np.cos(km))
            H[m, k] = -v[k] * v[m] * (G[k, m] * np.sin(km) + B[k, m] * np.cos(km))
            N[k, m] = v[k] * (G[k, m] * np.cos(km) + B[k, m] * np.sin(km))
            N[m, k] = v[m] * (G[k, m] * np.cos(km) - B[k, m] * np.sin(km))
            M[k, m] = -v[k] * v[m] * (G[k, m] * np.cos(km) + B[k, m] * np.sin(km))
            M[m, k] = -v[k] * v[m] * (G[k, m] * np.cos(km) - B[k, m] * np.sin(km))
            L[k, m] = v[k] * (G[k, m] * np.sin(km) - B[k, m] * np.cos(km))
            L[m, k] = -v[m] * (G[k, m] * np.sin(km) + B[k, m] * np.cos(km))
        if config.switches is not None:
            for line in config.switches.itertuples():
                if pos[line.Index]:
                    k = line.start - 1
                    m = line.end - 1
                    km = ang[k] - ang[m]
                    H[k, m] = v[k] * v[m] * (G[k, m] * np.sin(km) - B[k, m] * np.cos(km))
                    H[m, k] = -v[k] * v[m] * (G[k, m] * np.sin(km) + B[k, m] * np.cos(km))
                    N[k, m] = v[k] * (G[k, m] * np.cos(km) + B[k, m] * np.sin(km))
                    N[m, k] = v[m] * (G[k, m] * np.cos(km) - B[k, m] * np.sin(km))
                    M[k, m] = -v[k] * v[m] * (G[k, m] * np.cos(km) + B[k, m] * np.sin(km))
                    M[m, k] = -v[k] * v[m] * (G[k, m] * np.cos(km) - B[k, m] * np.sin(km))
                    L[k, m] = v[k] * (G[k, m] * np.sin(km) - B[k, m] * np.cos(km))
                    L[m, k] = -v[m] * (G[k, m] * np.sin(km) + B[k, m] * np.cos(km))
        J1 = np.hstack((H,N))
        J2 = np.hstack((M,L))
        J = np.vstack((J1,J2))
        if np.linalg.det(J) == 0:
            error = np.infty
            break
        DP = P - Pcalc
        DQ = Q - Qcalc
        DS = np.concatenate((DP,DQ))
        delta = np.linalg.solve(J,DS)
        print(delta)
    return v, ang, Pcalc, Qcalc,error,J,DP,DQ
