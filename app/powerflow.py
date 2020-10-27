# Author: Lucas Zenichi Terada
# Institution: University of Campinas

import numpy as np
from numpy.distutils.command.config import config
from app.circuit import *
from data import config


def newton_raphson(Y, K, maxiter, maxerror, pos):
    P = (config.bus.pg - config.bus.pd) / config.param['snom']
    Q = (config.bus.qg - config.bus.qd) / config.param['snom']
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
        H = H_calculation(v, ang, size, G, B, P, Q, Pcalc, Qcalc, pos)
        L = L_calculation(v, ang, size, G, B, P, Q, Pcalc, Qcalc, pos)
        M = M_calculation(v, ang, size, G, B, P, Q, Pcalc, Qcalc, pos)
        N = N_calculation(v, ang, size, G, B, P, Q, Pcalc, Qcalc, pos)
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
        error = mismatche(DP, DQ)
    return v, np.rad2deg(ang), Pcalc, Qcalc, error, iter
