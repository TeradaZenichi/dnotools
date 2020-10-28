# Author: Lucas Zenichi Terada
# Institution: University of Campinas

import numpy as np
from numpy.distutils.command.config import config
from app.circuit import *
from app.tools import *
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


def declouped(Y, K, maxiter, maxerror, pos):
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
        if np.linalg.det(H) == 0 or np.linalg.det(L) == 0:
            error = np.infty
            break
        (Pcalc, Qcalc) = power_calculation(v, ang, size, K, G, B)
        DP = P - Pcalc
        DQ = Q - Qcalc
        ang = ang + np.linalg.solve(H, DP)
        v = v + np.linalg.solve(L, DQ)
        (Pcalc, Qcalc) = power_calculation(v, ang, size, K, G, B)
        DP = P - Pcalc
        DQ = Q - Qcalc
        error = mismatche(DP, DQ)
    return v, np.rad2deg(ang), Pcalc, Qcalc, error, iter


def alternate_declouped(Y, K, maxiter, maxerror, pos):
    P = (config.bus.pg - config.bus.pd) / config.param['snom']
    Q = (config.bus.qg - config.bus.qd) / config.param['snom']
    size = len(config.bus)
    ang = np.asarray(config.bus.ang)
    v = np.asarray(config.bus.v)
    G = Y.real
    B = Y.imag
    error = np.infty
    iter = 0
    (KP,KQ,p,q) = (1,1,0,0)
    (Pcalc, Qcalc) = power_calculation(v, ang, size, K, G, B)
    while iter < maxiter:
        (Pcalc) = active_power_calculation(v, ang, size, K, G, B)
        DP = P - Pcalc
        if active_mismatche(DP) < maxerror:
            KP = 0
            if KQ == 0:
                break
        if active_mismatche(DP) >= maxerror:
            H = H_calculation(v, ang, size, G, B, P, Q, Pcalc, Qcalc, pos)
            ang = ang + np.linalg.solve(H,DP)
            p = p+1
            KQ = 1
        Qcalc = reactive_power_calculation(v, ang, size, K, G, B)
        DQ = Q - Qcalc
        if reactive_mismatche(DQ) < maxerror:
            KQ = 0
            if KP == 0:
                break
        if reactive_mismatche(DQ) >= maxerror:
            L = L_calculation(v, ang, size, G, B, P, Q, Pcalc, Qcalc, pos)
            v = v + np.linalg.solve(L,DQ)
            q = q+1
            KP = 1
        iter = np.minimum(p,q)
    return v, np.rad2deg(ang), Pcalc, Qcalc, error, iter

def fast_declouped(Y, K, maxiter, maxerror, pos):
    P = (config.bus.pg - config.bus.pd) / config.param['snom']
    Q = (config.bus.qg - config.bus.qd) / config.param['snom']
    size = len(config.bus)
    ang = np.asarray(config.bus.ang)
    v = np.asarray(config.bus.v)
    G = Y.real
    B = Y.imag
    error = np.infty
    iter = 0
    (KP, KQ, p, q) = (1, 1, 0, 0)
    (B1,B2) = B_calculation(size, pos, B)
    while iter < maxiter:
        (Pcalc) = active_power_calculation(v, ang, size, K, G, B)
        DP = P - Pcalc
        if active_mismatche(DP) < maxerror:
            KP = 0
            if KQ == 0:
                break
        if active_mismatche(DP) >= maxerror:
            ang = ang + np.linalg.solve(B1,DP)
            p = p+1
            KQ = 1
        Qcalc = reactive_power_calculation(v, ang, size, K, G, B)
        DQ = Q - Qcalc
        if reactive_mismatche(DQ) < maxerror:
            KQ = 0
            if KP == 0:
                break
        if reactive_mismatche(DQ) >= maxerror:
            v = v + np.linalg.solve(B2,DQ)
            q = q+1
            KP = 1
        iter = np.minimum(p,q)
    return v, np.rad2deg(ang), Pcalc, Qcalc, error, iter