import sympy as sp
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import scipy as spy


def A_ref_aerodinamico(R_ar, m_dot_3, T3, P3, perda_pressao_total, fator_perda_pressao):
    #perda_pressao_total = deltaP3-4/P3
    #fator_perda_pressao = deltaP3-4/qRef
    A_ref_aero = math.sqrt(R_ar * (math.pow((m_dot_3 * math.sqrt(T3)) / P3, 2))
                           * (fator_perda_pressao / perda_pressao_total))
    return A_ref_aero


def A_ref_quimico(m_dot_3, T3, P3, theta, D_ref, b):
    A_ref_q = (theta * m_dot_3) / (math.pow(P3, 1.75) *
                                   math.pow(D_ref, 0.75) * math.exp(T3 / b))
    return A_ref_q


def D_ref(A_ref, D_int):
    D_ref = (math.sqrt(((4 * A_ref)/math.pi) + math.pow(D_int, 2)) - D_int) / 2
    return D_ref


def A_ft(A_ref):
    A_ft = 0.7 * A_ref
    return A_ft


def D_ft(A_ft, D_int, D_ref):
    D_ft = A_ft / (math.pi * (D_int + D_ref))
    return D_ft


def phi_global(m_dot_comb, m_dot_3, phi_estq):
    phi_global = (m_dot_comb / m_dot_3) / phi_estq
    return phi_global


def phi_pobre(T3):
    phi_pobre = 0.67 - (0.0004*T3)
    return phi_pobre


def phi_rico(T3):
    phi_rico = 1.82 + (0.0006*T3)
    return phi_rico


def phi_zp(phi_global, air_zp_per_cent):
    phi_zp = phi_global / air_zp_per_cent
    return phi_zp
