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


def eta_zr(T3, P3):
    eta_zr = 0.83 + (0.17*np.tanh(1.5475*0.001*(T3 + 108*np.log(P3) - 1863)))
    return eta_zr


def eta_zp(T3, P3):
    eta_zp = 0.92 + (0.12*np.tanh(1.5475*0.001*(T3 + 108*np.log(P3) - 1863)))
    return eta_zp


def psi_T3(T3, P3, phi_zs, m_dot_comb, v_zs):
    psi_T3 = ((10**(-3.054*(phi_zs**-1.205)))*(T3 **
                                               (1.2327*(phi_zs**-1.205))))*(m_dot_comb/(v_zs*P3**(2*phi_zs)))
    return psi_T3


def D_asterisco(P3, deltaP):
    D_asterisco = 0.736 - (0.0173*(P3/deltaP))
    return D_asterisco


def eta_zs(T3, P3, phi_zs, m_dot_comb, v_zs, deltaP):
    eta_zs = 1/(10**10**(0.911*np.log(psi_T3(T3, P3, phi_zs, m_dot_comb, v_zs)) +
                8.02*phi_zs - 1.097 + D_asterisco(P3, deltaP)))
    return eta_zp


def tg_1(T_med_zr):
    Tg1_calculo = T_med_zr
    return Tg1_calculo


def tg_2(T_med_zr, T_saida_zp, L_zp, L_zr, x):
    Tg2_calculo = T_med_zr + \
        ((T_saida_zp - T_med_zr)/(L_zp - L_zr) * (x - L_zr))
    return Tg2_calculo


def tg_3(T_saida_zp, T_saida_zs, L_zs, L_zp, x):
    Tg3_calculo = T_saida_zp + ((T_saida_zs - T_saida_zp)/L_zs)*(x-L_zp)
    return Tg3_calculo


def tg_4(T_saida_zs, T_saida_zd, L_dz, x, L_zp, L_zs):
    Tg3_calculo = T_saida_zs + \
        ((T_saida_zd - T_saida_zs) / L_dz * (x - L_zp - L_zs))
    return Tg3_calculo
