import sympy as sp
import numpy as np
import math


def A_ref_aerodinamico(R_ar, m_dot_3, T3, P3, perda_pressao_total, fator_perda_pressao):
    # perda_pressao_total = deltaP3-4/P3
    # fator_perda_pressao = deltaP3-4/qRef
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
    eta_zr = 0.56 + (0.44*np.tanh(1.5475*0.001*(T3 + (108*np.log(P3)) - 1863)))
    return eta_zr


def eta_zp(T3, P3):
    eta_zp = 0.71 + (0.29*np.tanh(1.5475*0.001*(T3 + (108*np.log(P3)) - 1863)))
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
    eta_zs = 0.71 + (0.29)

    return eta_zs


def tg_zr(T_med_zr):
    tg_zr = T_med_zr
    return tg_zr


def tg_zp(T_out_zp, l_zp, T3, x):
    tg_zp = (((T_out_zp - T3) * x)/l_zp)+T3
    return tg_zp


def tg_zs(T_out_zp, T_out_zs, l_zs, l_zp, x):
    tg_zs = (((T_out_zs - T_out_zp) * x)/(l_zp+l_zs))+T_out_zp
    return tg_zs


def tg_zd(T_out_zs, T_out_zd, l_zd, l_cc, x):
    tg_zd = (((T_out_zd - T_out_zs)*x)/l_cc)+T_out_zs
    return tg_zd


def mg_zr(m_dot_zp):
    mg_zr = (3/4)*m_dot_zp
    return mg_zr


def mg_zp(mg_zr, m_dot_zp, l_zp, l_zr, i):
    mg_zp = mg_zr + ((m_dot_zp - mg_zr)*(i - l_zr)/(l_zp - l_zr))

    return mg_zp


def mg_zs(m_dot_zs, mg_zp, l_zp, l_zs, i):
    mg_zs = mg_zp + (((m_dot_zs - mg_zp)*(i - l_zp))/l_zs)
    return mg_zs


def mg_zd(mg_zs, m_dot_zd, l_zp, l_zs, l_zd, i):
    mg_zd = mg_zs + (((m_dot_zd-mg_zs)*(i-(l_zp+l_zs)))/l_zd)
    return mg_zd
