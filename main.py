import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline, BSpline
import sympy as sp

import equations as eq
import dataFramesPrinter as printer


P3_EM = 2000000  # MPA
T3_EM = 814  # K
m3_EM = 18.1  # kg/s
T4_EM = 1600  # K
mComb_EM = 0.427  # kg/s

P3_MA = 700000  # MPA
T3_MA = 707  # K
m3_MA = 6.8  # kg/s
T4_MA = 1387  # K
mComb_MA = 0.132  # kg/s

P3_CRU = 1800000  # MPA
T3_CRU = 1060  # K
m3_CRU = 14.2  # kg/s
T4_CRU = 1393  # K
mComb_CRU = 0.140  # kg/s

P3_IDLE = 150000  # MPA
T3_IDLE = 343  # K
m3_IDLE = 1.05  # kg/s
T4_IDLE = 703  # K
mComb_IDLE = 0.0091  # kg/s

R_ar = 143.5  # J/kg*K
theta = 73000000
D_int = 0.05
phi_zp = 1.05  # deprecated
phi_zs = 8
phi_estq = 0.06818
b = 170*(2-math.log(phi_zp))
m_dot_zp = 0.245

perda_pressao_total = 0.06  # deltaP3-4/P3
fator_perda_pressao = 20  # deltaP3-4/qRef

# Razões estequimetricas

phi_global_EM = eq.phi_global(mComb_EM, m3_EM, phi_estq)
phi_global_MA = eq.phi_global(mComb_MA, m3_MA, phi_estq)
phi_global_CRU = eq.phi_global(mComb_CRU, m3_CRU, phi_estq)
phi_global_IDLE = eq.phi_global(mComb_IDLE, m3_IDLE, phi_estq)

phi_rico_EM = eq.phi_rico(T3_EM)
phi_rico_MA = eq.phi_rico(T3_MA)
phi_rico_CRU = eq.phi_rico(T3_CRU)
phi_rico_IDLE = eq.phi_rico(T3_IDLE)

phi_pobre_EM = eq.phi_pobre(T3_EM)
phi_pobre_MA = eq.phi_pobre(T3_MA)
phi_pobre_CRU = eq.phi_pobre(T3_CRU)
phi_pobre_IDLE = eq.phi_pobre(T3_IDLE)

printer.print_phis(phi_pobre_EM, phi_pobre_MA, phi_pobre_CRU, phi_pobre_IDLE, phi_rico_EM, phi_rico_MA,
                   phi_rico_CRU, phi_rico_IDLE, phi_global_EM, phi_global_MA, phi_global_CRU,  phi_global_IDLE)


# Limites de equivalencias pobre/rico
limit_equivalence_pobre_EM = phi_global_EM/phi_rico_EM
limit_equivalence_pobre_MA = phi_global_MA/phi_rico_MA
limit_equivalence_pobre_CRU = phi_global_CRU/phi_rico_CRU
limit_equivalence_pobre_IDLE = phi_global_IDLE/phi_rico_IDLE

limit_equivalence_rico_EM = phi_global_EM/phi_pobre_EM
limit_equivalence_rico_MA = phi_global_MA/phi_pobre_MA
limit_equivalence_rico_CRU = phi_global_CRU/phi_pobre_CRU
limit_equivalence_rico_IDLE = phi_global_IDLE/phi_pobre_IDLE

# printer. print_limitEquivalence(limit_equivalence_pobre_EM, limit_equivalence_pobre_MA, limit_equivalence_pobre_CRU, limit_equivalence_pobre_IDLE,
#                                 limit_equivalence_rico_EM, limit_equivalence_rico_MA, limit_equivalence_rico_CRU, limit_equivalence_rico_IDLE)

# Definição da porcentagem de ar na zona primeira a partir dos criterios 1-3 pag 53
air_zp_max = phi_global_EM/1.05  # 0.266
air_zp_min = phi_global_EM/1.5  # 0.1898

# Valor minimo definido a parir limite de equivalencia rico (IDLE)
air_zp_percent = 0.238578

phi_zp_EM = eq.phi_zp(phi_global_EM, air_zp_percent)
phi_zp_MA = eq.phi_zp(phi_global_MA, air_zp_percent)
phi_zp_CRU = eq.phi_zp(phi_global_CRU, air_zp_percent)
phi_zp_IDLE = eq.phi_zp(phi_global_IDLE, air_zp_percent)

# printer.print_phisZP(phi_zp_EM, phi_zp_MA, phi_zp_CRU, phi_zp_IDLE)

# Determinado o maior phi_zp das iterações
phi_zp = phi_zp_EM  # 1.45


# Determinação de parametros de operação

A_ref_aero_EM = eq.A_ref_aerodinamico(
    R_ar, m3_EM, T3_EM, P3_EM, perda_pressao_total, fator_perda_pressao)
A_ref_aero_MA = eq.A_ref_aerodinamico(
    R_ar, m3_MA, T3_MA, P3_MA, perda_pressao_total, fator_perda_pressao)
A_ref_aero_CRU = eq.A_ref_aerodinamico(
    R_ar, m3_CRU, T3_CRU, P3_CRU, perda_pressao_total, fator_perda_pressao)
A_ref_aero_IDLE = eq.A_ref_aerodinamico(
    R_ar, m3_IDLE, T3_IDLE, P3_IDLE, perda_pressao_total, fator_perda_pressao)

D_ref_aero_EM = eq.D_ref(A_ref_aero_EM, D_int)
D_ref_aero_MA = eq.D_ref(A_ref_aero_MA, D_int)
D_ref_aero_CRU = eq.D_ref(A_ref_aero_CRU, D_int)
D_ref_aero_IDLE = eq.D_ref(A_ref_aero_IDLE, D_int)

A_ref_quim_EM = eq.A_ref_quimico(
    m3_EM, T3_EM, P3_EM, theta, D_ref_aero_EM, b)
A_ref_quim_MA = eq.A_ref_quimico(
    m3_MA, T3_MA, P3_MA, theta, D_ref_aero_MA, b)
A_ref_quim_CRU = eq.A_ref_quimico(
    m3_CRU, T3_CRU, P3_CRU, theta, D_ref_aero_CRU, b)
A_ref_quim_IDLE = eq.A_ref_quimico(
    m3_IDLE, T3_IDLE, P3_IDLE, theta, D_ref_aero_IDLE, b)

D_ref_quim_EM = eq.D_ref(A_ref_quim_EM, D_int)
D_ref_quim_MA = eq.D_ref(A_ref_quim_MA, D_int)
D_ref_quim_CRU = eq.D_ref(A_ref_quim_CRU, D_int)
D_ref_quim_IDLE = eq.D_ref(A_ref_quim_IDLE, D_int)

A_ft_aero_EM = eq.A_ft(A_ref_aero_EM)
A_ft_aero_MA = eq.A_ft(A_ref_aero_MA)
A_ft_aero_CRU = eq.A_ft(A_ref_aero_CRU)
A_ft_aero_IDLE = eq.A_ft(A_ref_aero_IDLE)

A_ft_quim_EM = eq.A_ft(A_ref_quim_EM)
A_ft_quim_MA = eq.A_ft(A_ref_quim_MA)
A_ft_quim_CRU = eq.A_ft(A_ref_quim_CRU)
A_ft_quim_IDLE = eq.A_ft(A_ref_quim_IDLE)

D_ft_aero_EM = eq.D_ft(A_ft_aero_EM, D_int, D_ref_aero_EM)
D_ft_aero_MA = eq.D_ft(A_ft_aero_MA, D_int, D_ref_aero_MA)
D_ft_aero_CRU = eq.D_ft(A_ft_aero_CRU, D_int, D_ref_aero_CRU)
D_ft_aero_IDLE = eq.D_ft(A_ft_aero_IDLE, D_int, D_ref_aero_IDLE)

D_ft_quim_EM = eq.D_ft(A_ft_quim_EM, D_int, D_ref_quim_EM)
D_ft_quim_MA = eq.D_ft(A_ft_quim_MA, D_int, D_ref_quim_MA)
D_ft_quim_CRU = eq.D_ft(A_ft_quim_CRU, D_int, D_ref_quim_CRU)
D_ft_quim_IDLE = eq.D_ft(A_ft_quim_IDLE, D_int, D_ref_quim_IDLE)

# printer. print_areas(A_ref_aero_EM, A_ref_aero_MA, A_ref_aero_CRU, A_ref_aero_IDLE, A_ref_quim_EM, A_ref_quim_MA, A_ref_quim_CRU, A_ref_quim_IDLE,
#                      D_ref_aero_EM, D_ref_aero_MA, D_ref_aero_CRU, D_ref_aero_IDLE,
#                      D_ref_quim_EM, D_ref_quim_MA, D_ref_quim_CRU, D_ref_quim_IDLE,
#                      A_ft_aero_EM, A_ft_aero_MA, A_ft_aero_CRU, A_ft_aero_IDLE,
#                      A_ft_quim_EM, A_ft_quim_MA, A_ft_quim_CRU, A_ft_quim_IDLE,
#                      D_ft_aero_EM, D_ft_aero_MA, D_ft_aero_CRU, D_ft_aero_IDLE,
#                      D_ft_quim_EM, D_ft_quim_MA, D_ft_quim_CRU, D_ft_quim_IDLE)

A_ref_maior = 0.1693
D_ref_maior = 0.2085
A_ft_maior = 0.1185
D_ft_maior = 0.1460

# print(f'A_ref: {A_ref_maior:.4f} \nD_ref: {D_ref_maior:.4f}\nA_ft:  {A_ft_maior:.4f}\nD_ft:  {D_ft_maior:.4f}')

# Determinação dos comprimentos da camara

l_zr = (1 / 2) * D_ft_maior
l_zp = (3 / 4) * D_ft_maior
l_zs = (1 / 2) * D_ft_maior
l_zd = (3 / 2) * D_ft_maior
l_cc = l_zp + l_zs + l_zd

# printer.print_comprimentos(l_zp, l_zs, l_zd, l_cc)


# Continuando apenas com os dados do ponto de projeto MaximoEmpuxo

# Porcentagem de Ar apra resfriamento deve dar em torno de 50%
# 51.4% | m_dot_Areef/m_dot_3
air_arrefecimento_percent = (((0.1*T3_EM) - 30)/100)
m_dot_arref = air_arrefecimento_percent*m3_EM

# Fluxo de ar na zona primaria
phi_zp
air_zp_percent
m_dot_zp = air_zp_percent*m3_EM


# Fluxo de ar na zona secundaria
phi_zs = 0.8  # definido pela literatura
air_zs_percent = (phi_global_EM/phi_zs) - air_zp_percent
m_dot_zs = air_zs_percent*m3_EM


# Fluxo de massa na zona de diluição

air_zd_percent = 1 - (air_zp_percent + air_zs_percent +
                      air_arrefecimento_percent)

phi_zd = phi_global_EM / air_zd_percent

m_dot_zd = air_zd_percent*m3_EM

massasDF = pd.DataFrame({'m_dot': [m_dot_zp, m_dot_zs, m_dot_zd, m_dot_arref],
                         })
massasDF.index = ['Zona Primaria', "Zona Secundaria",
                  "Zona Diluicao", "Arrefecimento"]
print('\n', massasDF, '\n')

bal_massa = m3_EM - (m_dot_zp + m_dot_zs + m_dot_zd + m_dot_arref)

# ZONA RECIRCULACAO
T_in_zr = T3_EM  # Temperatura de Entrada na ZR
n_zr = eq.eta_zr(T3_EM, P3_EM)
# Valor de 1692 T encontrado através do gráfico gerado pelo cantera, referente a interssecao com phi=1
delta_Tphi = 1692
T_out_zr = T3_EM + n_zr*delta_Tphi  # Temperatura de Saída na ZR
T_med_zr = (T_in_zr/3) + (2*T_out_zr/3)
print("Temp saida ZR", T_out_zr)
print("Temp media ZR", T_med_zr)

# ZONA PRIMARIA
T_in_zp = T_med_zr
# Valor de 1570.5 T encontrado através do gráfico gerado pelo cantera, referente a interssecao com phi_zp=1.45
delta_T_zp = 1570.5
n_zp = eq.eta_zp(T3_EM, P3_EM)
T_out_zp = T3_EM + n_zp*delta_T_zp
print("Temp saida ZP", T_out_zp)

#  ZONA SECUNDARIA
Tin_zs = T_out_zp
n_zs = 1/phi_zs
v_zs = A_ft_maior*l_zs
# Valor de 1557.98 T encontrado através do gráfico gerado pelo cantera, referente a interssecao com phi_zs
delta_T_zs = 1557.98
deltaP = P3_EM*0.06  # Retirado da igualdade deltaP34/P3 = 6%
eta_zs = eq.eta_zs(T3_EM, P3_EM, phi_zs, mComb_EM, v_zs, deltaP)
T_out_zs = T3_EM + n_zs*delta_T_zs
print("Temp saida ZS", T_out_zs)

#  ZONA DILUICAO
v_zd = A_ft_maior*l_zd
# Valor de  812.7 T encontrado através do gráfico gerado pelo cantera, referente a interssecao com phi_global
delta_T_zd = 812.7
deltaP = P3_EM*0.06  # Retirado da igualdade deltaP34/P3 = 6%
eta_zd = eq.eta_zs(T3_EM, P3_EM, phi_zd, mComb_EM, v_zd, deltaP)
T_out_zd = T4_EM
print("Temp saida ZD", T_out_zd)

x_t = np.linspace(1e-4, l_zp+l_zs+l_zd, 100)
Tg = np.zeros(len(x_t))
index = 0

for i in x_t:
    if i >= 0 and i <= l_zr:
        Tg_1 = T_med_zr
        Tg[index] = Tg_1
    if i > l_zr and i <= l_zp:
        Tg_2 = T_med_zr + ((T_out_zp-T_med_zr)/(l_zp-l_zr))*(i-l_zr)
        Tg[index] = Tg_2
    if i > l_zp and i <= l_zp+l_zs:
        Tg_3 = T_out_zp + ((T_out_zs - T_out_zp)/l_zs)*(i-l_zp)
        Tg[index] = Tg_3
    if i > l_zp+l_zs and i <= l_zp+l_zs+l_zd:
        Tg_4 = T_out_zs + ((T_out_zd - T_out_zs)/l_zd)*(i - l_zp - l_zs)
        Tg[index] = Tg_4
    index = index + 1

plt.figure(2, figsize=(12, 7), dpi=80)
plt.plot(x_t, Tg, 'r')
plt.title('Temperatura dos Gases ao Longo do Tubo de Chama')
plt.grid()
plt.ylabel('Temperatura (K)')
plt.xlabel('Distância da face do tubo de chama (mm)')
plt.legend()
plt.show()
