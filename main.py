import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
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
R_ar2 = 287  # J/kg*K
theta = 73000000
D_int = 0.05
phi_zp = 1.05  # deprecated
phi_zs = 8
phi_estq = 0.06818
b = 170*(2-math.log(phi_zp))
m_dot_zp = 0.238578

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

A_ref_maior = 0.05890
D_ref_maior = 0.1142
A_ft_maior = A_ref_maior * 0.7
D_ft_maior = 0.0799

print(f'A_ref: {A_ref_maior:.4f} \nD_ref: {D_ref_maior:.4f}\nA_ft:  {A_ft_maior:.4f}\nD_ft:  {D_ft_maior:.4f}')

# Determinação dos comprimentos da camara

l_zr = (1 / 2) * D_ft_maior
l_zp = (3 / 4) * D_ft_maior
l_zs = (1 / 2) * D_ft_maior
l_zd = (3 / 2) * D_ft_maior
l_cc = l_zp + l_zs + l_zd


printer.print_comprimentos(l_zp, l_zs, l_zd, l_cc)

# Calculos difusor

# Eq. 26
A3 = 0.096 / 1  # saida do compressor/número de combustores anular(1)

# eq. 22
V3 = (m3_EM * T3_EM * R_ar2)/(A3 * P3_EM)
print('Velocidade de entrada: ', V3)

# condicional
if V3 < 150:
    print("Não é necessário difusor porque a velocidade de entrada do ar não é superior a 75 m/s")


# CÁLCULOS TURBILHONADOR

# Eq. professor - nmr injetores
N_inj = round(math.pi * (D_int+D_ref_maior)/D_ft_maior)

# Eq. qtd pás (10 por injetor)
N_pas = N_inj * 10

# Eq. 31
D_sw = l_zr / 2

# Eq. 32 - para 5% do total no turbilhonador
m_dot_sw = 0.05 * m3_EM
m_dot_zr = 3 * m_dot_sw

# Eq. 36
A0 = m3_EM * (A_ref_maior - A_ft_maior) / (m3_EM - (0.5*m_dot_zp))  # Eq. 24
perda_pressao_entrada = 0.25*math.pow((A_ref_maior/A0), 2)  # Eq. 39
perda_pressao_difusor = 0
perda_pressao_turbilhonador = fator_perda_pressao - \
    perda_pressao_entrada - perda_pressao_difusor

# eq. 35 - ângulo de 55°
K_sw = 1.3
sec_2_Bsw = 3.0396
A_sw = sec_2_Bsw * ((K_sw*math.pow(m_dot_sw, 2)*math.pow(A_ref_maior, 2)) /
                    (perda_pressao_turbilhonador*math.pow(m3_EM, 2))+math.pow(A_ft_maior, 2))

A_sw_corrigido = A_sw*1.5

# eq. 41 - 10%
D_at = 0.1 * D_ref_maior
R_at = D_at/2

# eq.

print('O número de injetores é             ', N_inj)
print('O número de pás é                   ', N_pas)
print('O comprimento da ZR é               ', l_zr)
print('O diâmetro do turbilhonador é       ', D_sw)
print('A vazão massica do turbilhonador é  ', m_dot_sw)
print('A vazão massica da ZR é             ', m_dot_zr)
print('A área do turbilhonador é           ', A_sw)
print('A área corrigida do turbilhonador é ', A_sw_corrigido)
print('O diâmetro do atomizador é          ', D_at)
print('Relação diâmetro turbilhonador e diâmetro tubo de chama é ',
      D_sw/D_ft_maior*100, '%')


# A velocidade de entrada do ar na camara nao deve ultrapassar 75 m/s para que a combustao seja estavel numa grande faixa de razoes ar/combustivel


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
air_zs_percent = (phi_global_EM/phi_zs)
m_dot_zs = air_zs_percent*m3_EM


# Fluxo de massa na zona de diluição

air_zd_percent = 1 - (air_zp_percent + air_zs_percent +
                      air_arrefecimento_percent)

phi_zd = phi_global_EM

air_zd_percent = 1

m_dot_zd = air_zd_percent*m3_EM

massasDF = pd.DataFrame({'m_dot': [m_dot_zp, m_dot_zs, m_dot_zd, m_dot_arref],
                         })
massasDF.index = ['Zona Primaria', "Zona Secundaria",
                  "Zona Diluicao", "Arrefecimento"]
print('\n', massasDF, '\n')

bal_massa = m3_EM - (m_dot_zp + m_dot_zs + m_dot_zd + m_dot_arref)
print("🐍 File: projeto-combustao/main.py | Line: 279 | undefined ~ bal_massa", bal_massa)

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

x_t = np.linspace(1e-4, l_cc, 100)
Tg = np.zeros(len(x_t))
index = 0

# Versão professor
for i in x_t:
    if i >= 0 and i <= l_zr:
        Tg_1 = eq.tg_zr(T_med_zr)
        Tg[index] = Tg_1

    if i > l_zr and i <= l_zp:
        Tg_2 = eq.tg_zp(T_out_zp, l_zp, T3_EM, i)
        Tg[index] = Tg_2

    if i > l_zp and i <= l_zp+l_zs:
        Tg_3 = eq.tg_zs(T_out_zp, T_out_zs, l_zs, l_zp, i)
        Tg[index] = Tg_3

    if i > l_zp+l_zs and i <= l_zp+l_zs+l_zd:
        Tg_4 = eq.tg_zd(T_out_zs, T_out_zd, l_zd, l_cc, i)
        Tg[index] = Tg_4

    index = index + 1


Tg_smoothed = gaussian_filter1d(Tg, sigma=2)

# plt.figure(2, figsize=(12, 7), dpi=80)
# plt.plot(x_t, Tg_smoothed, 'r')
# plt.title('Temperatura dos Gases ao Longo do Tubo de Chama')
# plt.grid()
# plt.vlines(l_zr, 1000, 2900, colors='b', linestyles='--',
#            label='Limite Zona de Recirculação')
# plt.vlines((l_zp), 1000, 2900, colors='g',
#            linestyles='--', label='Limite Zona Primária')
# plt.vlines((l_zp+l_zs), 1000, 2900, colors='r',
#            linestyles='--', label='Limite Zona Secundária')
# plt.vlines((l_cc), 1000, 2900, colors='m' l_zp+l_zs+l_zd, ,
#            linestyles='--', label='Limite Zona de Diluição')
# plt.ylim(1000, 2900)
# plt.ylabel('Temperatura (K)')
# plt.xlabel('Distância da face do tubo de chama (mm)')
# plt.legend()
# plt.show()

mg = np.zeros(len(x_t))
m_an = np.zeros(len(x_t))
index = 0


for i in x_t:
    if i >= 0 and i <= l_zr:
        mg_zr = eq.mg_zr(m_dot_zp)
        mg[index] = mg_zr
        m_an[index] = m3_EM - mg_zr

    if i > l_zr and i <= l_zp:
        mg_zp = eq.mg_zp(mg_zr, m_dot_zp, l_zp, l_zr, i)
        mg[index] = mg_zp
        m_an[index] = m3_EM - mg_zp

    if i > l_zp and i <= l_zp+l_zs:
        mg_zs = eq.mg_zs(m_dot_zs, mg_zp, l_zp, l_zs, i)
        print("🐍 File: projeto-combustao/main.py | Line: 381 | undefined ~ mg_zp, mg_zs", mg_zp, " ", mg_zs)
        mg[index] = mg_zs
        m_an[index] = m3_EM - mg_zs

    if i > l_zp+l_zs and i <= l_zp+l_zs+l_zd:
        mg_zd = eq.mg_zd(mg_zs, m_dot_zd, l_zp, l_zs, l_zd, i)
        mg[index] = mg_zd
        m_an[index] = m3_EM - mg_zd

    index = index + 1

mg_smoothed = gaussian_filter1d(mg, sigma=1)
m_an_smoothed = gaussian_filter1d(m_an, sigma=2)


plt.figure(2, figsize=(12, 7), dpi=80)
plt.plot(x_t, mg_smoothed, 'r')
plt.title('Temperatura dos Gases ao Longo do Tubo de Chama')
plt.vlines(l_zr, 0, 6, colors='b', linestyles='--',
           label='Limite Zona de Recirculação')
plt.vlines((l_zp), 0, 6, colors='g',
           linestyles='--', label='Limite Zona Primária')
plt.vlines((l_zp+l_zs), 0, 6, colors='r',
           linestyles='--', label='Limite Zona Secundária')
plt.vlines((l_cc), 0, 6, colors='m',
           linestyles='--', label='Limite Zona de Diluição')
plt.grid()
plt.ylabel('Fluxo de Massa [kg/s]')
plt.xlabel('Comprimento da Câmara de Combustão [mm]')
plt.title('Fluxo de Massa por Zona')
plt.legend()
plt.show()

T_max = 1100

Tw1, Tw2, epsilon_g, Re, m_esp_flui, v_fluid, D_flux_tb, mi_fluid, x = sp.symbols(
    'Tw1, Tw2, epsilon_g, Re, m_esp_flui , v_fluid , D_flux_tb,  mi_fluid,x')

lb = 0.9 * D_ft_maior

# Re = (m_esp_flui * v_fluid * D_flux_tb) / mi_fluid

# # m_esp_flui => massa especifica do fluido
# # v_fluid => velocidade media do fluido
# # D_flux_tb => diametro para o fluxo do tubo
# # mi_fluid => viscosidade dinamica do fluido

# m = 0.0283 * sp.Pow(Re, (-0.2))

L = 1.7
q = air_zp_percent  # DUVIDA


epsilon_w = 0.7  # Adotado o Nimonic(Aluminio) como material da parede

sigma = 5.67*10**(-8)  # Boltzmann

kw = 26  # [W/m.K] Condutividade termica do material
tw = 0.0005  # espessura parede

D_an = 2 * (D_ref_maior-D_ft_maior)  # Tg ITA
A_an = (A_ref_maior - A_ft_maior)


Tw_in = np.zeros(len(x_t))

for i in range(50):
    epsilon_g = (1 - sp.exp(-0.29 * P3_EM * L *
                 ((q * lb)**(0.5)) * ((Tg[i])**(-1.5))))

    mi_ar = (0.03863 + 0.00749*T3_EM - 5.8564*(10**-6)*(T3_EM**2) +
             2.7769*(10**-9)*(T3_EM**3) - 4.600774*(10**-13)*(T3_EM**4))*(10**-5)
    mi_g = (0.03863 + 0.00749*Tg[i] - 5.8564*(10**-6)*(Tg[i]**2) + 2.7769 *
            (10**-9)*(Tg[i]**3) - 4.600774*(10**-13)*(Tg[i]**4))*(10**-5)

    kg = (5.92657 * 10**(-4)) + (9.80957 * 10**(-5)) * Tg[i] * (-4.89398 *
                                                                10 ** (-8)) * (Tg[i] ** 2) + 1.5011410 * 10**(-11) * Tg[i] ** 3
    ka = ((5.92657 * 10**(-4)) + (9.80957 * 10**(-5)) * T3_EM -
          (4.89398 * 10**(-8)) * (T3_EM ** 2) + 1.5011410 * 10**(-11)*(T3_EM ** 3))

    R1 = (0.5 * sigma * (1 + epsilon_w) * epsilon_g *
          (Tg[i]**1.5) * ((Tg[i]**2.5)-(Tw1**2.5)))
    R2 = (0.6 * sigma * (Tw2 ** 4 - T3_EM**4))  # Adotado para o alunimio
    K_12 = (kw/tw)*(Tw1 - Tw2)

    C1 = (0.02 * (kg / (D_ft_maior ** 0.2)) *
          ((mg[i] / (A_ft_maior * mi_g))**0.8) * (Tg[i] - Tw1))
    C1_zp = (0.017 * (kg / (D_ft_maior ** 0.2)) *
             ((mg[i] / (A_ft_maior * mi_g))**0.8) * (Tg[i] - Tw1))
    C2 = (0.02 * (ka / (D_an ** 0.2)) *
          ((m_an[i] / (A_an * mi_ar))**0.8) * (Tw2-T3_EM))

    eq = [R1 + C1 - K_12, R2 + C2 - K_12]
    Tw_ext = sp.nsolve(eq, (Tw1, Tw2), (1100, 50))

    Tw_in[i] = Tw_ext[0]

    print("🐍 File: projeto-combustao/main.py | Line: 473 | undefined ~ Tw_ext", Tw_ext)

Tw_in_smoothed = gaussian_filter1d(Tw_in, sigma=2)

# plt.figure(2, figsize=(12, 7), dpi=80)
# plt.plot(x_t, Tg_smoothed, 'b')
# plt.plot(x_t, Tw_in_smoothed, 'r')
# plt.title('Temperatura dos Gases ao Longo do Tubo de Chama')
# plt.grid()
# plt.ylabel('Temperatura (K)')
# plt.xlabel('Distância da face do tubo de chama (mm)')
# plt.legend()
# plt.show()
