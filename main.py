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

printer. print_areas(A_ref_aero_EM, A_ref_aero_MA, A_ref_aero_CRU, A_ref_aero_IDLE, A_ref_quim_EM, A_ref_quim_MA, A_ref_quim_CRU, A_ref_quim_IDLE,
                     D_ref_aero_EM, D_ref_aero_MA, D_ref_aero_CRU, D_ref_aero_IDLE,
                     D_ref_quim_EM, D_ref_quim_MA, D_ref_quim_CRU, D_ref_quim_IDLE,
                     A_ft_aero_EM, A_ft_aero_MA, A_ft_aero_CRU, A_ft_aero_IDLE,
                     A_ft_quim_EM, A_ft_quim_MA, A_ft_quim_CRU, A_ft_quim_IDLE,
                     D_ft_aero_EM, D_ft_aero_MA, D_ft_aero_CRU, D_ft_aero_IDLE,
                     D_ft_quim_EM, D_ft_quim_MA, D_ft_quim_CRU, D_ft_quim_IDLE)

# Valores errados
# A_ref_maior = 0.05890
# D_ref_maior = 0.1142
# A_ft_maior = A_ref_maior * 0.7
# D_ft_maior = 0.0799

A_ref_maior = A_ref_quim_IDLE
D_ref_maior = D_ref_quim_IDLE
A_ft_maior = A_ft_quim_IDLE
D_ft_maior = D_ft_quim_IDLE

print(f'A_ref: {A_ref_maior:.4f} \nD_ref: {D_ref_maior:.4f}\nA_ft:  {A_ft_maior:.4f}\nD_ft:  {D_ft_maior:.4f}')

# Determinação dos comprimentos da camara

l_zr = (1 / 2) * D_ft_maior
l_zp = (3 / 4) * D_ft_maior
l_zs = (1 / 2) * D_ft_maior
l_zd = (3 / 2) * D_ft_maior
l_cc = l_zp + l_zs + l_zd


printer.print_comprimentos(l_zr, l_zp, l_zs, l_zd, l_cc)

# CÁLCULOS DIFUSOR

# eq. 26
V3 = 150  # Nova velocidade de entrada
# eq. 22
A3 = (m3_EM * T3_EM * R_ar2)/(V3 * P3_EM)

# eq. 25
A_an = A_ref_maior - A_ft_maior

m_dot_an = m3_EM - (3/4*m_dot_zp)

# eq. 24
A0 = m3_EM * A_an/m_dot_an


# eq. 27
AR = A0 / A3

# eq. 3.22
# A0=(a^2+2ab+2ac+b^2+2bc+c^2)-(a^2+2ab-2ac+b^2-2bc+c^2)
D0 = A0/(3.14159*(D_int+D_ref_maior))


# eq. 3.23
# A0=(a^2+2ab+2ac+b^2+2bc+c^2)-(a^2+2ab-2ac+b^2-2bc+c^2)=> A0=4ac+4bc=4D3*(Dint+Dref) => A0*()
D3 = A3/(3.14159*(D_int+D_ref_maior))


# CONSIDERANDO A PERDA DE PRESSÃO DO DIFUSOR = 1% = dPdif/P3 - 3.25 (Filipe)
tan_phi = math.pow((0.01*(A3**2)) /
                   (1.75*R_ar2*((m3_EM*(T3_EM**(1/2))/P3_EM)**2)*((1-A3/A0)**2)), 1/1.22)

# 3.26
L_dif = (D0/2 - D3/2) / (tan_phi)

print('Area de entrada do difusor                   ', A3)
print('Diâmetro de entrada do difusor              ', D3)
print('Area de saída do difusor       ', A0)
print('Diâmetro de saída do difusor  ', D0)
print('Ângulo do difusor             ', tan_phi*180/3.141516)
print('Comprimento do difusor           ', L_dif)


# CÁLCULOS TURBILHONADOR

# eq. professor - nmr injetores
N_inj = round(math.pi * (D_int+D_ref_maior)/D_ft_maior)

# eq. qtd pás (10 por injetor)
N_pas = N_inj * 10

# eq. 31
l_zr = D_ft_maior / 2
D_sw = l_zr / 2

# eq. 32 - para 5% do total no turbilhonador
m_dot_sw = 0.05 * m3_EM
m_dot_zr = 3 * m_dot_sw

# eq. 36
A0 = m3_EM * (A_ref_maior - A_ft_maior) / (m3_EM - (0.5*m_dot_zp))  # eq. 24
perda_pressao_entrada = 0.25*math.pow((A_ref_maior/A0), 2)  # eq. 39
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
print('O diâmetro do turbilhonador é       ', D_sw)
print('A vazão massica do turbilhonador é  ', m_dot_sw)
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

bal_massa = m3_EM - (m_dot_zp + (m_dot_zs - m_dot_zp) + (m_dot_zd - m_dot_zs))
print("🐍 File: projeto-combustao/main.py | Line: 281 | undefined ~ bal_massa", bal_massa)

print("\n")


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
print("\n")

x_t = np.linspace(1e-4, l_cc, 100)
Tg = np.zeros(len(x_t))
index = 0
countzr = 0
countzp = 0
countzs = 0
countzd = 0

# Versão professor
for i in x_t:
    if i >= 0 and i <= l_zr:
        Tg_1 = eq.tg_zr(T_med_zr)
        Tg[index] = Tg_1
        countzr = countzr+1

    if i > l_zr and i <= l_zp:
        Tg_2 = eq.tg_zp(T_out_zp, l_zp, T3_EM, i)
        Tg[index] = Tg_2
        countzp = countzp+1

    if i > l_zp and i <= l_zp+l_zs:
        Tg_3 = eq.tg_zs(T_out_zp, T_out_zs, l_zs, l_zp, i)
        Tg[index] = Tg_3
        countzs = countzs+1

    if i > l_zp+l_zs and i <= l_zp+l_zs+l_zd:
        Tg_4 = eq.tg_zd(T_out_zs, T_out_zd, l_zd, l_cc, i)
        Tg[index] = Tg_4
        countzd = countzd+1

    index = index + 1


Tg_smoothed = gaussian_filter1d(Tg, sigma=2)

plt.figure(2, figsize=(12, 7), dpi=80)
plt.plot(x_t, Tg_smoothed, 'r')
plt.title('Temperatura dos Gases ao Longo do Tubo de Chama')
plt.grid()
plt.vlines(l_zr, 1000, 2900, colors='b', linestyles='--',
           label='Limite Zona de Recirculação')
plt.vlines((l_zp), 1000, 2900, colors='g',
           linestyles='--', label='Limite Zona Primária')
plt.vlines((l_zp+l_zs), 1000, 2900, colors='r',
           linestyles='--', label='Limite Zona Secundária')
plt.vlines((l_cc), 1000, 2900, colors='m',
           linestyles='--', label='Limite Zona de Diluição')
plt.ylim(1000, 2900)
plt.ylabel('Temperatura (K)')
plt.xlabel('Distância da face do tubo de chama (mm)')
plt.legend()
plt.show()


#   Fluxo de massas

mg = np.zeros(len(x_t))
m_an = np.zeros(len(x_t))
m_an_zr = np.zeros(len(x_t))
m_an_zp = np.zeros(len(x_t))
m_an_zs = np.zeros(len(x_t))
m_an_zd = np.zeros(len(x_t))


index = 0


for i in x_t:
    if i >= 0 and i <= l_zr:
        mg_zr = eq.mg_zr(m_dot_zp)
        mg[index] = mg_zr
        m_an[index] = m3_EM - mg_zr
        m_an_zr[index] = m3_EM - mg_zr

    if i > l_zr and i <= l_zp:
        mg_zp = eq.mg_zp(mg_zr, m_dot_zp, l_zp, l_zr, i)
        mg[index] = mg_zp
        m_an[index] = m3_EM - mg_zp
        m_an_zp[index] = m3_EM - mg_zp

    if i > l_zp and i <= l_zp+l_zs:
        mg_zs = eq.mg_zs(m_dot_zs, mg_zp, l_zp, l_zs, i)
        mg[index] = mg_zs
        m_an[index] = m3_EM - mg_zs
        m_an_zs[index] = m3_EM - mg_zs

    if i > l_zp+l_zs and i <= l_zp+l_zs+l_zd:
        mg_zd = eq.mg_zd(mg_zs, m_dot_zd, l_zp, l_zs, l_zd, i)
        mg[index] = mg_zd
        m_an[index] = m3_EM - mg_zd
        m_an_zd[index] = m3_EM - mg_zd

    index = index + 1

mg_smoothed = gaussian_filter1d(mg, sigma=1)
m_an_smoothed = gaussian_filter1d(m_an, sigma=2)


plt.figure(2, figsize=(12, 7), dpi=80)
plt.plot(x_t, mg_smoothed, 'r')
plt.plot(x_t, m_an_smoothed, 'r')

plt.title('Temperatura dos Gases ao Longo do Tubo de Chama')
plt.vlines(l_zr, 0, 18, colors='b', linestyles='--',
           label='Limite Zona de Recirculação')
plt.vlines((l_zp), 0, 18, colors='g',
           linestyles='--', label='Limite Zona Primária')
plt.vlines((l_zp+l_zs), 0, 18, colors='r',
           linestyles='--', label='Limite Zona Secundária')
plt.vlines((l_cc), 0, 18, colors='m',
           linestyles='--', label='Limite Zona de Diluição')
plt.grid()
plt.ylabel('Fluxo de Massa [kg/s]')
plt.xlabel('Comprimento da Câmara de Combustão [mm]')
plt.title('Fluxo de Massa por Zona')
plt.legend()
plt.show()

T_max = 1100

Tw1, Tw2, x, nr = sp.symbols('Tw1, Tw2,x, nr')

lb = 0.9 * D_ft_maior

L = 1.7  # definido pelo professor

# q = phi_zd  # DUVIDA

epsilon_w = 0.7  # Adotado o Nimonic(Aluminio) como material da parede

sigma = 5.67*10**(-8)  # Boltzmann

kw = 26  # [W/m.K] Condutividade termica do material
tw = 0.0005  # espessura parede

D_an = 2 * (D_ref_maior-D_ft_maior)  # Tg ITA
A_an = (A_ref_maior - A_ft_maior)


Tw_in = np.zeros(len(x_t))

ka = ((5.92657 * 10**(-4)) + (9.80957 * 10**(-5)) * T3_EM -
      (4.89398 * 10**(-8)) * (T3_EM ** 2) + 1.5011410 * 10**(-11)*(T3_EM ** 3))

for i in range(len(x_t)):

    if i < 27:
        q = phi_zp
    elif i < 45:
        q = phi_zs
    elif i < 60:
        q = phi_zd

    epsilon_g = (1 - sp.exp(-0.29 * P3_EM * L *
                            ((q * lb)**(0.5)) * ((Tg[i])**(-1.5))))

    mi_ar = (0.03863 + 0.00749*T3_EM - 5.8564*(10**-6)*(T3_EM**2) +
             2.7769*(10**-9)*(T3_EM**3) - 4.600774*(10**-13)*(T3_EM**4))*(10**-5)
    mi_g = (0.03863 + 0.00749*Tg[i] - 5.8564*(10**-6)*(Tg[i]**2) + 2.7769 *
            (10**-9)*(Tg[i]**3) - 4.600774*(10**-13)*(Tg[i]**4))*(10**-5)

    kg = (5.92657 * 10**(-4)) + (9.80957 * 10**(-5)) * Tg[i] * (-4.89398 *
                                                                10 ** (-8)) * (Tg[i] ** 2) + 1.5011410 * 10**(-11) * Tg[i] ** 3

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

    if i <= 28:
        eq = [R1 + C1_zp - K_12, R2 + C2 - K_12]
    else:
        eq = [R1 + C1 - K_12, R2 + C2 - K_12]

    Tw_ext = sp.nsolve(eq, (Tw1, Tw2), (1100, 700))

    Tw_in[i] = Tw_ext[0]


Tw_in_smoothed = gaussian_filter1d(Tw_in, sigma=2)


plt.figure(2, figsize=(12, 7), dpi=80)
plt.plot(x_t, Tg_smoothed, 'b')
plt.plot(x_t, Tw_in_smoothed, 'r')
plt.title('Temperatura dos Gases ao Longo do Tubo de Chama')
plt.vlines(l_zr, 800, 2800, colors='b', linestyles='--',
           label='Limite Zona de Recirculação')
plt.vlines((l_zp), 800, 2800, colors='g',
           linestyles='--', label='Limite Zona Primária')
plt.vlines((l_zp+l_zs), 800, 2800, colors='darkorange',
           linestyles='--', label='Limite Zona Secundária')
plt.vlines((l_cc), 800, 2800, colors='m',
           linestyles='--', label='Limite Zona de Diluição')
plt.ylim(800, 2800)
plt.grid()
plt.ylabel('Temperatura (K)')
plt.xlabel('Distância da face do tubo de chama (mm)')
plt.legend()
plt.show()


# Filme de resfriamento
s = 0.0001  # tabela gasturb 335prd_den_vel_an(m_f,A_f)
t = 0.0001
Tw_in_res = np.zeros(len(x_t))


A_f = (2 * 3.1415926535897932384626 * s * (D_ref_maior + D_ft_maior))

for x in range(len(x_t)):

    if x < 27:
        q = phi_zp
    elif x < 45:
        q = phi_zs
    elif x < 60:
        q = phi_zd

    epsilon_g = (1 - sp.exp(-0.29 * P3_EM * L *
                            ((q * lb)**(0.5)) * ((Tg[x])**(-1.5))))

    mi_ar = (0.03863 + 0.00749*T3_EM - 5.8564*(10**-6)*(T3_EM**2) +
             2.7769*(10**-9)*(T3_EM**3) - 4.600774*(10**-13)*(T3_EM**4))*(10**-5)
    mi_g = (0.03863 + 0.00749*Tg[x] - 5.8564*(10**-6)*(Tg[x]**2) + 2.7769 *
            (10**-9)*(Tg[x]**3) - 4.600774*(10**-13)*(Tg[x]**4))*(10**-5)

    kg = (5.92657 * 10**(-4)) + (9.80957 * 10**(-5)) * Tg[x] * (-4.89398 *
                                                                10 ** (-8)) * (Tg[x] ** 2) + 1.5011410 * 10**(-11) * Tg[x] ** 3

    m_f = m_an[x] * (A_f/A_an)

    product_density_velocity_an = m_f / A_f
    product_density_velocity_g = mg[x] / A_ft_maior

    mzinho = product_density_velocity_an / product_density_velocity_g

    Rex = ((product_density_velocity_an * i) / mi_ar)

    if mzinho >= 0.5 and mzinho <= 1.3:
        nr = (1.1 * (mzinho ** 0.65) * ((mi_ar / mi_g)**0.15)
              * ((i/s)**(-0.2)) * ((t/s)**(-0.2)))
        Tw_ad = Tg[x] - nr * (Tg[x] - T3_EM)

        C1 = (0.069 * (ka/i) * (Rex**0.7) * (Tw_ad - Tw1))

    elif mzinho > 1.3 and mzinho < 4:
        nr = (1.28 * ((mi_ar / mi_g)**0.15) *
              ((i/s)**(-0.2)) * ((t/s)**(-0.2)))
        Tw_ad = Tg[x] - nr * (Tg[x] - T3_EM)
        C1 = (0.1 * (ka/i) * (Rex**0.8) * ((i/s)**(-0.36)) * (Tw_ad - Tw1))

    else:
        nr = 1
        Tw_ad = Tg[x] - nr * (Tg[x] - T3_EM)
        C1 = (0.1 * (ka/i) * (Rex**0.8) * ((i/s)**(-0.36)) * (Tw_ad - Tw1))

    R1 = (0.5 * sigma * (1 + epsilon_w) * epsilon_g *
          (Tg[x]**1.5) * ((Tg[x]**2.5)-(Tw1**2.5)))
    R2 = (0.6 * sigma * (Tw2 ** 4 - T3_EM**4))  # Adotado para o alunimio
    K_12 = (kw/tw)*(Tw1 - Tw2)

    C1_zp = (0.017 * (kg / (D_ft_maior ** 0.2)) *
             ((mg[x] / (A_ft_maior * mi_g))**0.8) * (Tw_ad - Tw1))

    C2 = (0.02 * (ka / (D_an ** 0.2)) *
          ((m_an[x] / (A_an * mi_ar))**0.8) * (Tw2-T3_EM))

    if i <= 28:
        eq = [R1 + C1_zp - K_12, R2 + C2 - K_12]
    else:
        eq = [R1 + C1 - K_12, R2 + C2 - K_12]

    Tw_ext_res = sp.nsolve(eq, (Tw1, Tw2), (1100, 50))

    Tw_in_res[x] = Tw_ext_res[0]

Tw_in_res_smoothed = gaussian_filter1d(Tw_in_res, sigma=2)


plt.figure(2, figsize=(12, 7), dpi=80)
plt.plot(x_t, Tg_smoothed, 'r')
plt.plot(x_t, Tw_in_smoothed, 'b')
plt.plot(x_t, Tw_in_res_smoothed, 'darkorange')
plt.title('Temperatura dos Gases ao Longo do Tubo de Chama')
# plt.vlines(l_zr, 800, 2800, colors='b', linestyles='--',
#            label='Limite Zona de Recirculação')
# plt.vlines((l_zp), 800, 2800, colors='g',
#            linestyles='--', label='Limite Zona Primária')
# plt.vlines((l_zp+l_zs), 800, 2800, colors='darkorange',
#            linestyles='--', label='Limite Zona Secundária')
# plt.vlines((l_cc), 800, 2800, colors='m',
#            linestyles='--', label='Limite Zona de Diluição')
# plt.ylim(800, 2800)
plt.grid()
plt.ylabel('Temperatura (K)')
plt.xlabel('Distância da face do tubo de chama (mm)')
plt.legend()
plt.show()


# Calculo dos orificios de admissao

m_an_zr = m_an_zr[m_an_zr != 0]
m_an_zp = m_an_zp[m_an_zp != 0]
m_an_zs = m_an_zs[m_an_zs != 0]
m_an_zd = m_an_zd[m_an_zd != 0]


mf_zr = m_an_zr * (A_f/A_an)
mf_zp = m_an_zp * (A_f/A_an)
mf_zs = m_an_zs * (A_f/A_an)
mf_zd = m_an_zd * (A_f/A_an)

mh_zr = m_dot_zr - mf_zr
mh_zp = m_dot_zp - mf_zp
mh_zs = m_dot_zs - mf_zs
mh_zd = m_dot_zd - mf_zd

Cd_hRZ = np.zeros((100))
Cd_hPZ = np.zeros((100))
Cd_hSZ = np.zeros((100))
Cd_hDZ = np.zeros((100))
Cd_hRZ[0] = 0.3
Cd_hPZ[0] = 0.3
Cd_hSZ[0] = 0.3
Cd_hDZ[0] = 0.3
delta = 0.8
deltaPh = 0.06*P3_EM


#   ZONA RECIRCULACAO
m_h_RZdot = m_dot_zp*0.5
m_aRZ_dot = m_an_zr
n = 0
Delta = 1
while abs(Delta) > 0.0001:
    Ah_RZ = ((143.5*m_h_RZdot**2*T3_EM) /
             (Cd_hRZ[n]**2*P3_EM**2*deltaPh/P3_EM))**0.5
    alpha_RZ = Ah_RZ/A_an
    beta_RZ = m_h_RZdot/m_aRZ_dot[0]
    mi_RZ = beta_RZ/alpha_RZ
    K_RZ = 1 + delta**2*(2*mi_RZ**2 + (4*mi_RZ**4 +
                                       (mi_RZ**2/delta**2) * (4*beta_RZ - beta_RZ**2))**0.5)
    n = n + 1
    Cd_hRZ[n] = (K_RZ - 1)/(delta*(4*K_RZ**2 - K_RZ*(2 - beta_RZ)**2)**0.5)
    Delta = Cd_hRZ[n] - Cd_hRZ[n-1]

print("Área Total dos furos (Ah_ZR):", round(Ah_RZ*10**6, 1), "mm²")


#   ZONA PRIMARIA
m_h_PZdot = m_dot_zp - m_dot_zr
m_aPZ_dot = m_an_zp
n = 0
Delta = 1

while abs(Delta) > 0.0001:
    Ah_PZ = ((143.5*m_h_PZdot**2*T3_EM) /
             (Cd_hPZ[n]**2*P3_EM**2*deltaPh/P3_EM))**0.5
    alpha_PZ = Ah_PZ/A_an
    beta_PZ = m_h_PZdot/m_aPZ_dot[0]
    mi_PZ = beta_PZ/alpha_PZ
    K_PZ = 1 + delta**2*(2*mi_PZ**2 + (4*mi_PZ**4 + (mi_PZ**2/delta**2) *
                                       (4*beta_PZ - beta_PZ**2))**0.5)
    n = n + 1
    Cd_hPZ[n] = (K_PZ - 1)/(delta*(4*K_PZ**2 - K_PZ*(2 - beta_PZ)**2)**0.5)
    Delta = Cd_hPZ[n] - Cd_hPZ[n-1]


print("Área Total dos furos (Ah_PZ):", round(Ah_PZ*10**6, 1), "mm²")


#  ZONA SECUNDARIA
m_h_SZdot = m_dot_zs - (m_dot_zp - m_dot_zr)
m_aSZ_dot = m_an_zs
n = 0
Delta = 1
while abs(Delta) > 0.0001:
    Ah_SZ = ((143.5*m_h_SZdot**2*T3_EM) /
             (Cd_hSZ[n]**2*P3_EM**2*deltaPh/P3_EM))**0.5
    alpha_SZ = Ah_SZ/A_an
    beta_SZ = m_h_SZdot/m_aSZ_dot[0]
    mi_SZ = beta_SZ/alpha_SZ
    K_SZ = 1 + delta**2*(2*mi_SZ**2 + (4*mi_SZ**4 + (mi_SZ**2/delta**2) *
                                       (4*beta_SZ - beta_SZ**2))**0.5)
    n = n + 1
    Cd_hSZ[n] = (K_SZ - 1)/(delta*(4*K_SZ**2 - K_SZ*(2 - beta_SZ)**2)**0.5)
    Delta = Cd_hSZ[n] - Cd_hSZ[n-1]
print("Área Total dos furos (Ah_SZ):", round(Ah_SZ*10**6, 1), "mm²")


#   ZONA DILUICAO
m_h_DZdot = m_dot_zd - (m_dot_zs - (m_dot_zp - m_dot_zr))
m_aDZ_dot = m_an_zd
n = 0
Delta = 1
while abs(Delta) > 0.0001:
    Ah_DZ = ((143.5*m_h_DZdot**2*T3_EM) /
             (Cd_hDZ[n]**2*P3_EM**2*deltaPh/P3_EM))**0.5
    alpha_DZ = Ah_DZ/A_an
    beta_DZ = m_h_DZdot/m_aDZ_dot[0]
    mi_DZ = beta_DZ/alpha_DZ
    K_DZ = 1 + delta**2*(2*mi_DZ**2 + (4*mi_DZ**4 + (mi_DZ**2/delta**2) *
                                       (4*beta_DZ - beta_DZ**2))**0.5)
    n = n + 1
    Cd_hDZ[n] = (K_DZ - 1)/(delta*(4*K_DZ**2 - K_DZ*(2 - beta_DZ)**2)**0.5)
    Delta = Cd_hDZ[n] - Cd_hDZ[n-1]

print("Área Total dos furos (Ah_DZ):", round(Ah_DZ*10**6, 1), "mm²")
print("\n")

Ah_total = Ah_RZ+Ah_PZ+Ah_SZ+Ah_DZ

Nc_h = np.array([30, 18, 18, 18], dtype=int)  # Números em cada seção
Dc_h = np.array([(4*Ah_RZ/(math.pi*Nc_h[0]))**0.5, (4*Ah_PZ/(math.pi*Nc_h[1]))**0.5, (4*Ah_SZ /
                                                                                      (math.pi*Nc_h[2]))**0.5, (4*Ah_DZ/(math.pi*Nc_h[3])) ** 0.5])  # Diâmetro de cada furo, em cada fileira
print(f'Diâmetro total (somados) dos furos (mm):{ Dc_h*10**3}', '\n')
Pc = np.array([30, 80, 150, 350], dtype=int)  # Posição de cada fileira


# Nh_ext = ((3.1415926535897932384626 *
#            (D_int + D_ref_maior + D_ft_maior)) / (2 * Dc_h))
# for i in range(len(Nh_ext)):
#     print(
#         f'O numero maximo de orificios externos da zona {i+1} é: {math.floor(Nh_ext[i])}')
# print('\n')

# Nh_int = ((3.1415926535897932384626 *
#            (D_int + D_ref_maior - D_ft_maior)) / (2 * Dc_h))

# for i in range(len(Nh_int)):
#     print(
#         f'O numero maximo de orificios internos da zona {i+1} é: {math.floor(Nh_int[i])}')
# print('\n')

# A_an_ext = ((3.1415926535897932384626/4) * (((D_int + (2 * D_ref_maior)) ** 2) -
#             ((D_int+D_ref_maior+D_ft_maior)**2)))

# A_an_int = ((3.1415926535897932384626/4) * (((D_int+D_ref_maior+D_ft_maior)**2) - (D_int**2)))


# Ah_an_ext_zr = (A_an_ext / Ah_total) * Ah_RZ
# Ah_an_ext_zp = (A_an_ext / Ah_total) * Ah_PZ
# Ah_an_ext_zs = (A_an_ext / Ah_total) * Ah_SZ
# Ah_an_ext_zd = (A_an_ext / Ah_total) * Ah_DZ


# Ah_an_int_zr = (A_an_int / (A_an_ext+A_an_int)) * Ah_RZ
# Ah_an_int_zp = (A_an_int / (A_an_ext+A_an_int)) * Ah_PZ
# Ah_an_int_zs = (A_an_int / (A_an_ext+A_an_int)) * Ah_SZ
# Ah_an_int_zd = (A_an_int / (A_an_ext+A_an_int)) * Ah_DZ
