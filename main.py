import sympy  as sp
import numpy as np 
import pandas as pd 
import math
import matplotlib.pyplot as plt
import scipy as spy
import equations as eq

P3_EM = 2000000 # MPA
T3_EM = 814 # K
m3_EM = 18.1 # kg/s
T4_EM = 1600 # K
mComb_EM = 0.427 # kg/s

P3_MA = 700000 # MPA
T3_MA = 707 # K
m3_MA = 6.8 # kg/s
T4_MA = 1387 # K
mComb_MA = 0.132 # kg/s

P3_CRU = 1800000 # MPA
T3_CRU = 1060 # K
m3_CRU = 14.2 # kg/s
T4_CRU = 1393 # K
mComb_CRU = 0.140 # kg/s

P3_IDLE = 150000 # MPA
T3_IDLE = 343 # K
m3_IDLE = 1.05 # kg/s
T4_IDLE = 703 # K
mComb_IDLE = 0.0091 # kg/s

R_ar = 143.5 # J/kg*K
theta = 73000000
D_int = 0.01717
phi_zp = 1.05
phi_zs = 8
b = 170*(2-math.log(phi_zp))

perda_pressao_total =  0.06 #deltaP3-4/P3
fator_perda_pressao = 20    #deltaP3-4/qRef

A_ref_aero_EM = eq.A_ref_aerodinamico(R_ar, m3_EM, T3_EM, P3_EM, perda_pressao_total, fator_perda_pressao)
A_ref_aero_MA = eq.A_ref_aerodinamico(R_ar, m3_MA, T3_MA, P3_MA, perda_pressao_total, fator_perda_pressao)
A_ref_aero_CRU = eq.A_ref_aerodinamico(R_ar, m3_CRU, T3_CRU, P3_CRU, perda_pressao_total, fator_perda_pressao)
A_ref_aero_IDLE = eq.A_ref_aerodinamico(R_ar, m3_IDLE, T3_IDLE, P3_IDLE, perda_pressao_total, fator_perda_pressao)

D_ref_aero_EM = eq.D_ref( A_ref_aero_EM, D_int)
D_ref_aero_MA = eq.D_ref( A_ref_aero_MA, D_int)
D_ref_aero_CRU = eq.D_ref( A_ref_aero_CRU, D_int)
D_ref_aero_IDLE = eq.D_ref( A_ref_aero_IDLE, D_int)

A_ref_quim_EM = eq.A_ref_quimico( m3_EM, T3_EM, P3_EM, theta, D_ref_aero_EM, b)
A_ref_quim_MA = eq.A_ref_quimico( m3_MA, T3_MA, P3_MA, theta, D_ref_aero_MA, b)
A_ref_quim_CRU = eq.A_ref_quimico( m3_CRU, T3_CRU, P3_CRU, theta, D_ref_aero_CRU, b)
A_ref_quim_IDLE = eq.A_ref_quimico( m3_IDLE, T3_IDLE, P3_IDLE, theta, D_ref_aero_IDLE, b)


D_ref_quim_EM = eq.D_ref( A_ref_quim_EM, D_int)
D_ref_quim_MA = eq.D_ref( A_ref_quim_MA, D_int)
D_ref_quim_CRU = eq.D_ref( A_ref_quim_CRU, D_int)
D_ref_quim_IDLE = eq.D_ref( A_ref_quim_IDLE, D_int)




