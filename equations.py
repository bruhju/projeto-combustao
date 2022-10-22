import sympy  as sp
import numpy as np 
import pandas as pd 
import math  
import matplotlib.pyplot as plt
import scipy as spy

def A_ref_aerodinamico (R_ar, m_dot_3, T3, P3, perda_pressao_total, fator_perda_pressao):
    #perda_pressao_total = deltaP3-4/P3
    #fator_perda_pressao = deltaP3-4/qRef
    A_ref_aero = math.sqrt( R_ar * ( math.pow(( m_dot_3 * math.sqrt(T3)) / P3 ,2)) * ( fator_perda_pressao / perda_pressao_total ))
    return A_ref_aero

def A_ref_quimico (m_dot_3, T3, P3, theta, D_ref, b):
    A_ref_q = ( theta * m_dot_3 ) / (  math.pow(P3,1.75) *  math.pow(D_ref,0.75) * math.exp(T3 / b) )
    return A_ref_q

def D_ref (A_ref, D_int):
    D_ref = ( math.sqrt(((4 * A_ref)/math.pi) + math.pow(D_int,2)) - D_int ) / 2
    return D_ref

