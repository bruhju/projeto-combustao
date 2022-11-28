import pandas as pd


def print_phis(phi_pobre_EM, phi_pobre_MA, phi_pobre_CRU, phi_pobre_IDLE, phi_rico_EM, phi_rico_MA, phi_rico_CRU, phi_rico_IDLE, phi_global_EM, phi_global_MA, phi_global_CRU,  phi_global_IDLE):
    phisDF = pd.DataFrame(
        {'phi_pobre': [phi_pobre_EM, phi_pobre_MA, phi_pobre_CRU, phi_pobre_IDLE],
         'phi_rico': [phi_rico_EM, phi_rico_MA, phi_rico_CRU, phi_rico_IDLE],
         'phi_global': [phi_global_EM, phi_global_MA, phi_global_CRU,  phi_global_IDLE]
         })
    phisDF.index = ['Expuxo Maximo', "Maxima Altitude", "Cruzeiro", "IDLE"]
    print('\n', phisDF, '\n')


def print_limitEquivalence(limit_equivalence_pobre_EM, limit_equivalence_pobre_MA, limit_equivalence_pobre_CRU, limit_equivalence_pobre_IDLE,
                           limit_equivalence_rico_EM, limit_equivalence_rico_MA, limit_equivalence_rico_CRU, limit_equivalence_rico_IDLE):
    limitEquivalenceDF = pd.DataFrame(
        {'limit_equivalence_pobre': [limit_equivalence_pobre_EM, limit_equivalence_pobre_MA, limit_equivalence_pobre_CRU, limit_equivalence_pobre_IDLE],
         'limit_equivalence_rico': [limit_equivalence_rico_EM, limit_equivalence_rico_MA, limit_equivalence_rico_CRU, limit_equivalence_rico_IDLE]
         })
    limitEquivalenceDF.index = ['Expuxo Maximo',
                                "Maxima Altitude", "Cruzeiro", "IDLE"]
    print('\n', limitEquivalenceDF, '\n')


def print_phisZP(phi_zp_EM, phi_zp_MA, phi_zp_CRU, phi_zp_IDLE):
    phiZPDF = pd.DataFrame(
        {'phi_zp': [phi_zp_EM, phi_zp_MA, phi_zp_CRU, phi_zp_IDLE],
         })
    phiZPDF.index = ['Expuxo Maximo',
                     "Maxima Altitude", "Cruzeiro", "IDLE"]
    print('\n', phiZPDF, '\n')


def print_areas(A_ref_aero_EM, A_ref_aero_MA, A_ref_aero_CRU, A_ref_aero_IDLE, A_ref_quim_EM, A_ref_quim_MA, A_ref_quim_CRU, A_ref_quim_IDLE,
                D_ref_aero_EM, D_ref_aero_MA, D_ref_aero_CRU, D_ref_aero_IDLE,
                D_ref_quim_EM, D_ref_quim_MA, D_ref_quim_CRU, D_ref_quim_IDLE,
                A_ft_aero_EM, A_ft_aero_MA, A_ft_aero_CRU, A_ft_aero_IDLE,
                A_ft_quim_EM, A_ft_quim_MA, A_ft_quim_CRU, A_ft_quim_IDLE,
                D_ft_aero_EM, D_ft_aero_MA, D_ft_aero_CRU, D_ft_aero_IDLE,
                D_ft_quim_EM, D_ft_quim_MA, D_ft_quim_CRU, D_ft_quim_IDLE):

    areasDF = pd.DataFrame({'Aref Aero': [A_ref_aero_EM, A_ref_aero_MA, A_ref_aero_CRU, A_ref_aero_IDLE],
                            'Aref Quim': [A_ref_quim_EM, A_ref_quim_MA, A_ref_quim_CRU, A_ref_quim_IDLE],
                            'Dref Aero': [D_ref_aero_EM, D_ref_aero_MA, D_ref_aero_CRU, D_ref_aero_IDLE],
                            'Dref Quim': [D_ref_quim_EM, D_ref_quim_MA, D_ref_quim_CRU, D_ref_quim_IDLE],
                            'Aft Aero': [A_ft_aero_EM, A_ft_aero_MA, A_ft_aero_CRU, A_ft_aero_IDLE],
                            'Aft Quim': [A_ft_quim_EM, A_ft_quim_MA, A_ft_quim_CRU, A_ft_quim_IDLE],
                            'Dft Aero': [D_ft_aero_EM, D_ft_aero_MA, D_ft_aero_CRU, D_ft_aero_IDLE],
                            'Dft Quim': [D_ft_quim_EM, D_ft_quim_MA, D_ft_quim_CRU, D_ft_quim_IDLE],
                            })
    areasDF.index = ['Expuxo Maximo', "Maxima Altitude", "Cruzeiro", "IDLE"]
    print('\n', areasDF, '\n')


def print_comprimentos(l_zr, l_zp, l_zs, l_zd, l_cc):
    comprimentosDF = pd.DataFrame(
        {'Comprimentos': [l_zr, l_zp, l_zs, l_zd, l_cc]})
    comprimentosDF.index = ['Zona Recirculacao', 'Zona Primaria',
                            "Zona Secundaria", "Zona Diluicao", "CC"]
    print('\n', comprimentosDF, '\n')
