import numpy
import scipy.optimize as sc_o
import Psat

R = 8.314472


def G_RT_main_l(T, P_bar, Vmix, x1,Tc_1, Pc_1_bar, Tc_2, Pc_2_bar, m_1, m_2):
    Pc_1 = Pc_1_bar*100
    Pc_2 = Pc_2_bar*100
    P = P_bar*100    

    a_1 = Psat.a(T, Tc_1, Pc_1, m_1)
    a_2 = Psat.a(T, Tc_2, Pc_2, m_2)
    
    b_1 = (R*Tc_1)/(8*Pc_1)
    b_2 = (R*Tc_2)/(8*Pc_2)    
    
    
    amix = a_mix(a_1, a_2, x1)
    bmix = b_mix(b_1, b_2, x1)
    
    volumes_l = dV_l(R, P_bar, x1, Tc_1, Pc_1_bar, Tc_2, Pc_2_bar, m_1, m_2)
    
    Term_1 = 
    
def dV_l(T, P_bar, x1, Tc_1, Pc_1_bar, Tc_2, Pc_2_bar, m_1, m_2, Vmixture):
    single_1 = V_single(T, Tc, P_bar, Pc_1_bar, m_1)[0]
    single_2 = V_single(T, Tc, P_bar, Pc_2_bar, m_2)[0]
    x2 = 1 - x1
    sigma = single_1*x1 + single2*x2
    return [Vmixture - sigma, single_1, single_2]


def V_single(T, Tc, P_bar, Pc_bar, m):
    Pc = Pc_bar*100
    P = P_bar*100
    calc = Psat.Psat(T, Tc, Pc_bar, m)
    sat_P = calc[0]*100


    a_eq = Psat.a(T, Tc, Pc, m)
    b = (R*Tc)/(8*Pc)

    def eq(V):
        return abs(Psat.vdw(T, a_eq, b, V, 0) - P)
    
    if P > sat_P:
        V_id_l = sc_o.fsolve(eq, [0.99*calc[1]])
        return [V_id_l, 0]
    elif P < sat_P:
        V_id_v = sc_o.fsolve(eq, [1.01*calc[2]])
        return [0, V_id_v]
    else:
        V_id_l = calc[1]
        V_id_v = calc[2]
        return [V_id_l, V_id_v]

def x2(x1):
    return 1 - x1

def a_mix(a11, a22, x1):
    a12 = (a11 +a22)/2
    x2 = 1 - x1
    return a11*x1**2 + 2*a12*x1*x2 + a22*x2**2

def b_mix(b1, b2, x1):
    x2 = 1 - x1
    return b1*x1 + b2*x2

print V_single(500, 508.1, 50, 47.02, 0.975)[0]
