import numpy

import Psat

R = 8.314472


#def G_RT(T, P, Vmix, x1):
    
def dV(T, P, x1, Tc_1, Pc_1, Tc_2, Pc_2, m_1, m_2):

    sat_P_1 = Psat.Psat(T, Tc_1, Pc_1, m_1)
    a_1 = Psat.a(T, Tc_1, Pc_1, m_1)
    b_1 = (R*Tc_1)/(8*Pc_1)

    def eq_1(V):
        return Psat.vdw(T, a_1, b_1, V, 0)
    
    if P != sat_P_1:
        eq = Psat.vdw()

    sat_P_2 = Psat.Psat(T, Tc_2, Pc_2, m_2)
    a_2 = Psat.a(T, Tc_2, Pc_2, m_2)
    b_2 = (R*Tc_2)/(8*Pc_2)

def x2(x1):
    return 1 - x1