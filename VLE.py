import numpy
import scipy.optimize as sc_o
import Psat

R = 8.314472


#def G_RT(T, P, Vmix, x1):
    
#def dV(T, P, x1, Tc_1, Pc_1, Tc_2, Pc_2, m_1, m_2, Vmixture):



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
        V_id = sc_o.fsolve(eq, [0.99*calc[1]])
    elif P < sat_P:
        V_id = sc_o.fsolve(eq, [1.01*calc[2]])
    else:
        V_id = calc[1] + calc[2]
    return V_id

def x2(x1):
    return 1 - x1

print V_single(500, 508.1, 50, 47.02, 0.975)