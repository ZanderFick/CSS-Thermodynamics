import numpy
import scipy.optimize as sc_o
import Psat

R = 8.314472


#def G_RT(T, P, Vmix, x1):
    
def dV(T, P, x1, Tc_1, Pc_1, Tc_2, Pc_2, m_1, m_2, Vmixture):



def V_ideal(T, Tc, P, Pc, m)
    calc = Psat.Psat(T, Tc, Pc, m)
    sat_P = calc[0]

    a_eq = Psat.a(T, Tc, Pc, m)
    b = (R*Tc)/(8*Pc)

    def eq(V):
        return Psat.vdw(T, a_eq, b, V, 0)
    
    if P > sat_P:
        V_id = sc_o.fsolve(eq, [0.99*calc[1]])
    elif P < sat_P:
        V_id = sc_o.fsolve(eq, [1.01*calc[2]])
    else:
        V_id = calc[1] + calc[2]
    return V_id

def x2(x1):
    return 1 - x1

print V_ideal(400, 508.1, 5)