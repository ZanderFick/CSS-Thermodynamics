import numpy
import scipy.optimize as sc_o
import Psat

R = 8.314472


def G_RT_main(T, P_bar, x1,Tc_1, Pc_1_bar, m_1, Tc_2, Pc_2_bar, m_2):
    Pc_1 = Pc_1_bar*100
    Pc_2 = Pc_2_bar*100
    P = P_bar*100    
    x2 = 1-x1

    a_1 = Psat.a(T, Tc_1, Pc_1, m_1)
    a_2 = Psat.a(T, Tc_2, Pc_2, m_2)
    
    b_1 = (R*Tc_1)/(8*Pc_1)
    b_2 = (R*Tc_2)/(8*Pc_2)    
    
    amix = a_mix(a_1, a_2, x1)
    bmix = b_mix(b_1, b_2, x1)
    
    
    volumes = dV(T, P_bar, x1, Tc_1, Pc_1_bar, Tc_2, Pc_2_bar, m_1, m_2)
   
    delta_V = volumes[0]
    V_i = [volumes[1], volumes[2]]
    Vmix = volumes[3]
    
    term_1 = P*delta_V[0]/R*T

    term_2_1 = numpy.log((V_i[0]-b_1)/(Vmix-bmix))
    term_2_2 = numpy.log((V_i[1]-b_2)/(Vmix-bmix))
    term_2 = x1*term_2_1 + x2*term_2_2
    term_3_1 = a_1/(R*T*V_i[0])
    term_3_2 = a_2/(R*T*V_i[1])
    
    term_3 = x1*term_3_1 + x2*term_3_2
    
    term_4 = -amix/(R*T*Vmix)
    
    terms = [term_1, term_2, term_3, term_4]
    return sum(terms)

    
def dV(T, P_bar, x1, Tc_1, Pc_1_bar, Tc_2, Pc_2_bar, m_1, m_2):
    single_1 = V_single(T, Tc_1, P_bar, Pc_1_bar, m_1)[0]
    single_2 = V_single(T, Tc_2, P_bar, Pc_2_bar, m_2)[0]
    Vmixture =  V_mixture(x1, single_1, single_2)
    x2 = 1 - x1
    sigma = single_1*x1 + single_2*x2
    return [Vmixture - sigma, single_1, single_2, Vmixture]

def V_mixture(x1, V_1, V_2):
    x2 = 1-x1
    return x1*V_1 + x2*V_2


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
        return [V_id_l]
    elif P < sat_P:
        V_id_v = sc_o.fsolve(eq, [1.01*calc[2]])
        return [V_id_v]
    else:
        V_id_l = calc[1]
        V_id_v = calc[2]
        return [V_id_l, V_id_v]

def x2(x1):
    return 1 - x1

def a_mix(a11, a22, x1):
    a12 = numpy.sqrt(a11*a22)
    x2 = 1 - x1
    return a11*x1**2 + 2*a12*x1*x2 + a22*x2**2

def b_mix(b1, b2, x1):
    x2 = 1 - x1
    return b1*x1 + b2*x2

#           G_RT_main(T, P_bar, x1,Tc_1, Pc_1_bar, Tc_2, Pc_2_bar, m_1, m_2)


import pylab
Data_x = numpy.linspace(0,1)
Data_y_1 = G_RT_main(500, 20, Data_x, 508.1, 47.02, 0.967, 647.3, 220.64, 1.014)
Data_y_2=  G_RT_main(300, 20, Data_x, 508.1, 47.02, 0.967, 647.3, 220.64, 1.014)
pylab.plot(Data_x, Data_y_1, 'r')
pylab.plot(Data_x, Data_y_2, 'b')
pylab.show()


