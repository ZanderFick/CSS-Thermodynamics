import numpy
import scipy.optimize as sc_o
import Psat

R = 8.314472

def dG_RT(T, P_bar, x1, Tc_1, Pc_bar_1, m_1, Tc_2, Pc_bar_2, m_2):
    Pc_1 = Pc_bar_1*100
    Pc_2 = Pc_bar_2*100
    P = P_bar*100

    x2 = 1- x1

    b1 = (R*Tc_1)/(8*Pc_1)
    b2 = (R*Tc_2)/(8*Pc_2)
    
    bmix = b_mix(x1, x2, b1, b2)
    
    a1 = Psat.a(T, Tc_1, Pc_1, m_1)
    a2 = Psat.a(T, Tc_2, Pc_2, m_2)

    amix = a_mix(x1, x2, a1, a2)

    Volumes_1 = Volume_solve(T, P_bar, Pc_bar_1, Tc_1, m_1)
    Volumes_2 = Volume_solve(T, P_bar, Pc_bar_2, Tc_2, m_2)

#Liquid phase:
    V_l_1 = Volumes_1[0]
    V_l_2 = Volumes_2[0]
    print bmix
    
    if V_l_1 != 0 and V_l_2 != 0:
    
        Vmix_l = V_mix(x1, x2, V_l_1, V_l_2)
        deltaV_l = delta_V(x1, x2, V_l_1, V_l_2)
    
        term_1_l = P*deltaV_l/(R*T)

        term_2_l_a = numpy.log((V_l_1 - b1)/(Vmix_l - bmix))
        term_2_l_b = numpy.log((V_l_2 - b2)/(Vmix_l - bmix))
        term_2_l = x1*term_2_l_a + x2*term_2_l_b
    
        term_3_l_a = a1/(R*T*V_l_1)
        term_3_l_b = a2/(R*T*V_l_2)
        term_3_l = x1*term_3_l_a + x2*term_3_l_b
        
        term_4_l = amix/(R*T*Vmix_l)
    
        dG_l = term_1_l + term_2_l + term_3_l - term_4_l
        return dG_l
#Vapour Phase:
    V_v_1 = Volumes_1[1]
    V_v_2 = Volumes_2[1]

    if V_v_1 != 0 and V_v_2 != 0:
    
        Vmix_v = V_mix(x1, x2, V_v_1, V_v_2)
        deltaV_v = delta_V(x1, x2, V_v_1, V_v_2)
    
        term_1_v = P*deltaV_v/(R*T)
        term_2_v_a = numpy.log((V_v_1 - b1)/(Vmix_v - bmix))
        term_2_v_b = numpy.log((V_v_2 - b2)/(Vmix_v - bmix))
        term_2_v = x1*term_2_v_a + x2*term_2_v_b
    
        term_3_v_a = a1/(R*T*V_v_1)
        term_3_v_b = a2/(R*T*V_v_2)
        term_3_v = x1*term_3_v_a + x2*term_3_v_b
        
        term_4_v = amix/(R*T*Vmix_v)
        
    
        dG_v = term_1_v + term_2_v + term_3_v - term_4_v
        return dG_v

    

def Volume_solve(T, P_bar, Pc_bar, Tc, m):
    Pc = Pc_bar*100
    P = P_bar*100
    
    calc = Psat.Psat(T, Tc, Pc_bar, m)
    
    P_sat_bar = calc[0]
    
    P_sat = P_sat_bar*100
    a_eq = Psat.a(T, Tc, Pc, m)
    b_eq = (R*Tc)/(8*Pc)

    def V_root(V):
        return Psat.vdw(T, a_eq, b_eq, V, 0)-P
    if P > P_sat:
        V_guess = calc[1]*0.95
        V_root = sc_o.fsolve(V_root,V_guess)
        return [V_root[0], 0]
    elif P < P_sat:
        V_guess = calc[2]*1.05
        V_root = sc_o.fsolve(V_root,V_guess)
        return [0, V_root[0]]
    else:
        return[calc[1], calc[2]]

def V_mix(x1, x2, V1, V2):
    return V1*x1 + x2*V2

def delta_V(x1, x2, V1, V2):
    term_2_a = x1*V1
    term_2_b = x2*V2
    
    return V_mix(x1, x2, V1, V2) -(term_2_a + term_2_b)

def b_mix(x1, x2, b1, b2):
    return b1*x1 + b2*x2

def a_mix(x1, x2, a11, a22):
    a12 = numpy.sqrt(a11*a22)
    return a11*(x1**2) + 2*a12*x1*x2 + a22*(x2**2)
        
import pylab
Data_x = numpy.linspace(0,1)
Data_y = numpy.zeros(Data_x.size)
for k in range(Data_x.size):
    Data_y[k] = dG_RT(450, 15, Data_x[k], 507.9, 30.35, 0.969, 647.3, 220.64, 1.014)
pylab.plot(Data_x, Data_y, 'r')
pylab.show()