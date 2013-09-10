import numpy
import scipy.misc as sc_m
import scipy.optimize as sc_o
from lmfit import minimize, Parameters
import Psat
import pylab

R = 8.314472

def dG_RT(T, P_bar, x1, Tc_1, Pc_bar_1, m_1, Tc_2, Pc_bar_2, m_2, Go1, Go2):
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
    
    if x1 != 0 and x1 != 1:
        Gi = (x1*numpy.log(x1) + x2*numpy.log(x2))+(Go1*x1 + Go2*x2)/(R*T)
    else:
        Gi = 0
    
    if V_l_1 != 0 or V_l_2 != 0:
    
        Vmix_l = lin_mix(x1, x2, V_l_1, V_l_2)
        deltaV_l = delta_V(x1, x2, V_l_1, V_l_2)    
        dG_l = P*deltaV_l/(R*T) + x1*(numpy.log((V_l_1 - b1)/(Vmix_l - bmix))) + x2*(numpy.log((V_l_2 - b2)/(Vmix_l - bmix))) + x1*(a1/(R*T*V_l_1)) + x2*(a2/(R*T*V_l_2)) - amix/(R*T*Vmix_l)
        return dG_l+ Gi
#Vapour Phase:
    V_v_1 = Volumes_1[1]
    V_v_2 = Volumes_2[1]

    if V_v_1 != 0 or V_v_2 != 0:
    
        Vmix_v = lin_mix(x1, x2, V_v_1, V_v_2)
        deltaV_v = delta_V(x1, x2, V_v_1, V_v_2)
        dG_v = P*deltaV_v/(R*T) + x1*(numpy.log((V_v_1 - b1)/(Vmix_v - bmix))) + x2*(numpy.log((V_v_2 - b2)/(Vmix_v - bmix))) + x1*(a1/(R*T*V_v_1)) + x2*(a2/(R*T*V_v_2)) - amix/(R*T*Vmix_v)
        return  dG_v+ Gi

def deriv_dG_RT(T, P_bar, x1, Tc_1, Pc_bar_1, m_1, Tc_2, Pc_bar_2, m_2, Go1, Go2):
    def diff_func(x):
        return dG_RT(T, P_bar, x, Tc_1, Pc_bar_1, m_1, Tc_2, Pc_bar_2, m_2, Go1, Go2)
    return sc_m.derivative(diff_func, x1, dx=1e-5) 

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

def lin_mix(x1, x2, A1, A2):
    return A1*x1 + x2*A2

def delta_V(x1, x2, V1, V2):
    term_2_a = x1*V1
    term_2_b = x2*V2
    
    return lin_mix(x1, x2, V1, V2) -(term_2_a + term_2_b)

def b_mix(x1, x2, b1, b2):
    return b1*x1 + b2*x2

def a_mix(x1, x2, a11, a22):
    a12 = numpy.sqrt(a11*a22)
    print a11*(x1**2) + 2*a12*x1*x2 + a22*(x2**2)
    return a11*(x1**2) + 2*a12*x1*x2 + a22*(x2**2)

def tangent(T, P_bar, Tc_1, Pc_bar_1, m_1, Tc_2, Pc_bar_2, m_2, Go1, Go2):
    def Gmix(x):
        return dG_RT(T, P_bar, x, Tc_1, Pc_bar_1, m_1, Tc_2, Pc_bar_2, m_2, Go1, Go2)
    def d_Gmix(x):
        return deriv_dG_RT(T, P_bar, x, Tc_1, Pc_bar_1, m_1, Tc_2, Pc_bar_2, m_2, Go1, Go2)

    min_root = sc_o.fsolve(d_Gmix, 0.01)
    max_root = sc_o.fsolve(d_Gmix, 0.98)
    bnds = ((0.001, max_root[0]*0.8),(min_root[0]*1.1, 0.999))
    init = [min_root[0]*1.1, 0.999]
    
    def resid(x_vect):
        xa = x_vect[0]
        xb = x_vect[1]
        m1 = d_Gmix(xa)
        m2 = d_Gmix(xb)

        c1 = Gmix(xa) - xa*m1
        c2 = Gmix(xb) - xb*m2
        e1 = m1 + c1
        e2 = m2 + c2
        error = max([abs(e1-e2),abs(c1-c2)])
        return error

    out = sc_o.fmin(resid,init )
    return out


