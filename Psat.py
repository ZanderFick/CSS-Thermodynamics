import numpy
import scipy.optimize as sc_o

R = 8.314472

def  Psat(T, Tc, Pcbar, m):
    Pc = Pcbar*100

    ac = (27*R*R*Tc*Tc)/(64*Pc)
    b = (R*Tc)/(8*Pc)
    Vc = (8./9)*ac/(R*Tc) 

    aP = a(T, Tc, ac, m)
#derivative
    def Vdw_der(V):
        return d_vdw(T, aP, b, V)   
#initial guess
    V_extremes_i = sc_o.fsolve(Vdw_der,[0.1*R*T/Pc+0.1, R*T/Pc +1])  
    Psat_guess_i = vdw(T, aP, b, V_extremes_i[1], 0)

    def Vdw_eq(V):        
        return vdw(T, aP, b, V, Psat_guess_i) 

#Find real Psat
    def goal_func(Pguess):
        def Vdw_eq_goal(V):
           return vdw(T, aP, b, V, Pguess)
        V_extremes = sc_o.fsolve(Vdw_eq_goal,[0.1*R*T+0.1, R*T +1])  
        V_l = V_extremes[0]
        V_v = V_extremes[1]
        print Pguess
        return P_opt(V_v, V_l, T, Pguess, aP,  b)

    ao = sc_o.minimize(goal_func,[Psat_guess_i])

# Van der Waals EOS
def vdw(T, a, b, V, Pguess):
    term1 = R*T/(V-b)
    term2 = -a/(V**2)
    P = term1+term2 
    return P - Pguess

def d_vdw(T, a, b, V):
    term1 = -R*T/((V-b)**2)
    term2 = 2*a/(V**3)
    P = term1+term2 
    return P/100

#Function to Optomise:
def P_opt(Vv, Vl, T, Psat, a, b):
    return R*T*numpy.log((Vv-b)/(Vl-b)) + Psat*(Vl - Vv) + a*((1/Vv)-(1/Vl))

# a parameter function
def a(T, Tc, ac, m):
    Tr = T/Tc
    exponent = m*(1-Tr)
    a = ac*(numpy.e)**exponent
    return a

print Psat(306, 508.1, 47.02, 0.978)