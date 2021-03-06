import numpy
import scipy.optimize as sc_o

R = 8.314472

def  Psat(T, Tc, Pcbar, m):
    Pc = Pcbar*100
    b = (R*Tc)/(8*Pc)

    aP = a(T, Tc, Pc, m)
#derivative
    def Vdw_der(V):
        return d_vdw(T, aP, b, V)   

    if T > Tc:
        return 'Temperature above critical point'
    else:
#initial guess
# Find the V roots numerically: Volume Initial guess equations from the literature
        V_extremes_i = numpy.real(d_vdw_roots(T, aP, b))
        Pressures = vdw(T, aP, b, V_extremes_i, 0)
        Psum = sum(Pressures[Pressures>min(Pressures)])/2

        if Psum <0:
            Psat_guess_i = 0.1
        else:
            Psat_guess_i = Psum
#Find real Psat
        def goal_func(Pguess):
            def Vdw_eq_goal(V):
                return vdw(T, aP, b, V, Pguess)

            V_roots = vdw_roots(T, aP, b, Pguess)
            V_roots =  V_roots[V_roots >0]
            V_roots = numpy.real(V_roots[numpy.isreal(V_roots)])
            V_l = min(V_roots)
            V_v = max(V_roots)
           
            return P_opt(V_v, V_l, T, Pguess, aP, b)

        result = sc_o.fsolve(goal_func,Psat_guess_i)
        
        V_roots = vdw_roots(T, aP, b, result)
        V_roots =  V_roots[V_roots >0]
        V_roots = numpy.real(V_roots[numpy.isreal(V_roots)])
        V_l_ret = min(V_roots)
        V_v_ret = max(V_roots)
        
        return [result/100, V_l_ret, V_v_ret]

# Van der Waals EOS
def vdw(T, a, b, V, Pguess):
    term1 = R*T/(V-b)
    term2 = -a/(V**2)
    P = term1+term2 
    return P - Pguess

def vdw_roots(T, a, b, Pguess):
    C1 = Pguess
    C2 = -R*T - b*Pguess
    C3 = a 
    C4 = -a*b
    C_array = [C1, C2, C3, C4]
    return numpy.roots(C_array)

def d_vdw_roots(T, a, b):
    C1 = -R*T
    C2 = 2*a
    C3 = -4*a*b
    C4 = 2*a*b*b
    C_array = [C1, C2, C3, C4]
    return numpy.roots(C_array)


#Function to Optomise:
def P_opt(Vv, Vl, T, Psat, a, b):
    return R*T*numpy.log((Vv-b)/(Vl-b)) + Psat*(Vl - Vv) + a*((1/Vv)-(1/Vl))

# a parameter function
def a(T, Tc, Pc, m):
    ac = (27*R*R*Tc*Tc)/(64*Pc)
    Tr = T/Tc
    exponent = m*(1-Tr)
    a = ac*(numpy.e)**exponent
    return a