#Returns the Vl Vv ac and m optimised parameters
import csv
import numpy
from lmfit import minimize, Parameters
import Psat

R = 8.314472

def residual(Tc, Pc, Tb, Data):
    Data_x = Data[:,0]
    Data_y = Data[:,1]
    
    params = Parameters()
    
    params.add('T', value = Tb*1.1, min = Tb*1.05, max = 0.99*Tc)
    params.add('m', value = 0.9, min = 0.1)
    params.add('Tc', value = Tc, vary=False)
    params.add('Pc', value = Pc, vary=False)
    
    out = minimize(vdw_residual, params, args=(Data_x, Data_y))
    res = Data_y + out.residual
    import pylab
    pylab.plot(Data_x, Data_y, 'k+')
    pylab.plot(Data_x, res, 'r')
    pylab.show()

    return out.params['m'].value


def vdw_residual(params, x, data):
    T = params['T'].value
    m = params['m'].value
    Tc = params['Tc'].value
    Pc = params['Pc'].value    
    
    x_range = numpy.size(x, axis=0)
    model_data = numpy.zeros(x_range)
    for k in range(x_range):
        calc = Psat.Psat(x[k], Tc, Pc, m)
        model_data[k] = calc[0]*100
    return (model_data-data)/abs(model_data)
