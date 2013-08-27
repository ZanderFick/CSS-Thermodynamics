#Returns the Vl Vv ac and m optimised parameters
import csv
import numpy
from lmfit import minimize, Parameters
import Psat

R = numpy.float64(8.314472)

def residual(Tc, Pc, Tb, Data):
    import_data()
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


Test_data_acetone = numpy.array(   [[259.0,	4.220],
                        [267.9,	7.034],
                        [277.1,	11.459],
                        [286.5,	18.264],
                        [296.4,	28.514],
                        [306.5,	43.647],
                        [317.0,	65.574],
                        [327.9,	96.792],
                        [339.1,	140.505],
                        [350.7,	200.770],
                        [362.8,	282.663],
                        [375.2,	392.460],
                        [388.1,	537.855],
                        [401.3,	728.219],
                        [415.1,	974.919],
                        [429.3,	1291.736],
                        [444.0,	1695.453],
                        [459.3,	2206.747],
                        [475.0,	2851.724],
                        [491.3,	3665.109],
                        [506.4,	4586.851],
                        [508.1,	4702.000]])


print residual(508.1, 47.02, 329.22, Test_data_acetone)