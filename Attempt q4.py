import wx
from wx.lib.plot import PlotCanvas, PlotGraphics, PolyMarker, PolyLine
import numpy as np


class Psat_finder(wx.Frame):

    def __init__(self, parent, title):
        super(Psat_finder, self).__init__(parent, title=title, size=(800, 700))
        self.ButtonPanel = wx.Panel(self)


#Main window

        self.windowSizer = wx.BoxSizer(wx.VERTICAL)

        self.windowSizer.Add(self.ButtonPanel, 1, wx.EXPAND | wx.ALL , 1)

# Input boxes

        self.ac_val = wx.StaticText(self.ButtonPanel, label="ac = ")

        self.b_val = wx.StaticText(self.ButtonPanel, label="b = ")

        self.m_val = wx.StaticText(self.ButtonPanel, label="m = ")
        
        
        self.T_unit = wx.StaticText(self.ButtonPanel, label="Kelvin ")
        
        self.P_unit = wx.StaticText(self.ButtonPanel, label="Bar")

        self.Tc_input_label = wx.StaticText(self.ButtonPanel, label="Tc")
        self.Tc_input = wx.TextCtrl(self.ButtonPanel, size=(70, -1))

        self.Pc_input_label = wx.StaticText(self.ButtonPanel, label="Pc")
        self.Pc_input = wx.TextCtrl(self.ButtonPanel, size=(70, -1))

        self.m_input_label = wx.StaticText(self.ButtonPanel, label="m")
        self.m_input = wx.TextCtrl(self.ButtonPanel, size=(70, -1))

        self.T_slider = wx.Slider(self.ButtonPanel, -1, 273, 0, 1500, wx.DefaultPosition, (250,-1))
        self.T_slider_label = wx.StaticText(self.ButtonPanel, label="Isotherm Temperature = 273 K")

        self.Button_calc_and_plot = wx.Button(self.ButtonPanel, label='Calculate Parameters and Plot')

        self.PlotPanel = wx.Panel(self.ButtonPanel)
        
        self.canvas = PlotCanvas(self.PlotPanel)
        
        self.PlotPanelSizer = wx.BoxSizer()

        self.panelsizer = wx.GridBagSizer(9, 3)

        self.panelsizer.Add(self.Tc_input_label, (0, 0))
        self.panelsizer.Add(self.Tc_input, (1, 0))

        self.panelsizer.Add(self.Pc_input_label, (2, 0))
        self.panelsizer.Add(self.Pc_input, (3, 0))

        self.panelsizer.Add(self.m_input_label, (4, 0))
        self.panelsizer.Add(self.m_input, (5, 0))


        self.panelsizer.Add(self.ac_val, (0, 2))
        self.panelsizer.Add(self.T_unit, (1, 1))
        self.panelsizer.Add(self.b_val, (2, 2))
        self.panelsizer.Add(self.P_unit, (3, 1))
        self.panelsizer.Add(self.m_val, (4, 2))

        self.panelsizer.Add(self.Button_calc_and_plot, (5, 1))
        
        self.panelsizer.Add(self.T_slider_label, (6, 0), (1, 2), wx.EXPAND)
        self.panelsizer.Add(self.T_slider, (7, 0), (1, 3), wx.EXPAND)

        self.PlotPanelSizer.Add(self.canvas, 1, wx.EXPAND)

        self.panelsizer.Add(self.PlotPanel,(9, 0), (1, 3), wx.EXPAND)

        self.panelsizer.AddGrowableCol(2)
        self.panelsizer.AddGrowableRow(9)

#initialize sizers
        self.PlotPanel.SetSizerAndFit(self.PlotPanelSizer)

        self.ButtonPanel.SetSizerAndFit(self.panelsizer)

        self.SetSizerAndFit(self.windowSizer)

#Bind Events

        self.Button_calc_and_plot.Bind(wx.EVT_BUTTON, self.calculate_and_plot)
        self.T_slider.Bind(wx.EVT_SLIDER, self.calculate_and_plot)
        
#Initialise Window
        self.Centre()
        self.Show()

    def calculate_and_plot(self, event):
        R = 8.314472
        Tc = 0
        Pc = 0
        m = 0  
        T = float(self.T_slider.GetValue())
        T_Text = str(self.T_slider.GetValue())

        self.T_slider_label.SetLabel('Isotherm Temperature = ' + T_Text + ' K')

        Tc_text = self.Tc_input.GetValue()
        Pc_text = self.Pc_input.GetValue()
        m_text = self.m_input.GetValue()

        if (Tc_text != '') and (Pc_text != '') and (m_text != ''):
            Tc = float(Tc_text)
            Pc = float(Pc_text)*101.325
            m = float(m_text)

# Bereken aCrit
        if (Tc != 0) and (Pc != 0) and (m != 0):
            ac = (27*R*R*Tc*Tc)/(64*Pc)
            b = (R*Tc)/(8*Pc)

            self.ac_label = str(round(ac,4))
            self.b_label = str(round(b,4))

            self.ac_val.SetLabel('ac = ' + self.ac_label)
            self.b_val.SetLabel('b = ' + self.b_label)
            self.m_val.SetLabel('m = ' + m_text)

#Calculate datapoints for vdw EOS 
            a = self.a(ac, m, T, Tc)
            
            Data = np.zeros((1001))
            
            for k in range(0,1001):
                Data[k] = 1000*((0.1**(10-0.01*k))) + b+0.01
            Line_Data = self.vdw(T, Data, a, b)/(101.325)
            Plot_Data1 = np.vstack((Data, Line_Data)).T
            
            Isotherm = PolyLine(Plot_Data1, legend= 'Isotherm, T = ' + T_Text, colour='blue')     
               
            
            if T < Tc:
                maxminV = self.maxmin(Data, T, Tc, a, b) 
                Lmin = self.vdw(T, float(maxminV[0]), a, b)/(101.325)
                Rmax = self.vdw(T, float(maxminV[1]), a, b)/(101.325)
                
                Psatguess = (Lmin + Rmax)/2
                
                GuessV = self.roots(T, Data, Plot_Data1, Psatguess, a, b)
                
                Marker1 = self.vdw(T, float(GuessV[0]), a, b)/(101.325)
                Marker2 = self.vdw(T, float(GuessV[1]), a, b)/(101.325)
                Marker3 = self.vdw(T, float(GuessV[2]), a, b)/(101.325)  
                
                Markers = np.hstack((Marker1, Marker2, Marker3))
                
                Markerdata = np.vstack((GuessV, Markers)).T
                print Markerdata
                
                Markerplot = PolyMarker(Markerdata, colour= 'green')
                
                 
                
                self.canvas.Draw(PlotGraphics([Isotherm, Markerplot], '', ' V', 'P'))
            else:
                self.canvas.Draw(PlotGraphics([Isotherm], '', ' V', 'P'))    
            self.canvas.setLogScale((True, False))
            
            

# a parameter function
    def a(self, ac, m, T, Tc):
        Tr = T/Tc
        exponent = m*(1-Tr)
        a = ac*(np.e)**exponent
        return a

# Van der Waals EOS
    def vdw(self, T, V, a, b):
        R = 8.314472
        term1 = R*T/(V-b)
        term2 = -a/(V**2)
        P = term1+term2 
        return P 
        
    def d_vdw_dV(self, T, V, a, b):
        R = 8.314472
        term1 = -R*T/((V-b)**2)
        term2 = 2*a/(V**3)
        dPdV = term1+term2 
        return dPdV

    def maxmin(self, Data, T, Tc, a, b):
#Find the minimum and maximum point
        minV = min(Data)
        maxV = max(Data)
        
        if T < Tc:
        # Left side derivative change
            check = 1
            d = self.d_vdw_dV(T, minV, a, b)
            sign = abs(d)/d
            i = 0
            while (check == 1) and (i < 1000):
                d = self.d_vdw_dV(T, Data[i], a, b)
                newsign = abs(d)/d
                check = sign/newsign
                sign = newsign
                i += 1
            L = Data[i-1]
        # Right side derivative change
            check = 1
            d = self.d_vdw_dV(T, maxV, a, b)
            sign = abs(d)/d
            i = 1000
            while (check == 1) and (i > 0):
                d = self.d_vdw_dV(T, Data[i], a, b)
                newsign = abs(d)/d
                check = sign/newsign
                sign = newsign
                i += -1
            R = Data[i+1]
            
            return [L, R]


    def roots(self,T, Vdata, Pdata, Psatguess, a , b):
        minV = min(Vdata)
        maxV = max(Vdata)
        newdata = Pdata - Psatguess

        check = 1
        d = self.vdw(T, minV, a, b)
        sign = abs(d)/d
        i = 0
        while (check == 1) and (i < 1000):
            d = self.vdw(T, Vdata[i], a, b)
            newsign = abs(d)/d
            check = sign/newsign
            sign = newsign
            i += 1
        R1 = Vdata[i-1]

        check = 1
        d = self.vdw(T, Vdata[i], a, b)
        sign = abs(d)/d
        while (check == 1) and (i < 1000):
            d = self.vdw(T, Vdata[i], a, b)
            newsign = abs(d)/d
            check = sign/newsign
            sign = newsign
            i += 1
        R2 = Vdata[i-1]

        check = 1
        d = self.vdw(T, Vdata[i], a, b)
        sign = abs(d)/d
        while (check == 1) and (i < 1000):
            d = self.vdw(T, Vdata[i], a, b)
            newsign = abs(d)/d
            check = sign/newsign
            sign = newsign
            i += 1
        R3 = Vdata[i-1]
        
        return [R1, R2, R3]

if __name__ == '__main__':

    app = wx.App()
    Psat_finder(None, title='Psat finder')
    app.MainLoop()