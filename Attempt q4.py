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

            Vc = (8./9)*ac/(R*Tc)   
            Data = np.arange(b+0.1,1000.1+b,0.1)
            Line_Data = (self.vdw(T, Data, self.a(ac, m, T, Tc), b))/(101.325)
            Plot_Data = np.vstack((Data, Line_Data)).T
            
            Isotherm = PolyLine(Plot_Data, legend= 'Isotherm, T = ' + T_Text, colour='red')
            
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
        term2 = a/(V**2)
        P = term1-term2 
        return P        
        

if __name__ == '__main__':

    app = wx.App()
    Psat_finder(None, title='Psat finder')
    app.MainLoop()