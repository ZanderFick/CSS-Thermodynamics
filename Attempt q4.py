import wx
from wx.lib.plot import PlotCanvas, PlotGraphics, PolyMarker, PolyLine
import numpy as np
import Psat

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
        
        self.P_sat = wx.StaticText(self.ButtonPanel, label="Psat = ")
        
        
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
        self.panelsizer.Add(self.b_val, (1, 2))
        self.panelsizer.Add(self.P_unit, (3, 1))
        self.panelsizer.Add(self.m_val, (2, 2))
        self.panelsizer.Add(self.P_sat, (4, 2))

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

        self.SetSizer(self.windowSizer)

#Bind Events

        self.Button_calc_and_plot.Bind(wx.EVT_BUTTON, self.calculate_and_plot)
        self.T_slider.Bind(wx.EVT_SLIDER, self.calculate_and_plot)
        
#Initialise Window
        self.Centre()
        self.Show()
        

    def calculate_and_plot(self, event):
        print ''
        print 'go'
        print ''
        self.T_slider_label.SetLabel('Isotherm Temperature = ' + str(self.T_slider.GetValue()) + ' K')
        Tc = 0
        Pc = 0
        m = 0
        T = float(self.T_slider.GetValue())
        T_Text = str(self.T_slider.GetValue())       

        Tc_text = self.Tc_input.GetValue()
        Pc_text = self.Pc_input.GetValue()
        m_text = self.m_input.GetValue()

        if (Tc_text != '') and (Pc_text != '') and (m_text != ''):
            Tc = float(Tc_text)
            Pc = float(Pc_text)*100
            m = float(m_text)

# Bereken aCrit
        if (Tc != 0) and (Pc != 0) and (m != 0):
            R = Psat.R
            ac = (27*R*R*Tc*Tc)/(64*Pc)
            b = (R*Tc)/(8*Pc)
            
            calc = Psat.Psat(T, Tc, Pc/100, m)
            Psatval = calc[0]
            Pguesss = calc[1]

            self.ac_label = str(round(ac,4))
            self.b_label = str(round(b,4))
            if np.isreal(Psatval):
                self.Psat_label = str(round(Psatval,4))
            else:
                self.Psat_label = Psatval

            self.ac_val.SetLabel('ac = ' + self.ac_label)
            self.b_val.SetLabel('b = ' + self.b_label)
            self.m_val.SetLabel('m = ' + m_text)
            self.P_sat.SetLabel('Psat = '+ self.Psat_label)

#Calculate datapoints for vdw EOS 
            Data = np.zeros((1001))
            Line_Data = np.zeros((1001))
            
            for k in range(0,1001):
                Data[k] = 100000*((0.1**(10-0.01*k))) + b+0.01
                Value = max([min([Psat.vdw(T,Psat.a(T, Tc, Pc, m),b , Data[k],0)/100, 100]), -1000])
                Line_Data[k] = Value
            Plot_Data1 = np.vstack((Data, Line_Data)).T
            
            
            Isotherm = PolyLine(Plot_Data1, legend= 'Isotherm, T = ' + T_Text, colour='blue')  
            
            self.canvas.Draw(PlotGraphics([Isotherm], '', ' V', 'P'))
            self.canvas.setLogScale((True, False))


if __name__ == '__main__':

    app = wx.App()
    Psat_finder(None, title='Psat finder')
    app.MainLoop()