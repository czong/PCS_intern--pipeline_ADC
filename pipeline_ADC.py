import pdb
import sys
import numpy as np
from scipy.signal import *
from traits.api import *
from traitsui.api import *
from matplotlib.figure import Figure
import wx
from traitsui.wx.basic_editor_factory import BasicEditorFactory
import decimal
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
import matplotlib.pyplot as plt
from enable.api import ComponentEditor
from traitsui.wx.editor import Editor
from numpy.fft import fft,fftfreq,fftshift 
import matplotlib.pyplot as plt

class _MPLFigureEditor(Editor):
    scrollable = True

    def init(self,parent):
        self.control = self._create_canvas(parent)
        self.set_tooltip()

    def update_editor(self):
        pass

    def _create_canvas(self,parent):
        """ Create the MPL canvas."""
        # the panel lets us add additional controls.
        panel = wx.Panel(parent,-1,style=wx.CLIP_CHILDREN)
        sizer = wx.BoxSizer(wx.VERTICAL)
        panel.SetSizer(sizer)
        # matplotlib commands to create a canvas
        mpl_control = FigureCanvas(panel,-1,self.value)
        sizer.Add(mpl_control,1,wx.LEFT |wx.TOP |wx.GROW)
        toolbar = NavigationToolbar2Wx(mpl_control)
        sizer.Add(toolbar,0,wx.EXPAND)
        self.value.canvas.SetMinSize((10,10))
        return panel

class MPLFigureEditor(BasicEditorFactory):
    klass = _MPLFigureEditor    


class pipelineADC(HasTraits):
    # variables
    Fs = Expression(400e6)
    Fin = Expression(50.1e6)
    N = Expression(1e4)
     
    # constants
    G = 3.98 
    precision = 3    

    # view panel
    figure = plt.figure()
    traits_view = View(
            VGroup(
                Item('Fs'),
                Item('Fin'),
                Item('N'),
                Item('figure',editor = MPLFigureEditor(),show_label=False)
            ),
            width = 1200,
            height = 900,
            resizable = True,
            title = 'pipeline ADC'
    )
     
    def __init__(self,*args,**kw):
        super(pipelineADC,self).__init__(*args,**kw)
        # eval the expression variables
        self.Fs = self.toEngFormat(self.Fs)
        self.Fin = self.toEngFormat(self.Fin)
        self.N = self.toEngFormat(self.N)
        # dependent values 
        self.Ts = 1./float(self.Fs)
        self.n = np.arange(0,int(float(self.N)))
        self.t = np.asarray(self.n)*self.Ts
        self.on_trait_change(self.update,'Fs')
        self.on_trait_change(self.update,'Fin')
        self.on_trait_change(self.update,'N')    
        self.calculate(self.Fs,self.Fin,self.N)
        #self.axes1 = self.figure.add_axes()
        #self.axes1.plot(self.f,self.Yf,color='blue')
        #self.axes1.set_title('power spectrum')
        self.axes1=self.figure.add_subplot(111)
        self.axes1.plot(self.f,self.Yf,color='blue')
        self.axes1.hold(False)
        self.figure.suptitle('Power Spectrum')
    def calculate(self,Fs,Fin,N):
        #pdb.set_trace()
        N = int(float(N))
        Fin = float(Fin)
        Fs = float(Fs)
        Ts = 1./Fs
        n = np.arange(0,N)
        t = np.asarray(n)*Ts
        
        x0=np.zeros(N)
        x1=np.zeros(N)
        x2=np.zeros(N)
        x3=np.zeros(N)
        x4=np.zeros(N)
        x5=np.zeros(N)
        x6=np.zeros(N)
        x7=np.zeros(N)
        x8=np.zeros(N)
        x9=np.zeros(N)        
        x10=np.zeros(N)
        x11=np.zeros(N)
        OUT = np.zeros(N)
        
        d1=np.zeros(N)
        d2=np.zeros(N)
        d3=np.zeros(N)
        d4=np.zeros(N)
        d5=np.zeros(N)
        d6=np.zeros(N)
        d7=np.zeros(N)
        d8=np.zeros(N)
        d9=np.zeros(N)        
        d10=np.zeros(N)
        d11=np.zeros(N)

        z1=np.zeros(N)
        z2=np.zeros(N)
        z3=np.zeros(N)
        z4=np.zeros(N)
        z5=np.zeros(N)
        z6=np.zeros(N)
        z7=np.zeros(N)
        z8=np.zeros(N)
        z9=np.zeros(N)        
        z10=np.zeros(N)
        z11=np.zeros(N)
        #pdb.set_trace() 
        y = 1*np.sin(2*np.pi*Ts*n*Fin)
        x0=y
        
        for i in range(0,N-4):
            d1[i]=-6./8.
            z1[i]=0
            if x0[i]>-5./8.:
                d1[i]=-4./8.
                z1[i]=1
            
            if x0[i]>-3./8.:
                d1[i]=-2./8.
                z1[i]=2
            
            if x0[i]>-1./8.:
                d1[i]=0
                z1[i]=3

            if x0[i]>1./8.:
                d1[i]=2./8.
                z1[i]=4

            if x0[i]>3./8.:
                d1[i]=4./8.
                z1[i]=5
        
            if x0[i]>5./8.:
                d1[i]=6./8.
                z1[i]=6

            x1[i]=self.G*(x0[i]-d1[i])

            d2[i]=-6./8.
            z2[i]=0
           
            if x1[i]>-5./8.:
                d2[i]=-4./8.
                z2[i]=1
            
            if x1[i]>-3./8.:
                d2[i]=-2./8.
                z2[i]=2
            
            if x1[i]>-1./8.:
                d2[i]=0
                z2[i]=3

            if x1[i]>1./8.:
                d2[i]=2./8.
                z2[i]=4

            if x1[i]>3./8.:
                d2[i]=4./8.
                z2[i]=5
        
            if x1[i]>5./8.:
                d2[i]=6./8.
                z2[i]=6

            x2[i]=self.G*(x1[i]-d2[i])+.001
 
            d3[i]=-6./8.
            z3[i]=0
           
            if x2[i]>-5./8.:
                d3[i]=-4./8.
                z3[i]=1
            
            if x2[i]>-3./8.:
                d3[i]=-2./8.
                z3[i]=2
            
            if x2[i]>-1./8.:
                d3[i]=0
                z3[i]=3

            if x2[i]>1./8.:
                d3[i]=2./8.
                z3[i]=4

            if x2[i]>3./8.:
                d3[i]=4./8.
                z3[i]=5
        
            if x2[i]>5./8.:
                d3[i]=6./8.
                z3[i]=6

            x3[i]=self.G*(x2[i]-d3[i])-.003

            d4[i]=-6./8.
            z4[i]=0
           
            if x3[i]>-5./8.:
                d4[i]=-4./8.
                z4[i]=1
            
            if x3[i]>-3./8.:
                d4[i]=-2./8.
                z4[i]=2
            
            if x3[i]>-1./8.:
                d4[i]=0
                z4[i]=3

            if x3[i]>1./8.:
                d4[i]=2./8.
                z4[i]=4

            if x3[i]>3./8.:
                d4[i]=4./8.
                z4[i]=5
        
            if x3[i]>5./8.:
                d4[i]=6./8.
                z4[i]=6

            x4[i]=4*(x3[i]-d4[i])+.01

            d5[i]=-6./8.
            z5[i]=0
           
            if x4[i]>-5./8.:
                d5[i]=-4./8.
                z5[i]=1
            
            if x4[i]>-3./8.:
                d5[i]=-2./8.
                z5[i]=2
            
            if x4[i]>-1./8.:
                d5[i]=0
                z5[i]=3

            if x4[i]>1./8.:
                d5[i]=2./8.
                z5[i]=4

            if x4[i]>3./8.:
                d5[i]=4./8.
                z5[i]=5
        
            if x4[i]>5./8.:
                d5[i]=6./8.
                z5[i]=6

            x5[i]=self.G*(x4[i]-d5[i])-.01
            
            z6[i]=0
            if x5[i]>-2/3:
                z6[i]=1
            
            if x5[i]>0:
                z6[i]=2
            
            if x5[i]>2/3:
                z6[i]=3

            OUT[i+4]=512*z1[i]+128*z2[i]+32*z3[i]+8*z4[i]+2*z5[i]+z6[i]
        
        OUT -= 2048
        n2 = np.arange(-N/2,N/2)
        f = n2*Fs/N
        win = kaiser(N,9.6)
        Yf=OUT*win
        Yf = 20*np.log10(.25*N)+20*np.log10(np.abs(fftshift(fft(Yf))))-20*np.log10(2048)
        Y = np.power(10,Yf/10.)
        self.f=f
        self.Yf=Yf
          
        


    def update(self):   
        self.calculate(self.toEngFormat(self.Fs),self.toEngFormat(self.Fin),self.toEngFormat(self.N))
        #self.figure.clf() 
        #self.axes1=self.figure.add_subplot(111)        
        self.axes1.plot(self.f,self.Yf,color='blue')
        self.figure.canvas.draw() 
        
    def toEngFormat(self,value):
        decimal.getcontext().prec = self.precision
        temp = decimal.Decimal(value)
        value_eng = temp.normalize().to_eng_string()
        return value_eng

    def go(self):
        self.configure_traits()
    
if __name__ == "__main__":
    main = pipelineADC()
    main.go()



