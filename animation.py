import os,sys
import numpy as np
import scipy
import time
#from PyQt5 import QtGui,QtCore,QtWidgets
#from PyQt5.QtCore import Qt, QSize
#from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
#from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.pyplot import Circle
#from matplotlib.patches import Arrow
#import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec
#import math
from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.mlab as mlab
#from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib as mpl
import imageio
from pwlib import *





class anim():
    
    def __init__(self,balwan, signals):
        self.balwan = balwan 
        self.v = signals[0] 
        self.k = signals[1]
        self.rv = signals[2]
        self.period = self.v.period
        self.play()

    def play(self):
        
        
        self.figure = plt.figure(figsize=(18,10))
        size = len(self.v.x)
        ind = np.arange(0,2*size,20)
        with imageio.get_writer('mygif.gif', mode='I',fps=5) as writer:
            
            for j in ind:
                
                               
                #print j
                self.figure.clf()
                
                try:
                    gs = gridspec.GridSpec(ncols=5, nrows=4,figure=self.figure)
                except:		
                    gs = gridspec.GridSpec(ncols=5, nrows=4,wspace=0.5,hspace=0.5)

                
                
                phasename = '$Phase$'
                self.f_v = self.figure.add_subplot(gs[0,0:2])
                self.f_v.set_ylabel('$V$ $[mag]$',fontsize=14)
                #self.f_v.set_xlabel(phasename,fontsize=14)
                self.f_v.set_xlim(-0.1,2.1)
                f = 0.05*np.abs(np.min(self.v.y)-np.max(self.v.y))
                self.f_v.set_ylim(np.min(self.v.y)-f,np.max(self.v.y)+f)
                if j < size:
                    i = j
                    self.f_v.plot(self.v.phase,self.v.mag,'.',color='blue')
                    self.f_v.plot(self.v.x[0:i],self.v.y[0:i],'-',color='black')
                    self.f_v.plot(self.v.x[i:],self.v.y[i:],'-',color='silver')
                    self.f_v.plot(self.v.x[i],self.v.y[i],'o',color='red')
                    self.f_v.plot(self.v.x+1,self.v.y,'-',color='silver')
                    self.f_v.plot(self.v.phase+1,self.v.mag,'.',color='blue')
                else:
                    i = j-size
                    self.f_v.plot(self.v.phase+1,self.v.mag,'.',color='blue')
                    self.f_v.plot(self.v.x[0:i]+1,self.v.y[0:i],'-',color='black')
                    self.f_v.plot(self.v.x[i:]+1,self.v.y[i:],'-',color='silver')
                    self.f_v.plot(self.v.x[i]+1,self.v.y[i],'o',color='red')
                    self.f_v.plot(self.v.phase,self.v.mag,'.',color='blue')
                    self.f_v.plot(self.v.x,self.v.y,'-',color='black')
                    
                self.f_k = self.figure.add_subplot(gs[1,0:2],sharex=self.f_v)
                self.f_k.set_ylabel('$K$ $[mag]$',fontsize=14)
                #self.f_k.set_xlabel(phasename,fontsize=14)
                self.f_k.set_xlim(-0.1,2.1)
                f = 0.05*np.abs(np.min(self.k.y)-np.max(self.k.y))
                self.f_k.set_ylim(np.min(self.k.y)-f,np.max(self.k.y)+f)
                if j < size:
                    i = j
                    self.f_k.plot(self.k.phase,self.k.mag,'.',color='blue')
                    self.f_k.plot(self.k.x[0:i],self.k.y[0:i],'-',color='black')
                    self.f_k.plot(self.k.x[i:],self.k.y[i:],'-',color='silver')
                    self.f_k.plot(self.k.x[i],self.k.y[i],'o',color='red')
                    self.f_k.plot(self.k.phase+1,self.k.mag,'.',color='blue')
                    self.f_k.plot(self.k.x+1,self.k.y,'-',color='silver')
                else:
                    i = j-size
                    self.f_k.plot(self.k.phase+1,self.k.mag,'.',color='blue')
                    self.f_k.plot(self.k.x[0:i]+1,self.k.y[0:i],'-',color='black')
                    self.f_k.plot(self.k.x[i:]+1,self.k.y[i:],'-',color='silver')
                    self.f_k.plot(self.k.x[i]+1,self.k.y[i],'o',color='red')
                    self.f_k.plot(self.k.phase,self.k.mag,'.',color='blue')
                    self.f_k.plot(self.k.x,self.k.y,'-',color='black')

        

                self.f_vk = self.figure.add_subplot(gs[2,0:2],sharex=self.f_v)
                self.f_vk.set_ylabel('$(V-K)$ $[mag]$',fontsize=14)
                #self.f_vk.set_xlabel(phasename,fontsize=14)
                self.f_vk.set_xlim(-0.1,2.1)
                f = 0.05*np.abs(np.min(self.v.y-self.k.y)-np.max(self.v.y-self.k.y))

                
                self.f_vk.set_ylim(np.min(self.v.y-self.k.y)-f,np.max(self.v.y-self.k.y)+f)
                if j < size:
                    i = j
                    #self.f_vk.plot(self.period*self.v.phase,self.v.mag-self.k.mag,'.',color='blue')
                    self.f_vk.plot(self.v.x[0:i],self.v.y[0:i]-self.k.y[0:i],'-',color='black')
                    self.f_vk.plot(self.v.x[i:],self.v.y[i:]-self.k.y[i:],'-',color='silver')
                    self.f_vk.plot(self.v.x[i],self.v.y[i]-self.k.y[i],'o',color='red')
                    self.f_vk.plot(self.v.x+1,self.v.y-self.k.y,'-',color='silver')
                
                else:
                    i = j-size
                    self.f_vk.plot(self.v.x[0:i]+1,self.v.y[0:i]-self.k.y[0:i],'-',color='black')
                    self.f_vk.plot(self.v.x[i:]+1,self.v.y[i:]-self.k.y[i:],'-',color='silver')
                    self.f_vk.plot(self.v.x[i]+1,self.v.y[i]-self.k.y[i],'o',color='red')
                    self.f_vk.plot(self.v.x,self.v.y-self.k.y,'-',color='black')

                color_min = np.min(self.v.y-self.k.y)
                color_max = np.max(self.v.y-self.k.y)
                color_curr = 0+int(255*(self.v.y[i]-self.k.y[i]-color_min)/(color_max-color_min))
                #print color_curr
                color = (1,1,1.-float(color_curr/255.))


                self.f_rv = self.figure.add_subplot(gs[3,0:2])
                self.f_rv.set_ylabel('$V_r-V_{\gamma}$ $[km/s]$',fontsize=14)
                self.f_rv.set_xlabel(phasename,fontsize=14)	
                self.f_rv.set_xlim(-0.1,2.1)
                f = 0.05*np.abs(np.min(self.rv.y)-np.max(self.rv.y))

                self.f_rv.set_ylim(np.min(self.rv.y)-f,np.max(self.rv.y)+f)
                if j < size:
                    i = j
                    self.f_rv.plot(self.rv.phase,self.rv.mag,'.',color='blue')
                    self.f_rv.plot(self.rv.x[0:i],self.rv.y[0:i],'-',color='black')
                    self.f_rv.plot(self.rv.x[i:],self.rv.y[i:],'-',color='silver')
                    self.f_rv.plot(self.rv.x[i],self.rv.y[i],'o',color='red')
                    self.f_rv.plot(self.rv.x+1,self.rv.y,'-',color='silver')
                    self.f_rv.plot(self.rv.phase+1,self.rv.mag,'.',color='blue')
                else:
                    i = j-size
                    self.f_rv.plot(self.rv.phase+1,self.rv.mag,'.',color='blue')
                    self.f_rv.plot(self.rv.x[0:i]+1,self.rv.y[0:i],'-',color='black')
                    self.f_rv.plot(self.rv.x[i:]+1,self.rv.y[i:],'-',color='silver')
                    self.f_rv.plot(self.rv.x[i]+1,self.rv.y[i],'o',color='red')
                    self.f_rv.plot(self.rv.x,self.rv.y,'-',color='black')
                    self.f_rv.plot(self.rv.phase,self.rv.mag,'.',color='blue')

                
                
                self.ax = self.figure.add_subplot(gs[:,2:],projection='3d')
                

                '''self.f_star = self.figure.add_subplot(gs[:,2:])
                self.f_star.set_title('$BW$',fontsize=14)
                self.f_star.set_ylabel('$\Theta$ $[mas]$',fontsize=14)
                self.f_star.set_xlabel('$R[R_{\odot}]$',fontsize=14)
                
                self.f_rv.tick_params(bottom=True,top=True,left=True,right=True,labelsize=12,direction='in')
                self.f_v.tick_params(bottom=True,top=True,left=True,right=True,labelsize=12,direction='in')
                self.f_k.tick_params(bottom=True,top=True,left=True,right=True,labelsize=12,direction='in')
                self.f_vk.tick_params(bottom=True,top=True,left=True,right=True,labelsize=12,direction='in')
                self.f_star.tick_params(bottom=True,top=True,left=True,right=True,labelsize=12,direction='in')'''

                u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:50j]
                dst = 3.26*9.46*np.power(10.,15.)/self.balwan.plx
                #print dst
                rsun = 696340.*2.
                x = dst*self.balwan.fi_plot[i]*np.cos(u)*np.sin(v)/rsun
                y = dst*self.balwan.fi_plot[i]*np.sin(u)*np.sin(v)/rsun
                z = dst*self.balwan.fi_plot[i]*np.cos(v)/rsun
                
                
                #ax.plot_wireframe(x, y, z, color="black",linewidth=0.5)
                self.ax.plot_surface(x, y, z, color=color)
                self.ax.set_xlim(-1*dst*np.max(self.balwan.fi_plot)/rsun,dst*np.max(self.balwan.fi_plot)/rsun)
                self.ax.set_ylim(-1*dst*np.max(self.balwan.fi_plot)/rsun,dst*np.max(self.balwan.fi_plot)/rsun)
                self.ax.set_zlim(-1.*dst*np.max(self.balwan.fi_plot)/rsun,dst*np.max(self.balwan.fi_plot)/rsun)
                self.ax.set_ylabel('$Y[R_{\odot}]$',fontsize=14,color="white")
                self.ax.set_xlabel('$X[R_{\odot}]$',fontsize=14,color="white")
                self.ax.set_zlabel('$Z[R_{\odot}]$',fontsize=14,color="white")
                
                self.ax.set_facecolor('xkcd:black')
                self.ax.xaxis.label.set_color('white')
                self.ax.tick_params(axis='x', colors='white')
                self.ax.tick_params(axis='y', colors='white')
                self.ax.tick_params(axis='z', colors='white')
                
                self.ax.text2D(0.1,0.95,self.balwan.name+'\tP='+str(self.v.period)+'d \t Phase='+str("{:.1f}".format(self.v.x[i]))+'\t$R=$'+str("{:.2f}".format(dst*self.balwan.fi_plot[i]/rsun))+'$R_{\odot}$',color='white',fontsize=15, transform=self.ax.transAxes)
        
                self.f_v.invert_yaxis()
                self.f_k.invert_yaxis()
                self.f_vk.invert_yaxis()
                plt.savefig(str(i)+'.png')
                image = imageio.imread(str(i)+'.png')
                os.system('rm '+str(i)+'.png')
                writer.append_data(image)