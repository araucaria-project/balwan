#!/usr/bin/env python3
import os,sys
import numpy as np
import scipy
#from astropy.io import fits as pyfits
#from astropy import wcs
from PyQt5 import QtGui,QtCore,QtWidgets
from PyQt5.QtCore import Qt, QSize
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.pyplot import Circle
from matplotlib.patches import Arrow
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.mlab as mlab
from scipy.stats import norm
import matplotlib.pyplot as plt
from akimafit import *
from bw_lib import *
from bw_lib_pw import *
from pwlib import *
from animation import *
#1- wczytywanie V,K,rv,p
#2- fitowanie akima z mc z niepewnosciami
#3- usuwanie ewentualnych outlierow i refitowanie
#4- wybor przesuniecia fazowego, ekstynkcji i prawa poczerwienienia, zaleznosci irsb
#5- calkowanie rv z mc, wyznaczanie krzywej delta_fi z mc, fitowanie BW z mc
#6- usuwanie outlierow z BW, refitowanie BW

#irsb file is specified below, its format: name term0 term1 term2 ..., one relation in line, without empty lines, lines can be commented with #
#TODO zrobic przyciski do alternatywnego wybierania algorytmu fitowania metnk/bisect... dopracowac bisect, bo moja uproszczona wersja daje stabilne wyniki (zawsze bliskie metnk)

irsbfile = '/home/piotr/src/balwan/irsb.dat'	

def parse_file(plik):
	dane=os.popen('cat '+plik).read().split('\n')[:-1]
	pv = dane[0]
	pk = dane[1]
	prv =dane[2]
	period = dane[3].split('=')[1].split(',')
	
	for i,p in enumerate(period):
		period[i] = float(period[i])
	ebv = float(dane[4].split('=')[1])
	ebverr = float(dane[5].split('=')[1])
	plx = float(dane[6].split('=')[1])
	plxerr = float(dane[7].split('=')[1])
	psk= float(dane[8].split('=')[1])
	psrv = float(dane[9].split('=')[1])
	name = dane[10].split('=')[1]

	plikv = pv.split('[')[0]
	tv = int(pv.split('[')[1].split(']')[0].split(',')[0])
	mv =  int(pv.split('[')[1].split(']')[0].split(',')[1])
	ev =  int(pv.split('[')[1].split(']')[0].split(',')[2])	
	kv =  float(pv.split()[1])
	nv = float(pv.split()[2])

	plikk = pk.split('[')[0]
	tk = int(pk.split('[')[1].split(']')[0].split(',')[0])
	mk =  int(pk.split('[')[1].split(']')[0].split(',')[1])
	ek =  int(pk.split('[')[1].split(']')[0].split(',')[2])	
	kk =  float(pk.split()[1])
	nk = float(pk.split()[2])

	plikrv = prv.split('[')[0]
	trv = int(prv.split('[')[1].split(']')[0].split(',')[0])
	mrv =  int(prv.split('[')[1].split(']')[0].split(',')[1])
	erv =  int(prv.split('[')[1].split(']')[0].split(',')[2])	
	krv =  float(prv.split()[1])
	nrv = float(prv.split()[2])

	return [plikv,tv,mv,ev,kv,nv],[plikk,tk,mk,ek,kk,nk],[plikrv,trv,mrv,erv,krv,nrv],period,ebv,ebverr,plx,plxerr,psk,psrv,name

def readfile(plik):
	d = os.popen('cat ' + plik).read()
	if '\n' in d:
		dane = d.split('\n')[:-1]
	elif '\r' in d:
		dane = d.split('\n')[:-1]
	else:
		print('Cannot parse an input file\n')
		sys.exit()
	data = []
	nrow = len(dane)
	ncol = 0
	k=0	
	for i,linia in enumerate(dane):
		#print linia
		if '#' not in linia:
			if ',' in linia:
				line_data = linia.split(',')
			else:
				line_data = linia.split()
			for j,item in enumerate(line_data):
				if k == 0:
					ncol = ncol+1
					data.append([item])

				else:
					data[j].append(item)

			k=1

	return data,ncol,nrow
	


class syst_test():
	def __init__(self):
		self.param = ''
		self.unit = ''
		self.values = []
		self.p = []
		self.r0 = []
		self.sigma = []
			

class display(QtWidgets.QWidget):


	def __init__(self, parent=None):
		super(display,self).__init__(parent)

		self.initUI()

	def initUI(self):
		self.redd_v = 3.136
		self.redd_k = 0.363
		self.oldp = 0
		self.irsb = []
		self.suggested_colors = ['blue','green','red','cyan','magenta','yellow','black','pink']
		
		#self.check_buttons = []
		self.name = None
		self.signals = []
		self.bw = None	
			
		self.whichpt = None
		self.whichpt_isdel = False
		self.whichpt_isrej = False
		self.plot_model=False

		self.inittcol = 1
		self.initmagcol = 2
		self.initerrcol = 3
		self.initcolor = self.suggested_colors[0]
		self.initmodelcolor = self.suggested_colors[0]
		#----input file------
		#---button to load input file-----
		#self.select_starsList=QtWidgets.QFileDialog()
		
		self.b1 = QtWidgets.QPushButton(self)#image button
		self.b1.setText("Load data")
		self.b1.clicked.connect(self.b1_clicked)

		self.b1_2 = QtWidgets.QPushButton(self)#image button
		self.b1_2.setText("Load conffile")
		self.b1_2.clicked.connect(self.b1_2_clicked)

		self.b1_3 = QtWidgets.QPushButton(self)#image button
		self.b1_3.setText("Export data")
		self.b1_3.clicked.connect(self.b1_3_clicked)

		self.b2 = QtWidgets.QPushButton(self)#1st file button
		self.b2.setText("Model curves")
		self.b2.clicked.connect(self.b2_clicked)

		self.b2_2 = QtWidgets.QPushButton(self)#1st file button
		self.b2_2.setText("BaWAn")
		self.b2_2.clicked.connect(self.b2_2_clicked)

		self.b2_3 = QtWidgets.QPushButton(self)#1st file button
		self.b2_3.setText("Exit")
		self.b2_3.clicked.connect(self.b2_3_clicked)

		#----edit parameters-------
	
		self.period_label = QtWidgets.QLabel('Period:',self)
		self.period_field = QtWidgets.QLineEdit()
		self.period_field.setText(str(0))

		self.period_change_label = QtWidgets.QLabel('Period change (d/yr):',self)
		self.period_change_field = QtWidgets.QLineEdit()
		self.period_change_field.setText(str(0))

		self.period_change_t0_label = QtWidgets.QLabel('Period change HJD0:',self)
		self.period_change_t0_field = QtWidgets.QLineEdit()
		self.period_change_t0_field.setText(str(0))
		
		self.shift_label= QtWidgets.QLabel('Phase shift V-K:',self)
		self.shift_field = QtWidgets.QLineEdit()
		self.shift_field.setText(str(0))

		self.shift_rv_label= QtWidgets.QLabel('Phase shift RV:',self)
		self.shift_rv_field = QtWidgets.QLineEdit()
		self.shift_rv_field.setText(str(0))
	
		self.minph_label = QtWidgets.QLabel('Minimum phase:',self)
		self.minph_field = QtWidgets.QLineEdit()
		self.minph_field.setText(str(0.0))

		self.maxph_label = QtWidgets.QLabel('Maximum phase:',self)
		self.maxph_field = QtWidgets.QLineEdit()
		self.maxph_field.setText(str(1.0))

		self.ebv_label = QtWidgets.QLabel('E(B-V):',self)
		self.ebv_field = QtWidgets.QLineEdit()
		self.ebv_field.setText(str(0.))
		self.ebverr_label = QtWidgets.QLabel('+-',self)
		self.ebverr_field = QtWidgets.QLineEdit()
		self.ebverr_field.setText(str(0.))

		self.irsb_label = QtWidgets.QLabel('IRSB:',self)
		'''self.irsba_field = QtWidgets.QLineEdit()
		self.irsba_field.setText(str(-0.1336))
		self.irsb_label2 = QtWidgets.QLabel('*(V-K)+',self)
		self.irsbb_field = QtWidgets.QLineEdit()
		self.irsbb_field.setText(str(3.9530))'''
		self.irsb_list = QtWidgets.QComboBox()
		global irsbfile
		self.addirsb(irsbfile)

		self.v_k_label = QtWidgets.QLabel('V model knots:',self)
		self.v_k_field = QtWidgets.QLineEdit()
		self.v_k_field.setText(str(0.1))
		self.v_n_label = QtWidgets.QLabel('bins size:',self)
		self.v_n_field = QtWidgets.QLineEdit()
		self.v_n_field.setText(str(0.1))

		self.k_k_label = QtWidgets.QLabel('K/V-K/phi model knots:',self)
		self.k_k_field = QtWidgets.QLineEdit()
		self.k_k_field.setText(str(0.1))
		self.k_n_label = QtWidgets.QLabel('bins size:',self)
		self.k_n_field = QtWidgets.QLineEdit()
		self.k_n_field.setText(str(0.1))
		
		self.rv_k_label = QtWidgets.QLabel('RV model knots:',self)
		self.rv_k_field = QtWidgets.QLineEdit()
		self.rv_k_field.setText(str(0.1))
		self.rv_n_label = QtWidgets.QLabel('bins size:',self)
		self.rv_n_field = QtWidgets.QLineEdit()
		self.rv_n_field.setText(str(0.1))

		self.buttons_group1 = QtWidgets.QButtonGroup()
		self.b_plx = QtWidgets.QRadioButton("Parallax [mas]")
		self.b_plx.setChecked(True)
		
		self.plx_field = QtWidgets.QLineEdit()
		self.plx_field.setText(str(0.))
		self.plx_label = QtWidgets.QLabel('+-',self)
		self.plxerr_field = QtWidgets.QLineEdit()
		self.plxerr_field.setText(str(0.))
		
		self.b_pf = QtWidgets.QRadioButton("p-factor")
		self.pf_field = QtWidgets.QLineEdit()
		self.pf_field.setText(str(0.))
		self.pf_label = QtWidgets.QLabel('+-',self)
		self.pferr_field = QtWidgets.QLineEdit()
		self.pferr_field.setText(str(0.))
		#self.b_dst.setChecked(False)
		self.buttons_group1.addButton(self.b_plx)
		self.buttons_group1.addButton(self.b_pf)


		self.buttons_group_fit_method = QtWidgets.QButtonGroup()
		self.method_label = QtWidgets.QLabel('Fit method',self)
		self.b_drs = QtWidgets.QRadioButton("DR scaling")
		self.b_drs.setChecked(True)
		self.b_lsq = QtWidgets.QRadioButton("Line (OLS)")
		self.b_lsq.setChecked(False)
		self.b_bis = QtWidgets.QRadioButton("Line (bisector)")
		self.b_bis.setChecked(False)
	
		#
		self.buttons_group_fit_method.addButton(self.b_drs)
		self.buttons_group_fit_method.addButton(self.b_lsq)
		self.buttons_group_fit_method.addButton(self.b_bis)
		
		#text field with results and input info
		self.fit_results = QtWidgets.QLabel('-----Load data-----',self)
		self.fit_results.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
		self.fit_results.setStyleSheet("background-color: lightgreen; border: 1px solid black;")

		self.scroll = QtWidgets.QScrollArea()
		self.scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
		self.scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
		self.scroll.setWidgetResizable(True)
		self.scroll.setWidget(self.fit_results)

		#text field with point data
		self.pt_info = QtWidgets.QLabel('-----Select data point-----\n-----with doubleclick-----',self)
		self.pt_info.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
		self.pt_info.setStyleSheet("background-color: white; border: 1px solid black;")
		self.scroll2 = QtWidgets.QScrollArea()
		self.scroll2.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
		self.scroll2.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
		self.scroll2.setWidgetResizable(True)
		self.scroll2.setWidget(self.pt_info)
		
		self.mc_label = QtWidgets.QLabel('Monte Carlo simulation:',self)
		self.mc_field = QtWidgets.QLineEdit()
		self.mc_field.setText(str(1000))
		self.mc_button = QtWidgets.QCheckBox()
		self.mc_button.setChecked(False)

		self.b_testsys = QtWidgets.QPushButton()
		self.b_testsys.setText('Systematic test')
		self.b_testsys.clicked.connect(self.showdialog_testsyst)

		self.b_mcopt = QtWidgets.QPushButton()
		self.b_mcopt.setText('MC optimisation')
		self.b_mcopt.clicked.connect(self.showdialog_mcopt)
				
		self.b_anim = QtWidgets.QPushButton()
		self.b_anim.setText('Animation')
		self.b_anim.clicked.connect(self.b_anim_clicked)

		#------------------display image-----------------------------

		# image
		self.figure = Figure(figsize=(15,4),constrained_layout=True)
		#self.figure.tight_layout()
		self.canvas = FigureCanvas(self.figure)
		self.toolbar = NavigationToolbar(self.canvas,self)
		self.update_plot(0)
		self.scroll_plot = QtWidgets.QScrollArea()
		self.scroll_plot.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
		self.scroll_plot.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
		self.scroll_plot.setWidgetResizable(True)
		self.scroll_plot.setWidget(self.canvas)

		
		#---------------LOGO----------------------------------------------
		
		self.logo = QtWidgets.QLabel('',self)
		pxmap = QtGui.QPixmap('balwan.png')
		pxmap2 = pxmap.scaledToWidth(56)
		self.logo.setPixmap(pxmap2)

		
		self.logo2 = QtWidgets.QLabel('',self)
		pxmap = QtGui.QPixmap('ara_logo.png')
		pxmap3 = pxmap.scaledToWidth(128)
		self.logo2.setPixmap(pxmap3)

		
		
		#--------------GUI SETTINGS----------------------------------------
		

		hbox_main=QtWidgets.QHBoxLayout()#main
		vbox_left=QtWidgets.QVBoxLayout()#left
		vbox_center=QtWidgets.QVBoxLayout()#left
		vbox_right=QtWidgets.QVBoxLayout()#right
	
		hbox_logos=QtWidgets.QHBoxLayout()
		hbox_logos.addWidget(self.logo,1)
		hbox_logos.addWidget(self.logo2,2)

		vbox_left.addLayout(hbox_logos)
		
		h1 = QtWidgets.QHBoxLayout()
		h1_1 = QtWidgets.QHBoxLayout()
		h1_2 = QtWidgets.QHBoxLayout()
		h2 = QtWidgets.QHBoxLayout()
		h2_1 = QtWidgets.QHBoxLayout()
		h3 = QtWidgets.QHBoxLayout()
		h4 = QtWidgets.QHBoxLayout()
		h5 = QtWidgets.QHBoxLayout()
		h6 = QtWidgets.QHBoxLayout()
		h7 = QtWidgets.QHBoxLayout()

		h8 = QtWidgets.QHBoxLayout()
		h9 = QtWidgets.QHBoxLayout()
		h10 = QtWidgets.QHBoxLayout()
		h11 = QtWidgets.QHBoxLayout()

		h12 = QtWidgets.QHBoxLayout()
		h13 = QtWidgets.QHBoxLayout()
		
		hbox_buttons=QtWidgets.QHBoxLayout()#buttons
		hbox_buttons2=QtWidgets.QHBoxLayout()#buttons
		hbox_buttons3=QtWidgets.QHBoxLayout()
		hbox_buttons.addWidget(self.b1)
		hbox_buttons.addWidget(self.b1_2)
		hbox_buttons.addWidget(self.b1_3)
		vbox_left.addLayout(hbox_buttons,1)

		
		h1.addWidget(self.period_label)
		h1.addWidget(self.period_field)
		h1_1.addWidget(self.period_change_label)
		h1_1.addWidget(self.period_change_field)
		h1_2.addWidget(self.period_change_t0_label)
		h1_2.addWidget(self.period_change_t0_field)
		h2.addWidget(self.shift_label)
		h2.addWidget(self.shift_field)
		h2_1.addWidget(self.shift_rv_label)
		h2_1.addWidget(self.shift_rv_field)
		h3.addWidget(self.minph_label)
		h3.addWidget(self.minph_field)
		h4.addWidget(self.maxph_label)
		h4.addWidget(self.maxph_field)
		h5.addWidget(self.ebv_label)
		h5.addWidget(self.ebv_field)
		h5.addWidget(self.ebverr_label)
		h5.addWidget(self.ebverr_field)
		h6.addWidget(self.irsb_label)
		h6.addWidget(self.irsb_list)
		#h6.addWidget(self.irsba_field)
		#h6.addWidget(self.irsb_label2)
		#h6.addWidget(self.irsbb_field)
		h7.addWidget(self.mc_label)
		h7.addWidget(self.mc_button)
		h7.addWidget(self.mc_field)

		h8.addWidget(self.v_k_label)
		h8.addWidget(self.v_k_field)
		h8.addWidget(self.v_n_label)
		h8.addWidget(self.v_n_field)
		
		h9.addWidget(self.k_k_label)
		h9.addWidget(self.k_k_field)
		h9.addWidget(self.k_n_label)
		h9.addWidget(self.k_n_field)

		h10.addWidget(self.rv_k_label)
		h10.addWidget(self.rv_k_field)
		h10.addWidget(self.rv_n_label)
		h10.addWidget(self.rv_n_field)

		h11.addWidget(self.method_label)
		h11.addWidget(self.b_drs)
		h11.addWidget(self.b_lsq)
		h11.addWidget(self.b_bis)

		h12.addWidget(self.b_plx)
		h12.addWidget(self.plx_field)
		h12.addWidget(self.plx_label)
		h12.addWidget(self.plxerr_field)

		h13.addWidget(self.b_pf)
		h13.addWidget(self.pf_field)
		h13.addWidget(self.pf_label)
		h13.addWidget(self.pferr_field)


		hbox_buttons2.addWidget(self.b2)
		hbox_buttons2.addWidget(self.b2_2)
		hbox_buttons2.addWidget(self.b2_3)
		hbox_buttons3.addWidget(self.b_testsys)
		hbox_buttons3.addWidget(self.b_mcopt)
		hbox_buttons3.addWidget(self.b_anim)
		vbox_left.addLayout(h1,1)
		vbox_left.addLayout(h1_1,1)
		vbox_left.addLayout(h1_2,1)
		vbox_left.addLayout(h2,1)
		vbox_left.addLayout(h2_1,1)
		vbox_left.addLayout(h3,1)
		vbox_left.addLayout(h4,1)
		vbox_left.addLayout(h5,1)
		vbox_left.addLayout(h6,1)
		vbox_left.addLayout(h8,1)
		vbox_left.addLayout(h9,1)
		vbox_left.addLayout(h10,1)
		vbox_left.addLayout(h12,1)
		vbox_left.addLayout(h13,1)
		vbox_left.addLayout(h7,1)
		vbox_left.addLayout(h11,1)
		vbox_left.addLayout(hbox_buttons2,1)
		vbox_left.addLayout(hbox_buttons3,1)
		vbox_left.addWidget(self.scroll,3)
		vbox_left.addWidget(self.scroll2,3)

		

		vbox_right.addWidget(self.scroll_plot,1)
		vbox_right.addWidget(self.toolbar)
		
		
		
		#----gui settings-------
		
		hbox_main.addLayout(vbox_left,1)
		#hbox_main.addLayout(vbox_center,3)
		hbox_main.addLayout(vbox_right,5)
		self.setLayout(hbox_main)
		self.setGeometry(500,300,1400,1000)
		self.setWindowTitle("Baade+Wesselink Analysis (Ba+WAn)")
		self.setWindowIcon(QtGui.QIcon('balwan.png'))
		self.show()

	def colorchanged1(self,text):
		
		try:
			self.color_field.setStyleSheet("color: "+self.color_field.text())
		except:
			self.color_field.setStyleSheet("color: black")

	def addirsb(self,irsbfl):
		text = os.popen('cat '+irsbfl).read().split('\n')[:-1]
		
		for linia in text:
			terms = []
			if '#' not in linia:
				dane = linia.split()
				name = dane[0]
				for el in dane[1:]:
					terms.append(float(el))

			self.irsb.append(terms)
			self.irsb_list.addItem(name)
			self.wirsb()

	def wirsb(self):
		return self.irsb_list.currentIndex()
		

	def b1_clicked(self):
		
		self.showdialog_newband()

	def b_anim_clicked(self):
		
		anim(self.bw,self.signals)

	def b1_2_clicked(self):
		
		select_file = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', '.')
		if len(select_file[0]) > 1:

			self.signals = []

			f1,f2,f3,period,ebv,ebverr,plx,plxerr,psk,psrv,name = parse_file(select_file[0])

			self.signals.append(signal(f1[0],f1[1],f1[2],f1[3],period=period,k=f1[4],n=f1[5],name='',sc=None,hjd0=None,flux=True,ph=True))
			self.signals.append(signal(f2[0],f2[1],f2[2],f2[3],period=period,k=f2[4],n=f2[5],name='',sc=None,hjd0=self.signals[0].hjd0,flux=True,ph=True))
			self.signals.append(signal(f3[0],f3[1],f3[2],f3[3],period=period,k=f3[4],n=f3[5],name='',sc=None,hjd0=self.signals[0].hjd0,flux=False,ph=True,tomean=True))
			self.v_k_field.setText(str(f1[4]))
			self.v_n_field.setText(str(f1[5]))
			self.k_k_field.setText(str(f2[4]))
			self.k_n_field.setText(str(f2[5]))
			self.rv_k_field.setText(str(f3[4]))
			self.rv_n_field.setText(str(f3[5]))

			self.ebv_field.setText(str(ebv))
			self.ebverr_field.setText(str(ebverr))
			self.plx_field.setText(str(plx))
			self.plxerr_field.setText(str(plxerr))
			self.shift_field.setText(str(psk))
			self.shift_rv_field.setText(str(psrv))

			self.period_field.setText(str(period[0]))
			try:
				self.period_change_field.setText(str(period[1]))
			except:
				self.period_change_field.setText(str(0))

			try:
				self.period_change_t0_field.setText(str(period[2]))
			except:
				self.period_change_t0_field.setText(str(0))


			self.name = name

		self.b2_clicked()
		self.b2_2_clicked()
		self.update_plot(0)
		

	def b1_3_clicked(self):
		
		self.showdialog_export()

	def b2_clicked(self):
		#print self.signals[0].period
		if self.signals[0].period[0] != float(self.period_field.text()):
			self.signals[0].period[0] = float(self.period_field.text())
			self.signals[0].phasing()
			self.signals[1].period[0] = float(self.period_field.text())
			self.signals[1].phasing()
			self.signals[2].period[0] = float(self.period_field.text())
			self.signals[2].phasing()

		if self.signals[0].period[1] != float(self.period_change_field.text()):
			self.signals[0].period[1] = float(self.period_change_field.text())
			self.signals[0].phasing()
			self.signals[1].period[1] = float(self.period_change_field.text())
			self.signals[1].phasing()
			self.signals[2].period[1] = float(self.period_change_field.text())
			self.signals[2].phasing()

		if self.signals[0].period[2] != float(self.period_change_t0_field.text()):
			self.signals[0].period[2] = float(self.period_change_t0_field.text())
			self.signals[0].phasing()
			self.signals[1].period[2] = float(self.period_change_t0_field.text())
			self.signals[1].phasing()
			self.signals[2].period[2] = float(self.period_change_t0_field.text())
			self.signals[2].phasing()

		
			
		if float(self.v_k_field.text()) != self.signals[0].k:
			self.signals[0].k = float(self.v_k_field.text())

			
		if float(self.v_n_field.text()) != self.signals[0].n:
			self.signals[0].n = float(self.v_n_field.text())

		if float(self.k_k_field.text()) != self.signals[1].k:
			self.signals[1].k = float(self.k_k_field.text())

			
		if float(self.k_n_field.text()) != self.signals[1].n:
			self.signals[1].n = float(self.k_n_field.text())

			
		if float(self.rv_k_field.text()) != self.signals[2].k:
			self.signals[2].k = float(self.rv_k_field.text())

		if float(self.rv_n_field.text()) != self.signals[2].n:
			self.signals[2].n = float(self.rv_n_field.text())
		
		
		self.signals[0].fit()
		self.signals[1].ps=float(self.shift_field.text())
		self.signals[1].phasing()
		self.signals[1].fit()
		self.signals[2].ps=float(self.shift_rv_field.text())
		self.signals[2].phasing()
		self.signals[2].fit()
		if self.mc_button.isChecked():
			self.signals[0].mc(int(self.mc_field.text()),self.signals[1].phase)
			self.signals[1].mc(int(self.mc_field.text()),self.signals[1].phase)
			self.signals[2].mc(int(self.mc_field.text()),self.signals[1].phase)
			
		
		self.update_plot(0)


	def b2_2_clicked(self):
		if self.mc_button.isChecked():
			mc = int(self.mc_field.text())
		else:
			mc = 0	

		if self.b_drs.isChecked():
			method = 0

		if self.b_lsq.isChecked():
			method = 1

		if self.b_bis.isChecked():
			method = 2

		if self.b_plx.isChecked():

			self.bw=balwan(self.name,self.signals[0],self.signals[1],self.signals[2],float(self.ebv_field.text()),ebverr=float(self.ebverr_field.text()),R_V=self.redd_v,R_K=self.redd_k,plx=float(self.plx_field.text()),plxerr=float(self.plxerr_field.text()),p=None,perr=None,irsb=self.irsb[self.wirsb()],mc=mc,p0=float(self.minph_field.text()),p1=float(self.maxph_field.text()),method=method)

		else:
			self.bw=balwan(self.name,self.signals[0],self.signals[1],self.signals[2],float(self.ebv_field.text()),ebverr=float(self.ebverr_field.text()),R_V=self.redd_v,R_K=self.redd_k,plx=None,plxerr=None,p=float(self.pf_field.text()),perr=float(self.pferr_field.text()),irsb=self.irsb[self.wirsb()],mc=mc,p0=float(self.minph_field.text()),p1=float(self.maxph_field.text()),method=method)

		self.fit_results.setText(self.bw.text())

		self.oldp = self.bw.p
		self.update_plot(-1)

	def b2_3_clicked(self):
		self.close()


	def selectionchange(self,i):
		
		self.minx_field.setText(str(self.signals[i].lowxlimit))
		self.maxx_field.setText(str(self.signals[i].upxlimit))
		self.miny_field.setText(str(self.signals[i].lowylimit))
		self.maxy_field.setText(str(self.signals[i].upylimit))
		self.clip_field.setText(str(self.signals[i].nclipping))
		self.xshift_field.setText(str(self.signals[i].xshift))
		self.yshift_field.setText(str(self.signals[i].yshift))
		self.color_field.setText(str(self.signals[i].color))
		self.color_field.setStyleSheet("color: "+self.signals[i].color)
		self.fit_results.setText(self.signals[i].text)
		self.update_plot(i)

	

	def showdialog_newband(self):
		self.d = QtWidgets.QDialog()
		self.d.setWindowTitle("New signals")
		self.d.buttons_group1 = QtWidgets.QButtonGroup()
		self.d.buttons_group2 = QtWidgets.QButtonGroup()
		#self.d.setWindowModality(QtWidgets.ApplicationModal)
		
		
		self.d.name_label= QtWidgets.QLabel('Star name:',self)
		self.d.name_field = QtWidgets.QLineEdit()

		self.d.period_label= QtWidgets.QLabel('Pulsation period:',self)
		self.d.period_field = QtWidgets.QLineEdit()

		#------RV---------

		self.d.filerv_label= QtWidgets.QLabel('File RV:',self)
		self.d.filerv_field = QtWidgets.QLineEdit()
		#self.d.file_field.textChanged[str].connect(self.filechanged)
		self.d.b_select_filerv = QtWidgets.QPushButton(self)
		self.d.b_select_filerv.setText("Select file")
		self.d.b_select_filerv.clicked.connect(self.b_select_file_clicked_rv)

		self.d.rvt_label= QtWidgets.QLabel('time column:',self)
		self.d.rvt_field = QtWidgets.QLineEdit()
		self.d.rvt_field.setText(str(2))

		self.d.rv_label= QtWidgets.QLabel('RV column:',self)
		self.d.rv_field = QtWidgets.QLineEdit()
		self.d.rv_field.setText(str(8))

		self.d.rverr_label= QtWidgets.QLabel('error column:',self)
		self.d.rverr_field = QtWidgets.QLineEdit()
		self.d.rverr_field.setText(str(9))

		#------V---------
		
		self.d.filev_label= QtWidgets.QLabel('File V:',self)
		self.d.filev_field = QtWidgets.QLineEdit()
		#self.d.file_field.textChanged[str].connect(self.filechanged)
		self.d.b_select_filev = QtWidgets.QPushButton(self)
		self.d.b_select_filev.setText("Select file")
		self.d.b_select_filev.clicked.connect(self.b_select_file_clicked_v)

		self.d.vt_label= QtWidgets.QLabel('time column:',self)
		self.d.vt_field = QtWidgets.QLineEdit()
		self.d.vt_field.setText(str(2))

		self.d.v_label= QtWidgets.QLabel('V column:',self)
		self.d.v_field = QtWidgets.QLineEdit()
		self.d.v_field.setText(str(8))

		self.d.verr_label= QtWidgets.QLabel('error column:',self)
		self.d.verr_field = QtWidgets.QLineEdit()
		self.d.verr_field.setText(str(9))

		#------V---------
		
		self.d.filek_label= QtWidgets.QLabel('File K:',self)
		self.d.filek_field = QtWidgets.QLineEdit()
		#self.d.file_field.textChanged[str].connect(self.filechanged)
		self.d.b_select_filek = QtWidgets.QPushButton(self)
		self.d.b_select_filek.setText("Select file")
		self.d.b_select_filek.clicked.connect(self.b_select_file_clicked_k)

		self.d.kt_label= QtWidgets.QLabel('time column:',self)
		self.d.kt_field = QtWidgets.QLineEdit()
		self.d.kt_field.setText(str(2))

		self.d.k_label= QtWidgets.QLabel('K column:',self)
		self.d.k_field = QtWidgets.QLineEdit()
		self.d.k_field.setText(str(8))

		self.d.kerr_label= QtWidgets.QLabel('error column:',self)
		self.d.kerr_field = QtWidgets.QLineEdit()
		self.d.kerr_field.setText(str(9))

		self.d.b_save = QtWidgets.QPushButton(self)
		self.d.b_save.setText("Upload")
		self.d.b_save.clicked.connect(self.b_save_clicked)

		self.d.b_cancel = QtWidgets.QPushButton(self)
		self.d.b_cancel.setText("Cancel")
		self.d.b_cancel.clicked.connect(self.b_cancel_clicked)

		self.d.tabrv =  QtWidgets.QTableWidget()
		self.d.tabv =  QtWidgets.QTableWidget()
		self.d.tabk =  QtWidgets.QTableWidget()
		
		layout = QtWidgets.QHBoxLayout()
		layout_left = QtWidgets.QVBoxLayout()
		h1 = QtWidgets.QHBoxLayout()
		h1_2 = QtWidgets.QHBoxLayout()
		

		h2 = QtWidgets.QHBoxLayout()
		h3 = QtWidgets.QHBoxLayout()

		h4 = QtWidgets.QHBoxLayout()
		h5 = QtWidgets.QHBoxLayout()

		h6 = QtWidgets.QHBoxLayout()
		h7 = QtWidgets.QHBoxLayout()

		layout_right = QtWidgets.QVBoxLayout()
		hbox_radio_buttons1 = QtWidgets.QHBoxLayout()
		hbox_radio_buttons2 = QtWidgets.QHBoxLayout()
		hbuttons = QtWidgets.QHBoxLayout()

		h1.addWidget(self.d.name_label)
		h1.addWidget(self.d.name_field)

		h1_2.addWidget(self.d.period_label)
		h1_2.addWidget(self.d.period_field)

		h2.addWidget(self.d.filev_label)
		h2.addWidget(self.d.filev_field)
		h2.addWidget(self.d.b_select_filev)
		h3.addWidget(self.d.vt_label)
		h3.addWidget(self.d.vt_field)
		h3.addWidget(self.d.v_label)
		h3.addWidget(self.d.v_field)
		h3.addWidget(self.d.verr_label)
		h3.addWidget(self.d.verr_field)

		h4.addWidget(self.d.filek_label)
		h4.addWidget(self.d.filek_field)
		h4.addWidget(self.d.b_select_filek)
		h5.addWidget(self.d.kt_label)
		h5.addWidget(self.d.kt_field)
		h5.addWidget(self.d.k_label)
		h5.addWidget(self.d.k_field)
		h5.addWidget(self.d.kerr_label)
		h5.addWidget(self.d.kerr_field)

		h6.addWidget(self.d.filerv_label)
		h6.addWidget(self.d.filerv_field)
		h6.addWidget(self.d.b_select_filerv)
		h7.addWidget(self.d.rvt_label)
		h7.addWidget(self.d.rvt_field)
		h7.addWidget(self.d.rv_label)
		h7.addWidget(self.d.rv_field)
		h7.addWidget(self.d.rverr_label)
		h7.addWidget(self.d.rverr_field)
		
		hbuttons.addWidget(self.d.b_save)
		hbuttons.addWidget(self.d.b_cancel)
		
		layout_left.addLayout(h1)
		layout_left.addLayout(h1_2)
		layout_left.addLayout(h2)
		layout_left.addLayout(h3)
		layout_left.addLayout(h4)
		layout_left.addLayout(h5)
		layout_left.addLayout(h6)
		layout_left.addLayout(h7)
		layout_left.addLayout(hbuttons)

		layout_right.addWidget(self.d.tabv)
		layout_right.addWidget(self.d.tabk)
		layout_right.addWidget(self.d.tabrv)
		layout.addLayout(layout_left)
		layout.addLayout(layout_right)
		self.d.setLayout(layout)
		#self.setGeometry(500,300,700,500)
		self.d.setGeometry(500,300,1400,500)
		self.d.setWindowTitle("new population")
		self.d.exec_()

	def b_save_clicked(self):
		
		self.signals = []

		self.name=self.d.name_field.text()
		
		self.signals.append(signal(self.d.filev_field.text(),int(self.d.vt_field.text())-1,int(self.d.v_field.text())-1,int(self.d.verr_field.text())-1,period=float(self.d.period_field.text()),k=0.1,n=0.1,name='',sc=None,hjd0=None,flux=True,ph=True,downlimit=0,uplimit=100))
		self.signals.append(signal(self.d.filek_field.text(),int(self.d.kt_field.text())-1,int(self.d.k_field.text())-1,int(self.d.kerr_field.text())-1,period=float(self.d.period_field.text()),k=0.1,n=0.1,name='',sc=None,hjd0=self.signals[0].hjd0,flux=True,ph=True,downlimit=0,uplimit=100))
		self.signals.append(signal(self.d.filerv_field.text(),int(self.d.rvt_field.text())-1,int(self.d.rv_field.text())-1,int(self.d.rverr_field.text())-1,period=float(self.d.period_field.text()),k=0.1,n=0.1,name='',sc=None,hjd0=self.signals[0].hjd0,flux=False,ph=True,tomean=True))	

		self.update_plot(-1)

		#self.fit_results.setText(str(self.signals[-1].text))

	def b_select_file_clicked_rv(self):
		select_file = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file',  '.')
		self.d.filerv_field.setText(select_file[0])

		data,ncol,nrow = readfile(select_file[0])
		
		self.d.tabrv.setRowCount(nrow)
		self.d.tabrv.setColumnCount(ncol)
		for i,col in enumerate(data):
			for j,record in enumerate(col):
				newitem = QtWidgets.QTableWidgetItem(record)
				self.d.tabrv.setItem(j, i, newitem) 

	def b_select_file_clicked_v(self):
		select_file = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file',  '.')
		self.d.filev_field.setText(select_file[0])

		data,ncol,nrow = readfile(select_file[0])
		
		self.d.tabv.setRowCount(nrow)
		self.d.tabv.setColumnCount(ncol)
		for i,col in enumerate(data):
			for j,record in enumerate(col):
				newitem = QtWidgets.QTableWidgetItem(record)
				self.d.tabv.setItem(j, i, newitem) 

	def b_select_file_clicked_k(self):
		select_file = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', '.')
		self.d.filek_field.setText(select_file[0])

		data,ncol,nrow = readfile(select_file[0])
		
		self.d.tabk.setRowCount(nrow)
		self.d.tabk.setColumnCount(ncol)
		for i,col in enumerate(data):
			for j,record in enumerate(col):
				newitem = QtWidgets.QTableWidgetItem(record)
				self.d.tabk.setItem(j, i, newitem) 

		
				
		

	def b_cancel_clicked(self):
		self.d.close()

	def showdialog_export(self):

		self.e = QtWidgets.QDialog()
		self.e.setWindowTitle("Export data")
		self.e.file_label= QtWidgets.QLabel('File:',self)
		self.e.file_field = QtWidgets.QLineEdit()
		#self.d.file_field.textChanged[str].connect(self.filechanged)
		self.e.b_select_file = QtWidgets.QPushButton(self)
		self.e.b_select_file.setText("Select file")
		self.e.b_select_file.clicked.connect(self.b_select_savefile_clicked)

		self.e.progress_label = QtWidgets.QLabel('Saving progress:',self)
		self.e.progress = QtWidgets.QProgressBar(self)
		self.e.progress.setGeometry(200, 80, 250, 20)

		self.e.b_export = QtWidgets.QPushButton(self)
		self.e.b_export.setText("Export")
		self.e.b_export.clicked.connect(self.b_export_clicked)

		self.e.b_cancel = QtWidgets.QPushButton(self)
		self.e.b_cancel.setText("Cancel")
		self.e.b_cancel.clicked.connect(self.b_cancel_export_clicked)
		layout = QtWidgets.QVBoxLayout()
		h1 = QtWidgets.QHBoxLayout()
		h2 = QtWidgets.QHBoxLayout()
		h3 = QtWidgets.QHBoxLayout()
		h1.addWidget(self.e.file_label)
		h1.addWidget(self.e.file_field)	
		h1.addWidget(self.e.b_select_file)
		h2.addWidget(self.e.progress_label)
		h2.addWidget(self.e.progress)
		h3.addWidget(self.e.b_export)
		h3.addWidget(self.e.b_cancel)
		
		layout.addLayout(h1)
		layout.addLayout(h2)
		layout.addLayout(h3)
		self.e.setLayout(layout)
		self.e.exec_()

	def showdialog_mcopt(self):
		
		self.mcopt = QtWidgets.QDialog()
		self.mcopt.setWindowTitle("Test systematic")
		self.mcopt.file_label= QtWidgets.QLabel('File:',self)
		self.mcopt.file_field = QtWidgets.QLineEdit()
		#self.d.file_field.textChanged[str].connect(self.filechanged)
		
		self.mcopt.buttons_group1 = QtWidgets.QButtonGroup()
		self.mcopt.label = QtWidgets.QLabel("Select parameters to be tested")
		#Period
		self.mcopt.b_p = QtWidgets.QCheckBox("Period")
		self.mcopt.b_p.setChecked(False)
		self.mcopt.dp1_label = QtWidgets.QLabel("range size")
		self.mcopt.dp1_field = QtWidgets.QLineEdit("0.0001")
		
		#Period change
		self.mcopt.b_dp = QtWidgets.QCheckBox("Period change")
		self.mcopt.b_dp.setChecked(False)
		self.mcopt.ddp1_label = QtWidgets.QLabel("range size")
		self.mcopt.ddp1_field = QtWidgets.QLineEdit("0.0001")
	
		#Zero time for period
		self.mcopt.b_pt0 = QtWidgets.QCheckBox("Period HJD0")
		self.mcopt.b_pt0.setChecked(False)
		self.mcopt.dpt01_label = QtWidgets.QLabel("range size")
		self.mcopt.dpt01_field = QtWidgets.QLineEdit("0.0001")
	
		

		self.mcopt.nop_label = QtWidgets.QLabel("Number of draws")
		self.mcopt.nop_field = QtWidgets.QLineEdit("500")


		self.mcopt.results = QtWidgets.QLabel('-----Results-----',self)
		self.mcopt.results.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
		self.mcopt.results.setStyleSheet("background-color: lightgreen; border: 1px solid black;")

		self.mcopt.scroll = QtWidgets.QScrollArea()
		self.mcopt.scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
		self.mcopt.scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
		self.mcopt.scroll.setWidgetResizable(True)
		self.mcopt.scroll.setWidget(self.mcopt.results)

		self.mcopt.b_run= QtWidgets.QPushButton(self)
		self.mcopt.b_run.setText("Run")
		self.mcopt.b_run.clicked.connect(self.b_run_mcopt_clicked)

		self.mcopt.b_cancel = QtWidgets.QPushButton(self)
		self.mcopt.b_cancel.setText("Cancel")
		self.mcopt.b_cancel.clicked.connect(self.b_cancel_mcopt_clicked)

		#------------------display image-----------------------------

		# image
		self.mcopt.figure= Figure(figsize=(8,6),constrained_layout=True)
		self.mcopt.canvas = FigureCanvas(self.mcopt.figure)
		self.mcopt.toolbar = NavigationToolbar(self.mcopt.canvas,self)

		self.mcopt.buttons_group2 = QtWidgets.QButtonGroup()
		self.mcopt.b_plotp = QtWidgets.QRadioButton("Plot p")
		self.mcopt.b_plotp.setChecked(True)
		self.mcopt.b_plotp.clicked.connect(self.btnstate_radio)
		self.mcopt.b_plotr = QtWidgets.QRadioButton("Plot R")
		self.mcopt.b_plotr.setChecked(False)
		self.mcopt.b_plotr.clicked.connect(self.btnstate_radio)
		self.mcopt.b_plots = QtWidgets.QRadioButton("Plot sigma")
		self.mcopt.b_plots.setChecked(False)
		self.mcopt.b_plots.clicked.connect(self.btnstate_radio)
		
		self.mcopt.buttons_group2.addButton(self.mcopt.b_plotp)
		self.mcopt.buttons_group2.addButton(self.mcopt.b_plotr)
		self.mcopt.buttons_group2.addButton(self.mcopt.b_plots)

		self.update_mcopt_plot()
		
		layout = QtWidgets.QHBoxLayout()
		left = QtWidgets.QVBoxLayout()
		right = QtWidgets.QVBoxLayout()
		h0 = QtWidgets.QHBoxLayout()
		h1 = QtWidgets.QHBoxLayout()
		h2 = QtWidgets.QHBoxLayout()
		h3 = QtWidgets.QHBoxLayout()
		h4 = QtWidgets.QHBoxLayout()
		h5 = QtWidgets.QHBoxLayout()
		h6 = QtWidgets.QHBoxLayout()
		h7 = QtWidgets.QHBoxLayout()
		
		hright1 =  QtWidgets.QHBoxLayout()
		h0.addWidget(self.mcopt.label)

		h1.addWidget(self.mcopt.b_p)
		h1.addWidget(self.mcopt.dp1_label)
		h1.addWidget(self.mcopt.dp1_field)
		

		h2.addWidget(self.mcopt.b_dp)
		h2.addWidget(self.mcopt.ddp1_label)
		h2.addWidget(self.mcopt.ddp1_field)
		

		h3.addWidget(self.mcopt.b_pt0)
		h3.addWidget(self.mcopt.dpt01_label)
		h3.addWidget(self.mcopt.dpt01_field)
		

		
		h4.addWidget(self.mcopt.nop_label)
		h4.addWidget(self.mcopt.nop_field)

		
		
		h5.addWidget(self.mcopt.scroll)

		h6.addWidget(self.mcopt.b_run)
		h6.addWidget(self.mcopt.b_cancel)

		hright1.addWidget(self.mcopt.b_plotp)
		hright1.addWidget(self.mcopt.b_plotr)
		hright1.addWidget(self.mcopt.b_plots)
		left.addLayout(h0)
		left.addLayout(h1)
		left.addLayout(h2)
		left.addLayout(h3)
		left.addLayout(h4)
		left.addLayout(h5)
		left.addLayout(h6)
		
		
		right.addWidget(self.mcopt.canvas)
		right.addWidget(self.mcopt.toolbar)
		right.addLayout(hright1)		
		layout.addLayout(left,1)
		layout.addLayout(right,3)
		self.mcopt.setLayout(layout)
		self.mcopt.exec_()

	def update_mcopt_plot(self):
		print("ble")

	def b_cancel_mcopt_clicked(self):
		self.mcopt.close()

	def btnstate_radio_mcopt(self):
		if self.mcopt.b_plotp.isChecked():
			self.mcopt_what = 0

		if self.mcopt.b_plotr.isChecked():
			self.mcopt_what = 1

		if self.mcopt.b_plots.isChecked():
			self.mcopt_what = 2

		self.update_mcopt_plot()

	def b_run_mcopt_clicked(self):
		self.mcopt_results = syst_test()
		if self.b_drs.isChecked():
			method = 0

		if self.b_lsq.isChecked():
			method = 1

		if self.b_bis.isChecked():
			method = 2

		pold = self.signals[0].period
		for i in range(int(self.mcopt.nop_field.text())):
			
			pnew = []
			if self.mcopt.b_p.isChecked():
			
				pnew.append(np.random.uniform(low=pold[0]-float(self.mcopt.dp1_field.text()),high=pold[0]+float(self.mcopt.dp1_field.text())))

			else:
				pnew.append(pold[0])

			if self.mcopt.b_dp.isChecked():
			
				pnew.append(np.random.uniform(low=pold[1]-float(self.mcopt.ddp1_field.text()),high=pold[1]+float(self.mcopt.ddp1_field.text())))

			else:
				pnew.append(pold[1])

			if self.mcopt.b_pt0.isChecked():
			
				pnew.append(np.random.uniform(low=pold[2]-float(self.mcopt.dpt01_field.text()),high=pold[2]+float(self.mcopt.dpt01_field.text())))

			else:
				pnew.append(pold[2])

			#print(pnew)
			self.signals[0].period = pnew
			self.signals[1].period = pnew
			self.signals[2].period = pnew
			self.signals[0].phasing()
			self.signals[1].phasing()
			self.signals[2].phasing()
			self.signals[0].fit()
			self.signals[1].fit()
			self.signals[2].fit()
			bw=balwan(self.name,self.signals[0],self.signals[1],self.signals[2],float(self.ebv_field.text()),ebverr=float(self.ebverr_field.text()),R_V=self.redd_v,R_K=self.redd_k,plx=float(self.plx_field.text()),plxerr=float(self.plxerr_field.text()),p=None,perr=None,irsb=self.irsb[self.wirsb()],mc=0,p0=float(self.minph_field.text()),p1=float(self.maxph_field.text()),method=method)
			
			self.mcopt_results.values.append(pnew)
			self.mcopt_results.p.append(bw.p)
			self.mcopt_results.r0.append(bw.fi0*0.5*(bw.kmkpc/bw.plx))
			self.mcopt_results.sigma.append(bw.sigma)

		self.signals[0].period = pold
		self.signals[1].period = pold
		self.signals[2].period = pold
		self.signals[0].phasing()
		self.signals[1].phasing()
		self.signals[2].phasing()
		self.signals[0].fit()
		self.signals[1].fit()
		self.signals[2].fit()

		self.update_mcopt_plot(xlabel='$\Delta P[days]$')

	def update_mcopt_plot(self,xlabel='',p=0):
		try:
			self.mcopt.figure.clf()
			self.mcopt.ax = self.mcopt.figure.add_axes([0.1,0.1,0.8,0.8])
			periods = []
			period_changes = []
			p0times = []
			for pt in self.mcopt_results.values:
				periods.append(pt[0])
				period_changes.append(pt[1])
				p0times.append(pt[2])

		
			
			
			'''self.mcopt.ax.plot(periods,np.array(self.syst_test.r0)/695700,'-',color='blue')
				self.syst.ax.set_xlabel(xlabel,fontsize=16)
				self.syst.ax.set_ylabel('$R[R_{\odot}]$',fontsize=16)
				if float(self.syst.limup_field.text()) != 0 or float(self.syst.limdown_field.text()) != 0:
					self.syst.ax.vlines(float(self.syst.limup_field.text()),np.min(self.syst_test.r0)/695700,np.max(self.syst_test.r0)/695700,linestyles='dashed',colors='black')
					self.syst.ax.vlines(float(self.syst.limdown_field.text()),np.min(self.syst_test.r0)/695700,np.max(self.syst_test.r0)/695700,linestyles='dashed',colors='black')
					a,b,sa,sb,s = metnk(self.syst_test.values,self.syst_test.r0)
					self.syst.ax.fill_between(self.syst_test.values,(a*float(self.syst.limdown_field.text())+b)/695700,(a*float(self.syst.limup_field.text())+b)/695700,facecolor='lightgreen')
					self.syst.ax.hlines(b/695700,np.min(self.syst_test.values),np.max(self.syst_test.values),linestyles='solid',colors='green')

				self.syst.ax.set_title('$\sigma _R$='+str("{:.3f}".format(np.abs(ar*float(self.syst.limdown_field.text()))/695700.))+'$R_{\odot}$',fontsize=18)'''

			self.mcopt.ax.plot(p0times,self.mcopt_results.sigma,'o',color='blue')
			self.mcopt.ax.set_xlabel(xlabel,fontsize=16)
			self.mcopt.ax.set_ylabel('$\sigma$',fontsize=16)
			#if float(self.syst.limup_field.text()) != 0 or float(self.syst.limdown_field.text()) != 0:
			#		self.syst.ax.vlines(float(self.syst.limup_field.text()),np.min(self.syst_test.sigma),np.max(self.syst_test.sigma),linestyles='dashed',colors='black')
			#		self.syst.ax.vlines(float(self.syst.limdown_field.text()),np.min(self.syst_test.sigma),np.max(self.syst_test.sigma),linestyles='dashed',colors='black')


			self.mcopt.ax.tick_params(bottom=True,top=True,left=True,right=True,labelsize=14,direction='in')
			
			self.mcopt.canvas.draw()
		except:
			pass

#test systematic errors on fitted parameters
	def showdialog_testsyst(self):
		self.syst_what = 0
		self.syst_test = syst_test()
		self.syst = QtWidgets.QDialog()
		self.syst.setWindowTitle("Test systematic")
		self.syst.file_label= QtWidgets.QLabel('File:',self)
		self.syst.file_field = QtWidgets.QLineEdit()
		#self.d.file_field.textChanged[str].connect(self.filechanged)
		
		self.syst.buttons_group1 = QtWidgets.QButtonGroup()
		self.syst.label = QtWidgets.QLabel("Select parameter to be tested")
		#Period
		self.syst.b_p = QtWidgets.QRadioButton("Period")
		self.syst.b_p.setChecked(False)
		self.syst.dp1_label = QtWidgets.QLabel("dP from")
		self.syst.dp1_field = QtWidgets.QLineEdit("-0.0001")
		self.syst.dp2_label = QtWidgets.QLabel("to")
		self.syst.dp2_field = QtWidgets.QLineEdit("0.0001")
		#V
		self.syst.b_v = QtWidgets.QRadioButton("V")
		self.syst.b_v.setChecked(True)
		self.syst.dv1_label = QtWidgets.QLabel("dV from")
		self.syst.dv1_field = QtWidgets.QLineEdit("-0.1")
		self.syst.dv2_label = QtWidgets.QLabel("to")
		self.syst.dv2_field = QtWidgets.QLineEdit("0.1")
		#K
		self.syst.b_k = QtWidgets.QRadioButton("K")
		self.syst.b_k.setChecked(False)
		self.syst.dk1_label = QtWidgets.QLabel("dK from")
		self.syst.dk1_field = QtWidgets.QLineEdit("-0.1")
		self.syst.dk2_label = QtWidgets.QLabel("to")
		self.syst.dk2_field = QtWidgets.QLineEdit("0.1")

		#ebv
		self.syst.b_ebv = QtWidgets.QRadioButton("E(B-V)")
		self.syst.b_ebv.setChecked(False)
		self.syst.debv1_label = QtWidgets.QLabel("dE(B-V) from")
		self.syst.debv1_field = QtWidgets.QLineEdit("-0.1")
		self.syst.debv2_label = QtWidgets.QLabel("to")
		self.syst.debv2_field = QtWidgets.QLineEdit("0.1")

		#irsb
		self.syst.b_irsb = QtWidgets.QRadioButton("IRSB (zero term)")
		self.syst.b_irsb.setChecked(False)
		self.syst.dirsb1_label = QtWidgets.QLabel("dIRSB from")
		self.syst.dirsb1_field = QtWidgets.QLineEdit("-0.1")
		self.syst.dirsb2_label = QtWidgets.QLabel("to")
		self.syst.dirsb2_field = QtWidgets.QLineEdit("0.1")

		#plx
		self.syst.b_plx = QtWidgets.QRadioButton("Parallax")
		self.syst.b_plx.setChecked(False)
		self.syst.dplx1_label = QtWidgets.QLabel("dplx from")
		self.syst.dplx1_field = QtWidgets.QLineEdit("-0.1")
		self.syst.dplx2_label = QtWidgets.QLabel("to")
		self.syst.dplx2_field = QtWidgets.QLineEdit("0.1")

		#K phase shift
		self.syst.b_ks = QtWidgets.QRadioButton("K phase shift")
		self.syst.b_ks.setChecked(False)
		self.syst.dks1_label = QtWidgets.QLabel("dKshift from")
		self.syst.dks1_field = QtWidgets.QLineEdit("-0.1")
		self.syst.dks2_label = QtWidgets.QLabel("to")
		self.syst.dks2_field = QtWidgets.QLineEdit("0.1")

		#RV phase shift
		self.syst.b_rvs = QtWidgets.QRadioButton("RV phase shift")
		self.syst.b_rvs.setChecked(False)
		self.syst.drvs1_label = QtWidgets.QLabel("dRVshift from")
		self.syst.drvs1_field = QtWidgets.QLineEdit("-0.1")
		self.syst.drvs2_label = QtWidgets.QLabel("to")
		self.syst.drvs2_field = QtWidgets.QLineEdit("0.1")

		self.syst.nop_label = QtWidgets.QLabel("Number of points")
		self.syst.nop_field = QtWidgets.QLineEdit("50")

		self.syst.lim_label = QtWidgets.QLabel("Plot limits")
		self.syst.limup_field = QtWidgets.QLineEdit("0")
		self.syst.limdown_field = QtWidgets.QLineEdit("0")

		self.syst.buttons_group1 =  QtWidgets.QButtonGroup()
		self.syst.buttons_group1.addButton(self.syst.b_p)
		self.syst.buttons_group1.addButton(self.syst.b_v)
		self.syst.buttons_group1.addButton(self.syst.b_k)
		self.syst.buttons_group1.addButton(self.syst.b_ebv)
		self.syst.buttons_group1.addButton(self.syst.b_irsb)
		self.syst.buttons_group1.addButton(self.syst.b_plx)
		self.syst.buttons_group1.addButton(self.syst.b_ks)
		self.syst.buttons_group1.addButton(self.syst.b_rvs)

		'''self.syst.progress_label = QtWidgets.QLabel('Progress:',self)
		self.syst.progress = QtWidgets.QProgressBar(self)
		self.syst.progress.setGeometry(200, 80, 250, 20)'''

		self.syst.results = QtWidgets.QLabel('-----Results-----',self)
		self.syst.results.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
		self.syst.results.setStyleSheet("background-color: lightgreen; border: 1px solid black;")

		self.syst.scroll = QtWidgets.QScrollArea()
		self.syst.scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
		self.syst.scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
		self.syst.scroll.setWidgetResizable(True)
		self.syst.scroll.setWidget(self.syst.results)

		self.syst.b_run= QtWidgets.QPushButton(self)
		self.syst.b_run.setText("Run")
		self.syst.b_run.clicked.connect(self.b_run_syst_clicked)

		self.syst.b_cancel = QtWidgets.QPushButton(self)
		self.syst.b_cancel.setText("Cancel")
		self.syst.b_cancel.clicked.connect(self.b_cancel_syst_clicked)

		#------------------display image-----------------------------

		# image
		self.syst.figure= Figure(figsize=(8,6),constrained_layout=True)
		#self.figure.tight_layout()
		self.syst.canvas = FigureCanvas(self.syst.figure)
		self.syst.toolbar = NavigationToolbar(self.syst.canvas,self)

		self.syst.buttons_group2 = QtWidgets.QButtonGroup()
		self.syst.b_plotp = QtWidgets.QRadioButton("Plot p")
		self.syst.b_plotp.setChecked(True)
		self.syst.b_plotp.clicked.connect(self.btnstate_radio)
		self.syst.b_plotr = QtWidgets.QRadioButton("Plot R")
		self.syst.b_plotr.setChecked(False)
		self.syst.b_plotr.clicked.connect(self.btnstate_radio)
		self.syst.b_plots = QtWidgets.QRadioButton("Plot sigma")
		self.syst.b_plots.setChecked(False)
		self.syst.b_plots.clicked.connect(self.btnstate_radio)
		
		self.syst.buttons_group2.addButton(self.syst.b_plotp)
		self.syst.buttons_group2.addButton(self.syst.b_plotr)
		self.syst.buttons_group2.addButton(self.syst.b_plots)

		self.update_syst_plot()
		
		layout = QtWidgets.QHBoxLayout()
		left = QtWidgets.QVBoxLayout()
		right = QtWidgets.QVBoxLayout()
		h0 = QtWidgets.QHBoxLayout()
		h1 = QtWidgets.QHBoxLayout()
		h2 = QtWidgets.QHBoxLayout()
		h3 = QtWidgets.QHBoxLayout()
		h4 = QtWidgets.QHBoxLayout()
		h5 = QtWidgets.QHBoxLayout()
		h6 = QtWidgets.QHBoxLayout()
		h7 = QtWidgets.QHBoxLayout()
		h8 = QtWidgets.QHBoxLayout()
		h9 = QtWidgets.QHBoxLayout()
		h10 = QtWidgets.QHBoxLayout()
		h11 = QtWidgets.QHBoxLayout()
		h12 = QtWidgets.QHBoxLayout()
		hright1 =  QtWidgets.QHBoxLayout()
		h0.addWidget(self.syst.label)

		h1.addWidget(self.syst.b_p)
		h1.addWidget(self.syst.dp1_label)
		h1.addWidget(self.syst.dp1_field)
		h1.addWidget(self.syst.dp2_label)
		h1.addWidget(self.syst.dp2_field)

		h2.addWidget(self.syst.b_v)
		h2.addWidget(self.syst.dv1_label)
		h2.addWidget(self.syst.dv1_field)
		h2.addWidget(self.syst.dv2_label)
		h2.addWidget(self.syst.dv2_field)

		h3.addWidget(self.syst.b_k)
		h3.addWidget(self.syst.dk1_label)
		h3.addWidget(self.syst.dk1_field)
		h3.addWidget(self.syst.dk2_label)
		h3.addWidget(self.syst.dk2_field)

		h4.addWidget(self.syst.b_ebv)
		h4.addWidget(self.syst.debv1_label)
		h4.addWidget(self.syst.debv1_field)
		h4.addWidget(self.syst.debv2_label)
		h4.addWidget(self.syst.debv2_field)

		h5.addWidget(self.syst.b_irsb)
		h5.addWidget(self.syst.dirsb1_label)
		h5.addWidget(self.syst.dirsb1_field)
		h5.addWidget(self.syst.dirsb2_label)
		h5.addWidget(self.syst.dirsb2_field)
		
		h6.addWidget(self.syst.b_plx)
		h6.addWidget(self.syst.dplx1_label)
		h6.addWidget(self.syst.dplx1_field)
		h6.addWidget(self.syst.dplx2_label)
		h6.addWidget(self.syst.dplx2_field)
			
		h7.addWidget(self.syst.b_ks)
		h7.addWidget(self.syst.dks1_label)
		h7.addWidget(self.syst.dks1_field)
		h7.addWidget(self.syst.dks2_label)
		h7.addWidget(self.syst.dks2_field)

		h8.addWidget(self.syst.b_rvs)
		h8.addWidget(self.syst.drvs1_label)
		h8.addWidget(self.syst.drvs1_field)
		h8.addWidget(self.syst.drvs2_label)
		h8.addWidget(self.syst.drvs2_field)
		h9.addWidget(self.syst.nop_label)
		h9.addWidget(self.syst.nop_field)

		h10.addWidget(self.syst.lim_label)
		h10.addWidget(self.syst.limdown_field)
		h10.addWidget(self.syst.limup_field)
		
		h12.addWidget(self.syst.scroll)

		h11.addWidget(self.syst.b_run)
		h11.addWidget(self.syst.b_cancel)

		hright1.addWidget(self.syst.b_plotp)
		hright1.addWidget(self.syst.b_plotr)
		hright1.addWidget(self.syst.b_plots)
		left.addLayout(h0)
		left.addLayout(h1)
		left.addLayout(h2)
		left.addLayout(h3)
		left.addLayout(h4)
		left.addLayout(h5)
		left.addLayout(h6)
		left.addLayout(h7)
		left.addLayout(h8)
		left.addLayout(h9)
		left.addLayout(h10)
		left.addLayout(h12)
		left.addLayout(h11)
		
		right.addWidget(self.syst.canvas)
		right.addWidget(self.syst.toolbar)
		right.addLayout(hright1)		
		layout.addLayout(left,1)
		layout.addLayout(right,3)
		self.syst.setLayout(layout)
		self.syst.exec_()

	def btnstate_radio(self):
		if self.syst.b_plotp.isChecked():
			self.syst_what = 0

		if self.syst.b_plotr.isChecked():
			self.syst_what = 1

		if self.syst.b_plots.isChecked():
			self.syst_what = 2

		self.update_syst_plot()

	def update_syst_plot(self,xlabel='',p=0):
		self.syst.figure.clf()
		self.syst.ax = self.syst.figure.add_axes([0.1,0.1,0.8,0.8])
		if xlabel == '':
			if self.syst.b_v.isChecked():
				xlabel='$\Delta V[mag]$'

			elif self.syst.b_p.isChecked():
				xlabel='$\Delta P[days]$'
			
			elif self.syst.b_k.isChecked():
				xlabel='$\Delta K_s[mag]$'
			elif self.syst.b_ebv.isChecked():
				xlabel='$\Delta E(B-V)[mag]$'
			elif self.syst.b_irsb.isChecked():
				xlabel='$\Delta F_V[mag]$'
			elif self.syst.b_plx.isChecked():
				xlabel='$\Delta \omega [mas]$'
			elif self.syst.b_ks.isChecked():
				xlabel='$\Delta Kphase$'
			elif self.syst.b_rvs.isChecked():
				xlabel='$\Delta RVphase$'

			else:
				xlabel = ''
		if len(self.syst_test.p) > 0:
			ap,bp,sap,sbp,sp = metnk(self.syst_test.values,self.syst_test.p)
			ar,br,sar,sbr,sr = metnk(self.syst_test.values,self.syst_test.r0)
		else:
			ap,bp = 0,0
			ar,br=0,0
		if self.syst_what == 0:
			if float(self.syst.limup_field.text()) != 0 or float(self.syst.limdown_field.text()) != 0:
				self.syst.ax.vlines(float(self.syst.limup_field.text()),np.min(self.syst_test.p),np.max(self.syst_test.p),linestyles='dashed',colors='black')
				self.syst.ax.vlines(float(self.syst.limdown_field.text()),np.min(self.syst_test.p),np.max(self.syst_test.p),linestyles='dashed',colors='black')

				minp = None
				maxp = None
				for k,pt in enumerate(self.syst_test.values):
					if minp == None and pt > float(self.syst.limdown_field.text()):
						minp = 	self.syst_test.p[k]

					if minp != None and maxp == None and pt > float(self.syst.limup_field.text()):
						maxp = 	self.syst_test.p[k-1]
				
				a,b,sa,sb,s = metnk(self.syst_test.values[5:-5],self.syst_test.p[5:-5])
				if self.oldp != 0:
					bb = b
					b = self.oldp
				else:
					b = 0
					bb=0

				#print minp,maxp
				#self.syst.ax.fill_between(self.syst_test.values,a*float(self.syst.limdown_field.text())+b,a*float(self.syst.limup_field.text())+b,facecolor='lightgreen')'''
				self.syst.ax.fill_between(self.syst_test.values,minp,maxp,facecolor='lightgreen')
				self.syst.ax.hlines(b,np.min(self.syst_test.values),np.max(self.syst_test.values),linestyles='solid',colors='green')
				self.syst.ax.plot(self.syst_test.values,self.syst_test.p,'-',color='blue')
			self.syst.ax.set_xlabel(xlabel,fontsize=16)
			self.syst.ax.set_ylabel('$p-factor$',fontsize=16)
			self.syst.ax.set_title('$\sigma _p$='+str("{:.3f}".format(np.abs(ap*float(self.syst.limdown_field.text())))),fontsize=18)
			
			
		if self.syst_what == 1:
			self.syst.ax.plot(self.syst_test.values,np.array(self.syst_test.r0)/695700,'-',color='blue')
			self.syst.ax.set_xlabel(xlabel,fontsize=16)
			self.syst.ax.set_ylabel('$R[R_{\odot}]$',fontsize=16)
			if float(self.syst.limup_field.text()) != 0 or float(self.syst.limdown_field.text()) != 0:
				self.syst.ax.vlines(float(self.syst.limup_field.text()),np.min(self.syst_test.r0)/695700,np.max(self.syst_test.r0)/695700,linestyles='dashed',colors='black')
				self.syst.ax.vlines(float(self.syst.limdown_field.text()),np.min(self.syst_test.r0)/695700,np.max(self.syst_test.r0)/695700,linestyles='dashed',colors='black')
				a,b,sa,sb,s = metnk(self.syst_test.values,self.syst_test.r0)
				self.syst.ax.fill_between(self.syst_test.values,(a*float(self.syst.limdown_field.text())+b)/695700,(a*float(self.syst.limup_field.text())+b)/695700,facecolor='lightgreen')
				self.syst.ax.hlines(b/695700,np.min(self.syst_test.values),np.max(self.syst_test.values),linestyles='solid',colors='green')

			self.syst.ax.set_title('$\sigma _R$='+str("{:.3f}".format(np.abs(ar*float(self.syst.limdown_field.text()))/695700.))+'$R_{\odot}$',fontsize=18)

		if self.syst_what == 2:
			self.syst.ax.plot(self.syst_test.values,self.syst_test.sigma,'-',color='blue')
			self.syst.ax.set_xlabel(xlabel,fontsize=16)
			self.syst.ax.set_ylabel('$\sigma$',fontsize=16)
			if float(self.syst.limup_field.text()) != 0 or float(self.syst.limdown_field.text()) != 0:
				self.syst.ax.vlines(float(self.syst.limup_field.text()),np.min(self.syst_test.sigma),np.max(self.syst_test.sigma),linestyles='dashed',colors='black')
				self.syst.ax.vlines(float(self.syst.limdown_field.text()),np.min(self.syst_test.sigma),np.max(self.syst_test.sigma),linestyles='dashed',colors='black')


		self.syst.ax.tick_params(bottom=True,top=True,left=True,right=True,labelsize=14,direction='in')
		
		self.syst.canvas.draw()
		
	def b_run_syst_clicked(self):
		self.syst_test = syst_test()
		
		#test P
		if self.syst.b_p.isChecked():

			self.syst_test.values = np.arange(float(self.syst.dp1_field.text()),float(self.syst.dp2_field.text()),(float(self.syst.dp2_field.text())-float(self.syst.dp1_field.text()))/float(self.syst.nop_field.text()))
			period = self.signals[0].period
			for i,iks in enumerate(self.syst_test.values):
				
				self.signals[0].period = self.signals[0].period + iks
				self.signals[1].period = self.signals[0].period
				self.signals[2].period = self.signals[0].period
				self.signals[0].phasing()
				self.signals[1].phasing()
				self.signals[2].phasing()
				self.signals[0].fit()
				self.signals[1].fit()
				self.signals[2].fit()
				bw=balwan(self.name,self.signals[0],self.signals[1],self.signals[2],float(self.ebv_field.text()),ebverr=float(self.ebverr_field.text()),R_V=self.redd_v,R_K=self.redd_k,plx=float(self.plx_field.text()),plxerr=float(self.plxerr_field.text()),p=None,perr=None,irsb=self.irsb[self.wirsb()],mc=0,p0=float(self.minph_field.text()),p1=float(self.maxph_field.text()))
				
				self.signals[0].period = self.signals[0].period - iks
				self.signals[1].period = self.signals[0].period
				self.signals[2].period = self.signals[0].period
				self.signals[0].phasing()
				self.signals[1].phasing()
				self.signals[2].phasing()
				self.signals[0].fit()
				self.signals[1].fit()
				self.signals[2].fit()
				self.syst_test.p.append(bw.p)
				self.syst_test.r0.append(bw.fi0*0.5*(bw.kmkpc/bw.plx))
				self.syst_test.sigma.append(bw.sigma)
			self.update_syst_plot(xlabel='$\Delta P[days]$')

		#test V
		if self.syst.b_v.isChecked():

			self.syst_test.values = np.arange(float(self.syst.dv1_field.text()),float(self.syst.dv2_field.text()),(float(self.syst.dv2_field.text())-float(self.syst.dv1_field.text()))/float(self.syst.nop_field.text()))
			for i,iks in enumerate(self.syst_test.values):
				self.signals[0].mag = np.array(self.signals[0].mag)+iks
				self.signals[0].flux = self.signals[0].toflux(self.signals[0].mag)
				self.signals[0].fit()
				bw=balwan(self.name,self.signals[0],self.signals[1],self.signals[2],float(self.ebv_field.text()),ebverr=float(self.ebverr_field.text()),R_V=self.redd_v,R_K=self.redd_k,plx=float(self.plx_field.text()),plxerr=float(self.plxerr_field.text()),p=None,perr=None,irsb=self.irsb[self.wirsb()],mc=0,p0=float(self.minph_field.text()),p1=float(self.maxph_field.text()))
				
				self.signals[0].mag = np.array(self.signals[0].mag)-iks
				self.signals[0].flux = self.signals[0].toflux(self.signals[0].mag)
				self.signals[0].fit()
				self.syst_test.p.append(bw.p)
				self.syst_test.r0.append(bw.fi0*0.5*(bw.kmkpc/bw.plx))
				self.syst_test.sigma.append(bw.sigma)
			self.update_syst_plot(xlabel='$\Delta V[mag]$')
		#test K
		if self.syst.b_k.isChecked():
			self.syst_test.values = np.arange(float(self.syst.dk1_field.text()),float(self.syst.dk2_field.text()),(float(self.syst.dk2_field.text())-float(self.syst.dk1_field.text()))/float(self.syst.nop_field.text()))
			for i,iks in enumerate(self.syst_test.values):
				self.signals[1].mag = np.array(self.signals[1].mag)+iks
				self.signals[1].flux = self.signals[1].toflux(self.signals[1].mag)
				self.signals[1].fit()
				bw=balwan(self.name,self.signals[0],self.signals[1],self.signals[2],float(self.ebv_field.text()),ebverr=float(self.ebverr_field.text()),R_V=self.redd_v,R_K=self.redd_k,plx=float(self.plx_field.text()),plxerr=float(self.plxerr_field.text()),p=None,perr=None,irsb=self.irsb[self.wirsb()],mc=0,p0=float(self.minph_field.text()),p1=float(self.maxph_field.text()))
				self.signals[1].mag = np.array(self.signals[1].mag)-iks
				self.signals[1].flux = self.signals[1].toflux(self.signals[1].mag)
				self.signals[1].fit()
				
				self.syst_test.p.append(bw.p)
				self.syst_test.r0.append(bw.fi0*0.5*(bw.kmkpc/bw.plx))
				self.syst_test.sigma.append(bw.sigma)
			self.update_syst_plot(xlabel='$\Delta K_s[mag]$')

		#test E(B-V)
		if self.syst.b_ebv.isChecked():
			#print self.wirsb()
			self.syst_test.values = np.arange(float(self.syst.debv1_field.text()),float(self.syst.debv2_field.text()),(float(self.syst.debv2_field.text())-float(self.syst.debv1_field.text()))/float(self.syst.nop_field.text()))
			for i,iks in enumerate(self.syst_test.values):
				#print iks
				bw=balwan(self.name,self.signals[0],self.signals[1],self.signals[2],float(self.ebv_field.text())+iks,ebverr=float(self.ebverr_field.text()),R_V=self.redd_v,R_K=self.redd_k,plx=float(self.plx_field.text()),plxerr=float(self.plxerr_field.text()),p=None,perr=None,irsb=self.irsb[self.wirsb()],mc=0,p0=float(self.minph_field.text()),p1=float(self.maxph_field.text()))
				self.syst_test.p.append(bw.p)
				self.syst_test.r0.append(bw.fi0*0.5*(bw.kmkpc/bw.plx))
				self.syst_test.sigma.append(bw.sigma)
			self.update_syst_plot(xlabel='$\Delta E(B-V)[mag]$')

		#test irsb
		if self.syst.b_irsb.isChecked():
			self.syst_test.values = np.arange(float(self.syst.dirsb1_field.text()),float(self.syst.dirsb2_field.text()),(float(self.syst.dirsb2_field.text())-float(self.syst.dirsb1_field.text()))/float(self.syst.nop_field.text()))
			for i,iks in enumerate(self.syst_test.values):
				bw=balwan(self.name,self.signals[0],self.signals[1],self.signals[2],float(self.ebv_field.text()),ebverr=float(self.ebverr_field.text()),R_V=self.redd_v,R_K=self.redd_k,plx=float(self.plx_field.text()),plxerr=float(self.plxerr_field.text()),p=None,perr=None,irsb=self.irsb[self.wirsb()],dirsb=iks,mc=0,p0=float(self.minph_field.text()),p1=float(self.maxph_field.text()))
				self.syst_test.p.append(bw.p)
				self.syst_test.r0.append(bw.fi0*0.5*(bw.kmkpc/bw.plx))
				self.syst_test.sigma.append(bw.sigma)
			self.update_syst_plot(xlabel='$\Delta F_V[mag/mas^2]$')

		#test plx
		if self.syst.b_plx.isChecked():
			self.syst_test.values = np.arange(float(self.syst.dplx1_field.text()),float(self.syst.dplx2_field.text()),(float(self.syst.dplx2_field.text())-float(self.syst.dplx1_field.text()))/float(self.syst.nop_field.text()))
			for i,iks in enumerate(self.syst_test.values):
				bw=balwan(self.name,self.signals[0],self.signals[1],self.signals[2],float(self.ebv_field.text()),ebverr=float(self.ebverr_field.text()),R_V=self.redd_v,R_K=self.redd_k,plx=float(self.plx_field.text())+iks,plxerr=float(self.plxerr_field.text()),p=None,perr=None,irsb=self.irsb[self.wirsb()],mc=0,p0=float(self.minph_field.text()),p1=float(self.maxph_field.text()))
				self.syst_test.p.append(bw.p)
				self.syst_test.r0.append(bw.fi0*0.5*(bw.kmkpc/bw.plx))
				self.syst_test.sigma.append(bw.sigma)
			self.update_syst_plot(xlabel='$\Delta \omega [mas]$')

		#test kshift
		if self.syst.b_ks.isChecked():
			self.syst_test.values = np.arange(float(self.syst.dks1_field.text()),float(self.syst.dks2_field.text()),(float(self.syst.dks2_field.text())-float(self.syst.dks1_field.text()))/float(self.syst.nop_field.text()))
			for i,iks in enumerate(self.syst_test.values):
				self.signals[1].ps=float(iks)
				self.signals[1].phasing()
				self.signals[1].fit()
				bw=balwan(self.name,self.signals[0],self.signals[1],self.signals[2],float(self.ebv_field.text()),ebverr=float(self.ebverr_field.text()),R_V=self.redd_v,R_K=self.redd_k,plx=float(self.plx_field.text()),plxerr=float(self.plxerr_field.text()),p=None,perr=None,irsb=self.irsb[self.wirsb()],mc=0,p0=float(self.minph_field.text()),p1=float(self.maxph_field.text()))
				self.syst_test.p.append(bw.p)
				self.syst_test.r0.append(bw.fi0*0.5*(bw.kmkpc/bw.plx))
				self.syst_test.sigma.append(bw.sigma)
			self.update_syst_plot(xlabel='$\Delta Kphase$')
			self.signals[1].ps=0.
			self.signals[1].phasing()
			self.signals[1].fit()
	
		#test rvshift
		if self.syst.b_rvs.isChecked():
			self.syst_test.values = np.arange(float(self.syst.drvs1_field.text()),float(self.syst.drvs2_field.text()),(float(self.syst.drvs2_field.text())-float(self.syst.drvs1_field.text()))/float(self.syst.nop_field.text()))
			for i,iks in enumerate(self.syst_test.values):
				self.signals[2].ps=float(iks)
				self.signals[2].phasing()
				self.signals[2].fit()
				bw=balwan(self.name,self.signals[0],self.signals[1],self.signals[2],float(self.ebv_field.text()),ebverr=float(self.ebverr_field.text()),R_V=self.redd_v,R_K=self.redd_k,plx=float(self.plx_field.text()),plxerr=float(self.plxerr_field.text()),p=None,perr=None,irsb=self.irsb[self.wirsb()],mc=0,p0=float(self.minph_field.text()),p1=float(self.maxph_field.text()))
				
				self.syst_test.p.append(bw.p)
				self.syst_test.r0.append(bw.fi0*0.5*(bw.kmkpc/bw.plx))
				self.syst_test.sigma.append(bw.sigma)
			self.update_syst_plot(xlabel='$\Delta RVphase$')
			self.signals[2].ps=0.
			self.signals[2].phasing()
			self.signals[2].fit()

		ap,bp,sap,sbp,sp = metnk(self.syst_test.values,self.syst_test.p)
		if self.oldp != 0 :
			bp = self.oldp
		ar,br,sar,sbr,sr = metnk(self.syst_test.values,self.syst_test.r0)
		self.syst.results.setText('-----Results-----\np='+str("{:.3f}".format(bp))+'-'+str("{:.3f}".format(ap*float(self.syst.limdown_field.text())))+'+'+str("{:.3f}".format(ap*float(self.syst.limup_field.text())))+'\n\nR='+str("{:.3f}".format(br/695700.))+'-'+str("{:.3f}".format((ar*float(self.syst.limdown_field.text()))/695700.))+'+'+str("{:.3f}".format((ar*float(self.syst.limup_field.text()))/695700.)))
			

	def b_cancel_syst_clicked(self):
		self.syst.close()
		
	def b_select_savefile_clicked(self):
		name = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File')
		self.e.file_field.setText(name[0])

	def b_export_clicked(self):
		maxval = len(self.signals)
		increment = 100./maxval
		completed = 0.
		self.e.progress.setValue(completed)
		plik=open(self.e.file_field.text(),'w')
		for i,t in enumerate(self.signals):
			plik.write(self.signals[i].plik+'['+str(self.signals[i].tcol)+','+str(self.signals[i].mcol)+','+str(self.signals[i].ercol)+'] '+str(self.signals[i].k)+' '+str(self.signals[i].n)+'\n')
			
			
			for j,pt in enumerate(self.signals[i].x):
				time = self.signals[i].hjd0+pt*self.signals[i].period
				plik.write(str(time)+'\t'+str(self.signals[i].y[j])+'\n')

			completed = completed + increment
			self.e.progress.setValue(completed)
		
		

		plik.close()

	def b_cancel_export_clicked(self):
		self.e.close()
		
	
	

	

	#----------------PLOT--------------------

	def onkey(self,event):
		i = self.currband
		#print event.key,event.xdata,event.ydata

		x=event.xdata
		y=event.ydata
			
		dst_prev = 1000000.
			
		whichpt = -1
		for j,iks in enumerate(self.signals[i].periods):
			igrek = self.signals[i].mags[j]
			dst_curr = (float(x)-float(iks))**2.+(float(y)-float(igrek))**2.
			if (dst_curr<dst_prev):
				dst_prev = dst_curr
				whichpt = j
				whichpt_isdel = False
				whichpt_isrej = False

		if self.pltdel_button.isChecked():
			for j,iks in enumerate(self.signals[i].del_periods):
				igrek = self.signals[i].del_mags[j]
				dst_curr = (float(x)-float(iks))**2.+(float(y)-float(igrek))**2.
				if (dst_curr<dst_prev):
					dst_prev = dst_curr
					whichpt = j
					whichpt_isdel = True
					whichpt_isrej = False


			for j,iks in enumerate(self.signals[i].rej_periods):
				igrek = self.signals[i].rej_mags[j]
				dst_curr = (float(x)-float(iks))**2.+(float(y)-float(igrek))**2.
				if (dst_curr<dst_prev):
					dst_prev = dst_curr
					whichpt = j
					whichpt_isdel = False
					whichpt_isrej = True

		if event.key == 'd':
			if not whichpt_isdel:
				if not whichpt_isrej:
					self.signals[i].rmpt(whichpt)

		elif event.key == 'a':
			if whichpt_isdel:
				self.signals[i].addpt(whichpt)
			elif whichpt_isrej:
				self.signals[i].addrejpt(whichpt)

		elif event.key == 'escape':
			text = '-----Select data point-----\n-----with doubleclick-----'
			self.pt_info.setText(text)
			self.whichpt = None
			self.whichpt_isdel = False
			self.whichpt_isrej = False

		else:
			print('Function not known. Press d to delete point or a to add point')

		self.update_plot(i)

	def onclick(self,event):
		i = self.currband
		if event.dblclick:
			x=event.xdata
			y=event.ydata
			print(event.inaxes)
			
			#ra,dec = image2wcs(self.image_name,[x],[y])
			
			dst_prev = 1000000.
			
			self.whichpt = -1
			for j,iks in enumerate(self.signals[i].periods):
				igrek = self.signals[i].mags[j]
				dst_curr = (float(x)-float(iks))**2.+(float(y)-float(igrek))**2.
				if (dst_curr<dst_prev):
					dst_prev = dst_curr
					self.whichpt = j
					self.whichpt_isdel = False
					self.whichpt_isrej = False

			if self.pltdel_button.isChecked():
				for j,iks in enumerate(self.signals[i].del_periods):
					igrek = self.signals[i].del_mags[j]
					dst_curr = (float(x)-float(iks))**2.+(float(y)-float(igrek))**2.
					if (dst_curr<dst_prev):
						dst_prev = dst_curr
						self.whichpt = j
						self.whichpt_isdel = True
						self.whichpt_isrej = False

			
				for j,iks in enumerate(self.signals[i].rej_periods):
					igrek = self.signals[i].rej_mags[j]
					dst_curr = (float(x)-float(iks))**2.+(float(y)-float(igrek))**2.
					if (dst_curr<dst_prev):
						dst_prev = dst_curr
						self.whichpt = j
						self.whichpt_isdel = False
						self.whichpt_isrej = True

		try:
			if self.whichpt_isdel:
				text = 'Point deleted manually\n\nname: '+str(self.signals[i].del_names[self.whichpt])+'\nx: '+str(self.signals[i].del_periods[self.whichpt])+'\ny: '+str(self.signals[i].del_mags[self.whichpt])+'\ny_error: '+str(self.signals[i].del_magserr[self.whichpt])+'\nra: '+str(self.signals[i].del_ras[self.whichpt])+'\ndec: '+str(self.signals[i].del_decs[self.whichpt])+'\n'
			elif self.whichpt_isrej:
				text = 'Point rejected with 3-sigma clipping\n\nname: '+str(self.signals[i].rej_names[self.whichpt])+'\nlogP: '+str(self.signals[i].rej_periods[self.whichpt])+'\nmag: '+str(self.signals[i].rej_mags[self.whichpt])+'\nmag_error: '+str(self.signals[i].rej_magserr[self.whichpt])+'\nra: '+str(self.signals[i].rej_ras[self.whichpt])+'\ndec: '+str(self.signals[i].rej_decs[self.whichpt])+'\n'
			else:
				text = 'Accepted point\n\nname: '+str(self.signals[i].names[self.whichpt])+'\nx: '+str(self.signals[i].periods[self.whichpt])+'\ny: '+str(self.signals[i].mags[self.whichpt])+'\ny_error: '+str(self.signals[i].magserr[self.whichpt])+'\nra: '+str(self.signals[i].ras[self.whichpt])+'\ndec: '+str(self.signals[i].decs[self.whichpt])+'\n'

			self.pt_info.setText(text)
		except:
			pass
		self.update_plot(i)

	def update_plot(self,i):
		self.currband = i
		self.plot(i)

	def plot(self,i):
		self.figure.clf()
		try:
			gs = gridspec.GridSpec(ncols=5, nrows=3,figure=self.figure)
		except:		
			gs = gridspec.GridSpec(ncols=5, nrows=3,wspace=0.5,hspace=0.5)

		phasename = '$\phi$'
		self.f_v = self.figure.add_subplot(gs[0,0])
		self.f_v.set_ylabel('$V$ $[mag]$',fontsize=14)
		self.f_v.set_xlabel(phasename,fontsize=14)
		self.f_v.set_xlim(-0.1,1.1)

		self.f_k = self.figure.add_subplot(gs[1,0],sharex=self.f_v)
		self.f_k.set_ylabel('$K$ $[mag]$',fontsize=14)
		self.f_k.set_xlabel(phasename,fontsize=14)
		self.f_k.set_xlim(-0.1,1.1)

		self.f_vk = self.figure.add_subplot(gs[2,0],sharex=self.f_v)
		self.f_vk.set_ylabel('$(V-K)_0$ $[mag]$',fontsize=14)
		self.f_vk.set_xlabel(phasename,fontsize=14)
		self.f_vk.set_xlim(-0.1,1.1)

		self.f_rv = self.figure.add_subplot(gs[0,1])
		self.f_rv.set_ylabel('$V_r-V_{\gamma}$ $[km/s]$',fontsize=14)
		self.f_rv.set_xlabel(phasename,fontsize=14)	
		self.f_rv.set_xlim(-0.1,1.1)

		self.f_dr = self.figure.add_subplot(gs[1,1],sharex=self.f_rv)
		self.f_dr.set_ylabel('$\Delta R\'$ $[km]$',fontsize=14)
		self.f_dr.set_xlabel(phasename,fontsize=14)
		self.f_dr.set_xlim(-0.1,1.1)

		self.f_fi = self.figure.add_subplot(gs[2,1],sharex=self.f_rv)
		self.f_fi.set_ylabel('$\Theta$ $[mas]$',fontsize=14)
		self.f_fi.set_xlabel(phasename,fontsize=14)
		self.f_fi.set_xlim(-0.1,1.1)

		self.f_bw = self.figure.add_subplot(gs[:,2:])
		self.f_bw.set_title('$BW$',fontsize=14)
		self.f_bw.set_ylabel('$\Theta$ $[mas]$',fontsize=14)
		self.f_bw.set_xlabel('$2x \omega x \Delta R \' / 1AU$ $[mas]$',fontsize=14)
		
		self.f_rv.tick_params(bottom=True,top=True,left=True,right=True,labelsize=12,direction='in')
		self.f_v.tick_params(bottom=True,top=True,left=True,right=True,labelsize=12,direction='in')
		self.f_k.tick_params(bottom=True,top=True,left=True,right=True,labelsize=12,direction='in')
		self.f_dr.tick_params(bottom=True,top=True,left=True,right=True,labelsize=12,direction='in')
		self.f_vk.tick_params(bottom=True,top=True,left=True,right=True,labelsize=12,direction='in')
		self.f_fi.tick_params(bottom=True,top=True,left=True,right=True,labelsize=12,direction='in')
		self.f_bw.tick_params(bottom=True,top=True,left=True,right=True,labelsize=12,direction='in')
		self.canvas.setFocusPolicy( QtCore.Qt.ClickFocus )
		self.canvas.setFocus()
		self.canvas.mpl_connect('button_press_event', self.onclick)
		self.canvas.mpl_connect('key_press_event', self.onkey)

			
		if self.bw != None:
			if self.mc_button.isChecked():
				self.f_v.plot(self.signals[0].x,self.signals[0].y,'-',color='green')
				#self.f_v.plot(self.signals[0].x,np.add(self.signals[0].y,self.signals[0].yerr),'-',color='lightgreen')
				#self.f_v.plot(self.signals[0].x,np.subtract(self.signals[0].y,self.signals[0].yerr),'-',color='lightgreen')
				self.f_k.plot(self.signals[1].x,self.signals[1].y,'-',color='green')
				#self.f_k.plot(self.signals[1].x,np.add(self.signals[1].y,self.signals[1].yerr),'-',color='lightgreen')
				#self.f_k.plot(self.signals[1].x,np.subtract(self.signals[1].y,self.signals[1].yerr),'-',color='lightgreen')
				self.f_rv.plot(self.signals[2].x,self.signals[2].y,'-',color='green')
				#self.f_rv.plot(self.signals[2].x,np.add(self.signals[2].y,self.signals[2].yerr),'-',color='lightgreen')
				#self.f_rv.plot(self.signals[2].x,np.subtract(self.signals[2].y,self.signals[2].yerr),'-',color='lightgreen')
				self.f_v.errorbar(self.signals[0].phase,self.signals[0].mag,yerr=self.signals[0].magerr,fmt='o',color='blue',markersize=4)
				self.f_k.errorbar(self.signals[1].phase,self.signals[1].mag,yerr=self.signals[1].magerr,fmt='o',color='blue',markersize=4)
				self.f_rv.errorbar(self.signals[2].phase,self.signals[2].mag,yerr=self.signals[2].magerr,fmt='o',color='blue',markersize=4)

				self.f_v.plot(self.signals[0].x,self.signals[0].y,'-',color='black')
				self.f_k.plot(self.signals[1].x,self.signals[1].y,'-',color='black')
				self.f_rv.plot(self.signals[2].x,self.signals[2].y,'-',color='black')
				if len(self.bw.vk_err) > 0:
					self.f_vk.errorbar(self.bw.k.phase,self.bw.vk,yerr=self.bw.vk_err,fmt='o',color='blue',markersize=4)	
					self.f_dr.errorbar(self.bw.phase,self.bw.dr,yerr=self.bw.dr_err,fmt='o',color='blue',markersize=4)
					self.f_fi.errorbar(self.bw.phase,(180.*3600000./np.pi)*self.bw.fi,yerr=(180.*3600000./np.pi)*self.bw.fi_err,fmt='o',color='blue',markersize=4)
					self.f_dr.errorbar(self.bw.phase_del,self.bw.dr_del+self.bw.mean_rad,yerr=self.bw.dr_err_del,fmt='o',color='red',markersize=4)
					self.f_fi.errorbar(self.bw.phase_del,(180.*3600000./np.pi)*self.bw.fi_del,yerr=(180.*3600000./np.pi)*self.bw.fi_err_del,fmt='o',color='red',markersize=4)
				else:
					self.f_vk.plot(self.bw.k.phase,self.bw.vk,'o',color='blue',markersize=4)
					self.f_dr.plot(self.bw.phase,self.bw.dr,'o',color='blue',markersize=4)
					self.f_fi.plot(self.bw.phase,(180.*3600000./np.pi)*self.bw.fi,'o',color='blue',markersize=4)
					self.f_dr.plot(self.bw.phase_del,self.bw.dr_del,'o',color='red',markersize=4)
					self.f_fi.plot(self.bw.phase_del,(180.*3600000./np.pi)*self.bw.fi_del,'o',color='red',markersize=4)

				self.f_vk.plot(self.bw.phase_plot,self.bw.vk_plot,'-',color='black')
				self.f_dr.plot(self.bw.phase_plot,self.bw.dr_plot,'-',color='black')
				ind_min = np.argmin(self.bw.dr_plot)
				val_min = np.min(self.bw.dr_plot)
				ind_max = np.argmax(self.bw.dr_plot)
				val_max = np.max(self.bw.dr_plot)
				self.f_dr.vlines(self.bw.phase_plot[ind_min],val_min,val_max,linestyles='dotted',color='red')
				self.f_dr.vlines(self.bw.phase_plot[ind_max],val_min,val_max,linestyles='dotted',color='red')
				#self.f_dr.hlines(self.bw.mean_rad,0,1,linestyles='dashed',color='green')

				self.f_fi.plot(self.bw.phase_plot,(180.*3600000./np.pi)*self.bw.fi_plot,'-',color='black')
				self.f_fi.plot(self.bw.phase_plot,(180.*3600000./np.pi)*self.bw.fi_plot2,'--',color='green')
				val_min = np.min(self.bw.fi_plot)
				val_max = np.max(self.bw.fi_plot)
				self.f_fi.vlines(self.bw.phase_plot[ind_min],(180.*3600000./np.pi)*val_min,(180.*3600000./np.pi)*val_max,linestyles='dotted',color='red')
				self.f_fi.vlines(self.bw.phase_plot[ind_max],(180.*3600000./np.pi)*val_min,(180.*3600000./np.pi)*val_max,linestyles='dotted',color='red')
				#print self.bw.dr_err
				if len(self.bw.fi_err) > 0: #uncomment this to plot errorbars
					self.f_bw.errorbar((180.*3600000./np.pi)*np.multiply(self.bw.dr,2.*self.bw.plx/self.bw.kmkpc),(180.*3600000./np.pi)*np.array(self.bw.fi),xerr=(180.*3600000./np.pi)*np.multiply(self.bw.dr_err,2.*self.bw.plx/self.bw.kmkpc),yerr=(180.*3600000./np.pi)*np.array(self.bw.fi_err),fmt='o',color='blue')
				#	self.f_bw.errorbar((180.*3600000./np.pi)*np.multiply(self.bw.dr_del,2.*self.bw.plx/self.bw.kmkpc),(180.*3600000./np.pi)*np.array(self.bw.fi_del),xerr=(180.*3600000./np.pi)*np.multiply(self.bw.dr_del,2.*self.bw.plx/self.bw.kmkpc),yerr=(180.*3600000./np.pi)*np.array(self.bw.fi_err_del),fmt='o',color='red')
				#else:
				self.f_bw.plot((180.*3600000./np.pi)*np.multiply(self.bw.dr,2.*self.bw.plx/self.bw.kmkpc),(180.*3600000./np.pi)*np.array(self.bw.fi),'o',color='blue')
				self.f_bw.plot((180.*3600000./np.pi)*np.multiply(self.bw.dr_del,2.*self.bw.plx/self.bw.kmkpc),(180.*3600000./np.pi)*np.array(self.bw.fi_del),'o',color='red')

				self.f_bw.plot((180.*3600000./np.pi)*np.multiply(self.bw.dr,2.*self.bw.plx/self.bw.kmkpc),(180.*3600000./np.pi)*np.add(np.multiply(self.bw.p,np.multiply(self.bw.dr,2.*self.bw.plx/self.bw.kmkpc)),np.array(self.bw.fi0)),'-')

				
				self.f_bw.text((180.*3600000./np.pi)*np.mean(np.multiply(self.bw.dr,2.*self.bw.plx/self.bw.kmkpc)),(180.*3600000./np.pi)*np.min(self.bw.fi),'%s\n p-factor= %s $\pm$ %s\n R= %s $\pm$ %s $R_{\odot}$' % (self.bw.name,str("{:.2f}".format(self.bw.p)), str("{:.2f}".format(self.bw.sigmap)),str("{:.2f}".format(self.bw.fi0*0.5*(self.bw.kmkpc/self.bw.plx)/695700)),str("{:.2f}".format(self.bw.sigmafi0*0.5*(self.bw.kmkpc/self.bw.plx)/695700))),fontstyle='italic',fontsize=16)
				out = open('balwan_out.dat','w')
				try:
					for i,punkt in enumerate(self.bw.dr):
						out.write(str((180.*3600000./np.pi)*np.multiply(self.bw.dr[i],2.*self.bw.plx/self.bw.kmkpc))+'\t'+str((180.*3600000./np.pi)*np.array(self.bw.fi[i]))+'\n')
	
					out.close()
				except:
					pass
			else:
				self.f_v.errorbar(self.signals[0].phase,self.signals[0].mag,yerr=self.signals[0].magerr,fmt='o',color='blue',markersize=4)
				self.f_k.errorbar(self.signals[1].phase,self.signals[1].mag,yerr=self.signals[1].magerr,fmt='o',color='blue',markersize=4)
				self.f_rv.errorbar(self.signals[2].phase,self.signals[2].mag,yerr=self.signals[2].magerr,fmt='o',color='blue',markersize=4)

				
				self.f_v.errorbar(self.signals[0].phase,self.signals[0].mag,yerr=self.signals[0].magerr,fmt='o',color='blue',markersize=4)
				self.f_k.errorbar(self.signals[1].phase,self.signals[1].mag,yerr=self.signals[1].magerr,fmt='o',color='blue',markersize=4)
				self.f_rv.errorbar(self.signals[2].phase,self.signals[2].mag,yerr=self.signals[2].magerr,fmt='o',color='blue',markersize=4)
				self.f_v.plot(self.signals[0].x,self.signals[0].y,'-',color='black')
				self.f_k.plot(self.signals[1].x,self.signals[1].y,'-',color='black')
				self.f_rv.plot(self.signals[2].x,self.signals[2].y,'-',color='black')

				self.f_v.plot(self.signals[0].x,self.signals[0].y,'-',color='black')
				self.f_k.plot(self.signals[1].x,self.signals[1].y,'-',color='black')
				self.f_rv.plot(self.signals[2].x,self.signals[2].y,'-',color='black')

				self.f_vk.plot(self.bw.k.phase,self.bw.vk,'o',color='blue',markersize=4)	
				self.f_vk.plot(self.bw.phase_plot,self.bw.vk_plot,'-',color='black')
				self.f_dr.plot(self.bw.phase,self.bw.dr,'o',color='blue',markersize=4)
				self.f_dr.plot(self.bw.phase_del,self.bw.dr_del,'o',color='red',markersize=4)
				self.f_dr.plot(self.bw.phase_plot,self.bw.dr_plot,'-',color='black')
				ind_min = np.argmin(self.bw.dr_plot)
				val_min = np.min(self.bw.dr_plot)
				ind_max = np.argmax(self.bw.dr_plot)
				val_max = np.max(self.bw.dr_plot)
				self.f_dr.vlines(self.bw.phase_plot[ind_min],val_min,val_max,linestyles='dotted',color='red')
				self.f_dr.vlines(self.bw.phase_plot[ind_max],val_min,val_max,linestyles='dotted',color='red')
				#self.f_dr.hlines(0,0,1,linestyles='dashed',color='green')
				self.f_fi.plot(self.bw.phase,(180.*3600000./np.pi)*self.bw.fi,'o',color='blue',markersize=4)
				self.f_fi.plot(self.bw.phase_del,(180.*3600000./np.pi)*self.bw.fi_del,'o',color='red',markersize=4)
				self.f_fi.plot(self.bw.phase_plot,(180.*3600000./np.pi)*self.bw.fi_plot,'-',color='black')
				self.f_fi.plot(self.bw.phase_plot,(180.*3600000./np.pi)*self.bw.fi_plot2,'--',color='green')
				val_min = np.min(self.bw.fi_plot)
				val_max = np.max(self.bw.fi_plot)
				self.f_fi.vlines(self.bw.phase_plot[ind_min],(180.*3600000./np.pi)*val_min,(180.*3600000./np.pi)*val_max,linestyles='dotted',color='red')
				self.f_fi.vlines(self.bw.phase_plot[ind_max],(180.*3600000./np.pi)*val_min,(180.*3600000./np.pi)*val_max,linestyles='dotted',color='red')
				self.f_bw.plot((180.*3600000./np.pi)*np.multiply(self.bw.dr,2.*self.bw.plx/self.bw.kmkpc),(180.*3600000./np.pi)*np.array(self.bw.fi),'o',color='blue')
				self.f_bw.plot((180.*3600000./np.pi)*np.multiply(self.bw.dr_del,2.*self.bw.plx/self.bw.kmkpc),(180.*3600000./np.pi)*np.array(self.bw.fi_del),'o',color='red')	
				self.f_bw.plot((180.*3600000./np.pi)*np.multiply(self.bw.dr,2.*self.bw.plx/self.bw.kmkpc),(180.*3600000./np.pi)*np.add(np.multiply(self.bw.p,np.multiply(self.bw.dr,2.*self.bw.plx/self.bw.kmkpc)),self.bw.fi0),'-')
				if self.b_plx.isChecked():
					self.f_bw.text((180.*3600000./np.pi)*np.mean(np.multiply(self.bw.dr,2.*self.bw.plx/self.bw.kmkpc)),(180.*3600000./np.pi)*np.min(self.bw.fi),'%s\n p-factor= %s $\pm$ %s\n R= %s $\pm$ %s $R_{\odot}$' % (self.bw.name,str("{:.2f}".format(self.bw.p)), str("{:.2f}".format(self.bw.sigmap)),str("{:.2f}".format(self.bw.fi0*0.5*(self.bw.kmkpc/self.bw.plx)/695700)),str("{:.2f}".format(self.bw.sigmafi0*0.5*(self.bw.kmkpc/self.bw.plx)/695700))),fontstyle='italic',fontsize=16)
				#else:

		elif len(self.signals) == 3:
			self.f_v.errorbar(self.signals[0].phase,self.signals[0].mag,yerr=self.signals[0].magerr,fmt='o',color='blue',markersize=4)
			self.f_k.errorbar(self.signals[1].phase,self.signals[1].mag,yerr=self.signals[1].magerr,fmt='o',color='blue',markersize=4)
			self.f_rv.errorbar(self.signals[2].phase,self.signals[2].mag,yerr=self.signals[2].magerr,fmt='o',color='blue',markersize=4)

			self.f_v.plot(self.signals[0].x,self.signals[0].y,'-',color='black')
			self.f_k.plot(self.signals[1].x,self.signals[1].y,'-',color='black')
			self.f_rv.plot(self.signals[2].x,self.signals[2].y,'-',color='black')

		self.canvas.mpl_connect('button_press_event', self.onclick)

		
		self.f_v.invert_yaxis()
		self.f_k.invert_yaxis()
		self.f_vk.invert_yaxis()
		self.canvas.draw()	

'''	def plot(self,i):
		self.figure.clf()
		
		gs = gridspec.GridSpec(ncols=5, nrows=3, figure=self.figure)

		self.f_v = self.figure.add_subplot(gs[0,0])
		self.f_v.set_ylabel('$V$ $[mag]$')
		self.f_v.set_xlabel('$phase$')
		self.f_v.set_xlim(-0.1,1.1)

		self.f_k = self.figure.add_subplot(gs[1,0],sharex=self.f_v)
		self.f_k.set_ylabel('$K$ $[mag]$')
		self.f_k.set_xlabel('$phase$')
		self.f_k.set_xlim(-0.1,1.1)

		self.f_vk = self.figure.add_subplot(gs[2,0],sharex=self.f_v)
		self.f_vk.set_ylabel('$V-K$ $[mag]$')
		self.f_vk.set_xlabel('$phase$')
		self.f_vk.set_xlim(-0.1,1.1)

		self.f_rv = self.figure.add_subplot(gs[0,1])
		self.f_rv.set_ylabel('$V_r$ $[km/s]$')
		self.f_rv.set_xlabel('$phase$')	
		self.f_rv.set_xlim(-0.1,1.1)

		self.f_dr = self.figure.add_subplot(gs[1,1],sharex=self.f_rv)
		self.f_dr.set_ylabel('$\int V_r dt$ $[km]$')
		self.f_dr.set_xlabel('$phase$')
		self.f_dr.set_xlim(-0.1,1.1)

		self.f_fi = self.figure.add_subplot(gs[2,1],sharex=self.f_rv)
		self.f_fi.set_ylabel('$\Theta$ $[rad]$')
		self.f_fi.set_xlabel('$phase$')
		self.f_fi.set_xlim(-0.1,1.1)

		self.f_bw = self.figure.add_subplot(gs[:,2:])
		self.f_bw.set_title('$BW$')
		self.f_bw.set_ylabel('$\Theta$ $[mas]$')
		self.f_bw.set_xlabel('$2x \omega x \int V_r dt / kpc$ $[mas]$')
		
		self.f_rv.tick_params(bottom=True,top=True,left=True,right=True,labelsize=12,direction='in')
		self.f_v.tick_params(bottom=True,top=True,left=True,right=True,labelsize=12,direction='in')
		self.f_k.tick_params(bottom=True,top=True,left=True,right=True,labelsize=12,direction='in')
		self.f_dr.tick_params(bottom=True,top=True,left=True,right=True,labelsize=12,direction='in')
		self.f_vk.tick_params(bottom=True,top=True,left=True,right=True,labelsize=12,direction='in')
		self.f_fi.tick_params(bottom=True,top=True,left=True,right=True,labelsize=12,direction='in')
		self.f_bw.tick_params(bottom=True,top=True,left=True,right=True,labelsize=12,direction='in')
		self.canvas.setFocusPolicy( QtCore.Qt.ClickFocus )
		self.canvas.setFocus()
		self.canvas.mpl_connect('button_press_event', self.onclick)
		self.canvas.mpl_connect('key_press_event', self.onkey)

			
		if self.bw != None:
			if self.mc_button.isChecked():
				self.f_v.plot(self.signals[0].x,self.signals[0].y,'-',color='green')
				#self.f_v.plot(self.signals[0].x,np.add(self.signals[0].y,self.signals[0].yerr),'-',color='lightgreen')
				#self.f_v.plot(self.signals[0].x,np.subtract(self.signals[0].y,self.signals[0].yerr),'-',color='lightgreen')
				self.f_k.plot(self.signals[1].x,self.signals[1].y,'-',color='green')
				#self.f_k.plot(self.signals[1].x,np.add(self.signals[1].y,self.signals[1].yerr),'-',color='lightgreen')
				#self.f_k.plot(self.signals[1].x,np.subtract(self.signals[1].y,self.signals[1].yerr),'-',color='lightgreen')
				self.f_rv.plot(self.signals[2].x,self.signals[2].y,'-',color='green')
				#self.f_rv.plot(self.signals[2].x,np.add(self.signals[2].y,self.signals[2].yerr),'-',color='lightgreen')
				#self.f_rv.plot(self.signals[2].x,np.subtract(self.signals[2].y,self.signals[2].yerr),'-',color='lightgreen')
				self.f_v.errorbar(self.signals[0].phase,self.signals[0].mag,yerr=self.signals[0].magerr,fmt='o',color='blue',markersize=4)
				self.f_k.errorbar(self.signals[1].phase,self.signals[1].mag,yerr=self.signals[1].magerr,fmt='o',color='blue',markersize=4)
				self.f_rv.errorbar(self.signals[2].phase,self.signals[2].mag,yerr=self.signals[2].magerr,fmt='o',color='blue',markersize=4)

				self.f_v.plot(self.signals[0].x,self.signals[0].y,'-',color='black')
				self.f_k.plot(self.signals[1].x,self.signals[1].y,'-',color='black')
				self.f_rv.plot(self.signals[2].x,self.signals[2].y,'-',color='black')
				if len(self.bw.vk_err) > 0:
					self.f_vk.errorbar(self.bw.k.phase,self.bw.vk,yerr=self.bw.vk_err,fmt='o',color='blue',markersize=4)	
					self.f_dr.errorbar(self.bw.k.phase,self.bw.dr,yerr=self.bw.dr_err,fmt='o',color='blue',markersize=4)
					self.f_fi.errorbar(self.bw.k.phase,self.bw.fi,yerr=self.bw.fi_err,fmt='o',color='blue',markersize=4)
				else:
					self.f_vk.plot(self.bw.k.phase,self.bw.vk,'o',color='blue',markersize=4)
					self.f_dr.plot(self.bw.k.phase,self.bw.dr,'o',color='blue',markersize=4)
					self.f_fi.plot(self.bw.k.phase,self.bw.fi,'o',color='blue',markersize=4)

				self.f_vk.plot(self.bw.phase,self.bw.vk_plot,'-',color='black')
				self.f_dr.plot(self.bw.phase,self.bw.dr_plot,'-',color='black')
				ind_min = np.argmin(self.bw.dr_plot)
				val_min = np.min(self.bw.dr_plot)
				ind_max = np.argmax(self.bw.dr_plot)
				val_max = np.max(self.bw.dr_plot)
				self.f_dr.vlines(self.bw.phase[ind_min],val_min,val_max,linestyles='dotted',color='red')
				self.f_dr.vlines(self.bw.phase[ind_max],val_min,val_max,linestyles='dotted',color='red')

				self.f_fi.plot(self.bw.phase,self.bw.fi_plot,'-',color='black')
				#self.f_fi.plot(self.bw.phase,self.bw.fi_plot2,'-',color='green')
				val_min = np.min(self.bw.fi_plot)
				val_max = np.max(self.bw.fi_plot)
				self.f_fi.vlines(self.bw.phase[ind_min],val_min,val_max,linestyles='dotted',color='red')
				self.f_fi.vlines(self.bw.phase[ind_max],val_min,val_max,linestyles='dotted',color='red')
				if len(self.bw.fi_err) > 0:
					self.f_bw.errorbar((180.*3600000./np.pi)*np.multiply(self.bw.dr,2.*self.bw.plx/self.bw.kmkpc),(180.*3600000./np.pi)*self.bw.fi,yerr=(180.*3600000./np.pi)*self.bw.fi_err,fmt='o',color='blue')
				else:
					self.f_bw.plot((180.*3600000./np.pi)*np.multiply(self.bw.dr,2.*self.bw.plx/self.bw.kmkpc),(180.*3600000./np.pi)*self.bw.fi,'o',color='blue')
				self.f_bw.plot((180.*3600000./np.pi)*np.multiply(self.bw.dr,2.*self.bw.plx/self.bw.kmkpc),np.add(np.multiply(self.bw.p,np.multiply(self.bw.dr,2.*self.bw.plx/self.bw.kmkpc)),(180.*3600000./np.pi)*self.bw.fi0),'-')
				self.f_bw.text(np.mean(np.multiply(self.bw.dr,2.*self.bw.plx/self.bw.kmkpc)),np.min(self.bw.fi),'%s\n p-factor= %s $\pm$ %s' % (self.bw.name,str("{:.3f}".format(self.bw.p)),str("{:.3f}".format(self.bw.sigmap))),fontstyle='italic')
			else:
				self.f_v.errorbar(self.signals[0].phase,self.signals[0].mag,yerr=self.signals[0].magerr,fmt='o',color='blue',markersize=4)
				self.f_k.errorbar(self.signals[1].phase,self.signals[1].mag,yerr=self.signals[1].magerr,fmt='o',color='blue',markersize=4)
				self.f_rv.errorbar(self.signals[2].phase,self.signals[2].mag,yerr=self.signals[2].magerr,fmt='o',color='blue',markersize=4)

				self.f_v.plot(self.signals[0].x,self.signals[0].y,'-',color='black')
				self.f_k.plot(self.signals[1].x,self.signals[1].y,'-',color='black')
				self.f_rv.plot(self.signals[2].x,self.signals[2].y,'-',color='black')
				self.f_v.errorbar(self.signals[0].phase,self.signals[0].mag,yerr=self.signals[0].magerr,fmt='o',color='blue',markersize=4)
				self.f_k.errorbar(self.signals[1].phase,self.signals[1].mag,yerr=self.signals[1].magerr,fmt='o',color='blue',markersize=4)
				self.f_rv.errorbar(self.signals[2].phase,self.signals[2].mag,yerr=self.signals[2].magerr,fmt='o',color='blue',markersize=4)

				self.f_v.plot(self.signals[0].x,self.signals[0].y,'-',color='black')
				self.f_k.plot(self.signals[1].x,self.signals[1].y,'-',color='black')
				self.f_rv.plot(self.signals[2].x,self.signals[2].y,'-',color='black')

				self.f_vk.plot(self.bw.k.phase,self.bw.vk,'o',color='blue',markersize=4)	
				self.f_vk.plot(self.bw.phase,self.bw.vk_plot,'-',color='black')
				self.f_dr.plot(self.bw.k.phase,self.bw.dr,'o',color='blue',markersize=4)
				self.f_dr.plot(self.bw.phase,self.bw.dr_plot,'-',color='black')
				ind_min = np.argmin(self.bw.dr_plot)
				val_min = np.min(self.bw.dr_plot)
				ind_max = np.argmax(self.bw.dr_plot)
				val_max = np.max(self.bw.dr_plot)
				self.f_dr.vlines(self.bw.phase[ind_min],val_min,val_max,linestyles='dotted',color='red')
				self.f_dr.vlines(self.bw.phase[ind_max],val_min,val_max,linestyles='dotted',color='red')
				self.f_fi.plot(self.bw.k.phase,self.bw.fi,'o',color='blue',markersize=4)
				self.f_fi.plot(self.bw.phase,self.bw.fi_plot,'-',color='black')
				self.f_fi.plot(self.bw.phase,self.bw.fi_plot2,'-',color='green')
				val_min = np.min(self.bw.fi_plot)
				val_max = np.max(self.bw.fi_plot)
				self.f_fi.vlines(self.bw.phase[ind_min],val_min,val_max,linestyles='dotted',color='red')
				self.f_fi.vlines(self.bw.phase[ind_max],val_min,val_max,linestyles='dotted',color='red')
				self.f_bw.plot((180.*3600000./np.pi)*np.multiply(self.bw.dr,2.*self.bw.plx/self.bw.kmkpc),(180.*3600000./np.pi)*self.bw.fi,'o',color='blue')	
				self.f_bw.plot((180.*3600000./np.pi)*np.multiply(self.bw.dr,2.*self.bw.plx/self.bw.kmkpc),(180.*3600000./np.pi)*np.add(np.multiply(self.bw.p,np.multiply(self.bw.dr,2.*self.bw.plx/self.bw.kmkpc)),self.bw.fi0),'-')
				self.f_bw.text(np.mean(np.multiply(self.bw.dr,2.*self.bw.plx/self.bw.kmkpc)),np.min(self.bw.fi),'%s\n p-factor= %s $\pm$ %s' % (self.bw.name,str("{:.3f}".format(self.bw.p)),str("{:.3f}".format(self.bw.sigmap))),fontstyle='italic')
		elif len(self.signals) == 3:
			self.f_v.errorbar(self.signals[0].phase,self.signals[0].mag,yerr=self.signals[0].magerr,fmt='o',color='blue',markersize=4)
			self.f_k.errorbar(self.signals[1].phase,self.signals[1].mag,yerr=self.signals[1].magerr,fmt='o',color='blue',markersize=4)
			self.f_rv.errorbar(self.signals[2].phase,self.signals[2].mag,yerr=self.signals[2].magerr,fmt='o',color='blue',markersize=4)

			self.f_v.plot(self.signals[0].x,self.signals[0].y,'-',color='black')
			self.f_k.plot(self.signals[1].x,self.signals[1].y,'-',color='black')
			self.f_rv.plot(self.signals[2].x,self.signals[2].y,'-',color='black')

		#self.canvas.mpl_connect('button_press_event', self.onclick)

		
		self.f_v.invert_yaxis()
		self.f_k.invert_yaxis()
		self.f_vk.invert_yaxis()
		self.canvas.draw()'''

def main():

	app=QtWidgets.QApplication(sys.argv)

	ex=display()
	sys.exit(app.exec_())


if __name__=='__main__':

	main()
