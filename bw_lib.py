#! /usr/bin/env python
#
#======================================
#  Marek Gorski, Bartek Zgirski
#  05.02.2018
#======================================
#

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import Akima1DInterpolator


#fazowanie
def phase(Ojd,P,T0):
	try:
		Ojd=numpy.array(Ojd)
		phase=[(x-float(T0))/float(P)%1. for x in Ojd]
	except:
		phase=(Ojd-float(T0))/float(P)%1.
	return phase


#fourier series
def fourier(x, *a):
	ret = a[0]/2 + a[1]* np.cos(2*np.pi*x) + a[2]* np.sin(2*np.pi*x)
	i=2
	for deg in range(3, len(a)-1,2):
		ret += a[deg]*np.cos(2*i*np.pi*x)+a[deg+1]*np.sin(2*i*np.pi*x)
		i+=1
	return ret

#error of fit in a point
def errorfourier(x, *s):
	ret=s[0]**2+(s[1]* np.cos(2*np.pi*x))**2 + (s[2]* np.sin(2*np.pi*x))**2
	i=2
	for deg in range(3, len(s)-1,2):
		ret+=(s[deg]*np.cos(2*i*np.pi*x))**2+(s[deg+1]*np.sin(2*i*np.pi*x))**2
		i+=1
	return np.sqrt(ret)

#integrated fourier series
def intfourier(x, *a):
	ret = a[1]*np.sin(2*np.pi*x)/(2*np.pi) - a[2]* np.cos(2*np.pi*x)/(2*np.pi)
	i=2
	for deg in range(3, len(a)-1,2):
		ret += a[deg]*np.sin(2*i*np.pi*x)/(2*i*np.pi) - a[deg+1]*np.cos(2*i*np.pi*x)/(2*i*np.pi)
		i+=1
	return ret


#Linear bisector fit (OLS bisector)
#slope + error & intercept coefficient + error
#the result is line that bisects lines of OLS fit (Y|X) and OLS fit (X|Y)
#based on work Isobe et al. (1990)

def lin_bisect(x, y):

	xm=np.mean(x)
	ym=np.mean(y)
	y0=y-ym
	x0=x-xm
	n = x.size
	
	sxx=np.sum(np.power(x-xm,2))
	syy=np.sum(np.power(y-ym,2))
	sxy=np.sum(np.multiply(x-xm,y-ym))
	

	b1=sxy/sxx
	b2=syy/sxy
	b3=np.power((b1+b2),-1)*(b1*b2-1+np.sqrt((1+b1**2)*(1+b2**2)))

	varb1=np.sum(np.multiply(np.power(x0,2),np.power(y0-b1*x0,2)))/sxx**2
	varb2=np.sum(np.multiply(np.power(y0,2),np.power(y0-b2*x0,2)))/sxy**2
	
	covb1b2=np.power(b1*sxx**2,-1)*np.sum(x0*y0*(y0-b1*x0)*(y0-b2*x0))

	varb3=((1+b2**2)**2*varb1+2*(1+b1**2)*(1+b2**2)*covb1b2+(1+b1**2)**2*varb2)*b3**2/((b1+b2)**2*(1+b1**2)*(1+b2**2))

	a1 = ym-b1*xm
	a2 = ym-b2*xm
	a3 = ym-b3*xm

	g1 = b3*((b1+b2)*np.sqrt((1+b1**2)*(1+b2**2)))**(-1)
	g13 = g1*(1+b2**2)
	g23 = g1*(1+b1**2)

	var_a3 = np.sum((y0-b3*x0-n*xm*((g13/sxx)*x0*(y0-b1*x0)+g23/sxy*y0*(y0-b2*x0)))**2)/n**2
	var_a1 = np.sum((y0-b1*x0-n*xm*((x0/sxx)*(y0-b1*x0)))**2)/n**2

	lista = []
	lista.append(b3)
	lista.append(a3)
	lista.append(varb3)
	lista.append(var_a3)

	return b3,a3,varb3,var_a3,var_a3

#Akima spline (smooth spline of the 3rd degree that does not oscillate between data points)
#Before that data is cleared using 3sigma rule (two times) on residuals of fit of a regular cubic spline
#Akima spline is interpolation on means calculated for 'moving binning', k-separation of bins and n-size of an interval from which data are taken to calculate means of x and y - are free parameters

def Akimaspline(x,y,k=0.05,n=0.05):
	xleft=[]
	for i in x:
		xleft.append(i-1)
	xright=[]
	for i in x:
		i=i+1
		xright.append(i)
	z=np.concatenate((xleft,x))
	x=np.concatenate((z,xright))
	z=np.concatenate((y,y))
	y=np.concatenate((y,z))
	
	#sortowanie po x
	indices=np.argsort(x)
	newy=[]
	for i in indices:
		newy.append(y[i])
	x.sort()
	y=np.asarray(newy)

	xnew = []
	ynew = []
	for i,iks in enumerate(x):
		
		if x[i] != x[i-1]:
			
			xnew.append(x[i])
			ynew.append(y[i])

	x = xnew
	y = ynew
	
	spline=UnivariateSpline(x,y,bbox=[-1,2],k=3,s=5)
	sigma=np.sqrt(spline.get_residual())

	i=0
	while i<len(x):
		if 3*sigma < np.fabs(spline(x[i])-y[i]):
			del x[i]
			del y[i]
		i+=1

	spline=UnivariateSpline(x,y,bbox=[-1,2],k=3,s=5)
	sigma=np.sqrt(spline.get_residual())
	#stad mozna miec blad dopasowania spline

	i=0
	while i<len(x):
		if 3*sigma < np.fabs(spline(x[i])-y[i]):
			del x[i]
			del y[i]
		i+=1
	

	#krok sredniej ruchomej:
	#k=0.02
	#przedzial fazy brany do sredniej:
	#n=0.1

	i=-1.
	xmean=[]
	ymean=[]

	while i<=2.0:
		j=0
		xbin=[]
		ybin=[]
		while j<len(x):
			if x[j]<i+n/2. and x[j]>i-n/2.:
				xbin.append(x[j])
				ybin.append(y[j])	
			j+=1
		if len(xbin)>0 and np.mean(xbin) not in xmean:
			xmean.append(np.mean(xbin))
			ymean.append(np.mean(ybin))
		i+=k
	akima=Akima1DInterpolator(xmean,ymean)

	return akima


#zaznaczanie fazy
def what_phase(data,m0,P):
	import ephem
	now = ephem.julian_date(data)
 # now=now-2450000
	faza1=(float(now)-float(m0))/float(P)%1
	faza2=(float(now+0.5)-float(m0))/float(P)%1
	return faza1,faza2

#liczenie odchylenia standardowego punktow od funkcji
def std_dev(f, x, y):
	np.vectorize(f)
	y_mean=f(x)
	std_d=np.sqrt(np.sum(np.power(np.subtract(y,y_mean),2)))
	return std_d

	
