#!/usr/bin/env python
from scipy.optimize import curve_fit
from scipy.stats import norm
import os,sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib
matplotlib.use('TKAgg')

def abl_right(x,a,b):
	ab = np.power(10.,0.2*(a*(np.log10(x))+b))

	return ab

def abl_left(parallax, mag):
	ab = parallax*np.power(10.,np.subtract(np.multiply(0.2,mag),2.))
	return ab

def absmag(parallax,mag):
	am = mag - 5.*np.log10(1000./parallax)+5.
	return am


def random(mean,sigma):
	return np.random.normal(loc=mean,scale=sigma)

def calka(x,p):
	calka = 0.
	if len(x) == len(p):
		for i,xi in enumerate(x):
			if i > 0:
				calka = calka + (xi-x[i-1])*0.5*(p[i]+p[i-1])
			
	elif len(x) == len(p)+1:
		for i,xi in enumerate(x):
			if i > 0:
				calka = calka + (xi-x[i-1])*p[i-1]
	#print(calka)
	return calka

class param():
	def __init__(self,name):
		self.name = name
		self.values = []
		
	def srednia(self):
		#self.mean = np.mean(self.values)
		#self.std = np.std(self.values)
		self.mean, self.std = norm.fit(self.values)
		#print(self.name+'='+str(self.mean)+' +- '+str(self.std))

	def hist(self,nbins):
		
		n, bins, patches = plt.hist(self.values,bins=nbins)#,normed=1)
		chist = calka(bins,n)
		
		xmin, xmax = plt.xlim()
		x = np.linspace(xmin, xmax, nbins)
		p = norm.pdf(x, self.mean, self.std)#*np.power(2.*np.pi,0.5)
		c = calka(x,p)
		#plt.vlines(self.mean,0,np.max(p*chist),color = 'orange',alpha=0.8)
		#plt.plot(x, chist*p, 'r--', linewidth=2,alpha=0.7,color='orange')
		plt.xlabel(self.name,fontsize=16)
		plt.ylabel('$N$',fontsize=16)
		plt.gca().tick_params(bottom=True,top=True,left=True,right=True,labelsize=14,direction='in')

		y= norm.pdf(bins,self.mean,self.std)*chist/c
		l = plt.plot(bins,y,'r--',linewidth=2)
		plt.show()
		
		
class stars():
	def __init__(self,plik,namecol=0,plxcol=1,plxerrcol=2,percol=3,magcol=4,magerrcol=5,redcol=6,rederrcol=7,redcoef = 0.):
		self.names = []
		self.plxs = []
		self.plxerrs = []
		self.periods = []
		self.mags = []
		self.magerrs = []
		self.read(plik,namecol,plxcol,plxerrcol,percol,magcol,magerrcol,redcol,rederrcol,redcoef)

		self.plxs_random = []
		self.mags_random = []
		

	def read(self,plik,namecol,plxcol,plxerrcol,percol,magcol,magerrcol,redcol,rederrcol,redcoef):
		if 'posix' in os.name:
			dane = os.popen('cat '+plik).read().split('\n')[:-1]
		else:
			dane=open(plik,'r')
		for linia in dane:
			if '#' not in linia:
				if ',' in linia:
					s = linia.split(',')
				else:
					s = linia.split()

				self.names.append(s[namecol])
				self.plxs.append(float(s[plxcol]))
				self.plxerrs.append(float(s[plxerrcol]))
				self.periods.append(float(s[percol]))
				self.mags.append(float(s[magcol])-redcoef*float(s[redcol]))
				self.magerrs.append(np.power(np.power(float(s[magerrcol]),2.)+np.power(float(s[rederrcol]),2.),0.5))

		magerrs_plx_brighter = (np.log10(np.array(abl_left(np.add(self.plxs,self.plxerrs),self.mags))) - np.log10( np.array(abl_left(self.plxs,self.mags))) )* 5.
		magerrs_plx_fainter = (np.log10(np.array(abl_left(np.subtract(self.plxs,self.plxerrs),self.mags))) - np.log10( np.array(abl_left(self.plxs,self.mags))) )* 5.
		self.magerrs_total = [np.power(np.power(magerrs_plx_brighter,2.)+np.power(self.magerrs,2.),0.5),np.power(np.power(magerrs_plx_fainter,2.)+np.power(self.magerrs,2.),0.5)]
				
	def mc(self):
		self.plxs_random = []
		self.mags_random = []
		for i,p in enumerate(self.plxs):
			self.plxs_random.append(random(self.plxs[i],self.plxerrs[i]))
			self.mags_random.append(random(self.mags[i],self.magerrs[i]))

		#plt.plot(self.plxs,np.array(self.plxs) - np.array(self.plxs_random),'.')
		#plt.plot(self.mags,np.array(self.mags) - np.array(self.mags_random),'.')
		#plt.plot(np.log10(self.periods),absmag(np.array(self.plxs_random),np.array(self.mags_random)),'.')
		#plt.plot(np.log10(self.periods),absmag(np.array(self.plxs),np.array(self.mags)),'.')
		#plt.show()

	def fit_abl(self):
		popt,pcov = curve_fit(abl_right,self.periods,abl_left(self.plxs_random,self.mags_random))
		self.a = popt[0]
		self.b = popt[1]
		#print popt,pcov
		return popt 
		
		
def main(lst):

	try:
		filtr = 'mag'
		if 'J' in lst:
			magcol = 10
			magerrcol = 11
			redcoef = 0.892
			filtr = 'J'
		if 'H' in lst:
			magcol = 12
			magerrcol = 13
			redcoef = 0.553
			filtr = 'H'
		if 'K' in lst:
			magcol = 14
			magerrcol = 15
			redcoef = 0.363
			filtr = 'K'

		if '-n' in lst:
			ind = lst.index('-n')
			n = int(lst[ind+1])
		else:
			n = 10000

		print('Number of simulations '+str(n))

		if '-b' in lst:
			ind = lst.index('-b')
			bins = int(lst[ind+1])
		else:
			bins = 100
		print('Number of histogram bins '+str(bins))
		
		gwiazdy = stars(lst[1],namecol=0,plxcol=4,plxerrcol=5,percol=3,magcol=magcol,magerrcol=magerrcol,redcol=16,rederrcol=17,redcoef = redcoef)
		a = param('a')
		b = param ('b')
		for i in range(0,n,1):
			gwiazdy.mc()
			gwiazdy.fit_abl()
			a.values.append(gwiazdy.a)
			b.values.append(gwiazdy.b)

		if True:
			plt.errorbar(np.log10(gwiazdy.periods),absmag(np.array(gwiazdy.plxs),np.array(gwiazdy.mags)),yerr=gwiazdy.magerrs_total,fmt='.')
			x = np.arange(0,1.5,0.1)
			
			a.srednia()
			b.srednia()
			plt.plot(x,a.mean*x+b.mean)
			plt.xlabel('log(P)')
			plt.ylabel(filtr)
			if '-names' in lst:
				sort = np.argsort(gwiazdy.periods)
				for i,gw in enumerate(sort):
					plt.text(np.log10(gwiazdy.periods[gw])-0.01,absmag(np.array(gwiazdy.plxs[gw]),np.array(gwiazdy.mags[gw]))-0.3*np.power(-1,i),gwiazdy.names[gw],rotation='vertical',fontsize=8,color='black',alpha=0.8)
			plt.gca().invert_yaxis()
			plt.show()
			a.hist(bins)
			b.hist(bins)
			
	except:
		print('Usage: abl_fit.py datafile filter -n <num_of_simulations> -b <number_of_bins>\n\npwielgor@camk.edu.pl, 04.2021\n')
		sys.exit()

if __name__ == '__main__':
	main(sys.argv)
