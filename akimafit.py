#!/usr/bin/env python
import scipy
from bw_lib import Akimaspline
import os,sys
import matplotlib.pyplot as plt
import numpy as np
from random import sample
from abl_fit import param,random
from scipy.optimize import curve_fit
	

def toflux(x):
	return np.power(10.,np.array(x)/-2.5)

def tomag(x):
	return -2.5*np.log10(np.array(x))

def nottoflux(x):
	return x
def nottomag(x):
	return x



class signal():
	def __init__(self,plik,tcol,mcol,ercol,period,k,n,shift=0.,name='',sc=None,hjd0=None,flux=True,ph=True,tomean=False,ps=0.,downlimit=-500.,uplimit=500,template_file = ''):
		self.plik = plik
		self.tcol=tcol
		self.mcol=mcol
		self.ercol=ercol
		self.period = period
		self.shift = shift
		print(self.period)
		self.name = name
		self.k = k#knots step
		self.n = n#size of bin
		self.time = []#hjds
		self.mag = []#mags/rvs
		self.magerr = []#errors
		self.x = None#modeled lc phase
		self.y = None#modeled lc mags or rv
		self.yerr = None
		self.phasepts = []#array of objects of "parameter" class used to calculate mean value and error in each phase
		self.phasepts_plot = []
		self.phase = None#phase matrix
		self.sources = []#sources of photometry
		self.flux_orig = []
		self.mag_orig = []
		self.ps=ps
		self.tomean=tomean
		self.template_file = template_file.split('[')[0]
		try:
			self.template_tcol = int(template_file.split('[')[1].split(',')[0])
			self.template_mcol = int(template_file.split('[')[1].replace(']','',1).split(',')[1])
		except:
			self.template_tcol = 0
			self.template_tcol = 1

		
		self.read(tcol,mcol,ercol,sc=sc,downlimit=downlimit,uplimit=uplimit)
		
		if hjd0 != None:
			self.findhjd = False
			self.hjd0 = float(hjd0)
			
		else:
			self.findhjd = True
			self.hjd0 = self.time[np.argmin(self.mag)]
		if ph == True:
			self.phasing(self.period)
		else:
			self.phase=self.time
		
		if len(self.template_file) > 0:
			self.read_template()
			self.template_fit()

		
		
		
		
		if flux == True:
			self.toflux = toflux
			self.tomag = tomag
			self.fluxerr = np.power(10.,np.array(self.mag)/-2.5)*np.log(10.)*np.array(self.magerr)
		else:
			self.toflux = nottoflux
			self.tomag = nottomag
			self.fluxerr = self.magerr

		self.flux = self.toflux(self.mag)
		
		self.x,self.y,self.akima = self.fit(self.phase,self.flux,self.k,self.n)
		self.meanval = self.tomag(self.mean(self.akima))
		

	def read(self,tc,mc,ec,sc=None,downlimit=-500.,uplimit=500):
		plik = os.popen('cat '+self.plik).read().split('\n')[:-1]
		for linia in plik:
			if '#' not in linia:
				
				if ',' in linia:
					dane = linia.split(',')
				else:
					dane = linia.split()

				if float(dane[mc]) > downlimit and float(dane[mc]) < uplimit:
					self.time.append(float(dane[tc]))
					self.mag.append(float(dane[mc])+self.shift)
					self.magerr.append(float(dane[ec]))
					
					#MINIMAL PHOTOMETRIC ERROR IS SET TO 0.02!!!!!!!!!!!!!!!!
					#if self.magerr[-1] < 0.02:
					#	self.magerr[-1] = 0.02

					if sc != None:
						self.sources.append(dane[sc])
					else:
						self.sources.append('')

	def f_shift(self,x,a,b,c):#return template model but shifted to fit the data
		f= np.multiply(np.add(self.template_model(np.add(x,b)),a),c)
		return f
	
	def template_fit(self):
		p0 = (1.,np.mean(self.mag)-np.mean(self.template_mag),1.)
		popt,pcov=curve_fit(self.f_shift,self.phase,self.mag,p0=p0)
		self.templ_mag_shift = popt[0]
		self.templ_phase_shift = popt[1]
		self.templ_ampl = popt[2]
		if self.templ_phase_shift > 1.:
			self.templ_phase_shift = self.templ_phase_shift - 1.

		print(self.templ_mag_shift,self.templ_phase_shift,self.templ_ampl)

	def read_template(self):
		txt = os.popen('cat '+self.template_file).read().split('\n')[:-1]
		self.template_time = []
		self.template_mag = []
		for linia in txt:
			
			if '#' not in linia and ',' in linia:
				self.template_time.append(float(linia.split(',')[self.template_tcol]))	
				self.template_mag.append(float(linia.split(',')[self.template_mcol]))
			elif '#' not in linia:
				self.template_time.append(float(linia.split()[self.template_tcol]))	
				self.template_mag.append(float(linia.split()[self.template_mcol]))

		
		self.template_phasing(self.period)
		self.template_model = Akimaspline(self.template_phase,self.template_mag,k=0.04, n=0.04)
		sigma = np.std(np.subtract(self.template_mag,self.template_model(self.template_phase)))
		print("sigma",sigma)
		time = []
		mag = []
		for i,pt in enumerate(self.template_phase):
			if np.abs(np.subtract(self.template_mag[i],self.template_model(self.template_phase[i]))) < 2.*sigma:
				time.append(self.template_time[i])
				mag.append(self.template_mag[i])

		self.template_time = time
		self.template_mag = mag
		self.template_phasing(self.period)
		self.template_model = Akimaspline(self.template_phase,self.template_mag,k=0.04, n=0.05)
		plt.plot(self.template_phase,self.template_mag,'.',color='black')
		phases = np.arange(0,1,0.01)
		plt.plot(phases,self.template_model(phases),'-',color='green')
		plt.gca().invert_yaxis()
		plt.gca().tick_params(bottom=True,top=True,left=True,right=True,labelsize=14,direction='in')
		plt.xlabel('$Phase$',fontsize=16)
		plt.ylabel('$V_{ASAS}$',fontsize=16)
		plt.show()

	
	def template_phasing(self,period=None,ps=0.):
		self.template_phase = []
		if period==None:
			period=self.period

		print('period',period)
		try:
			for i,t in enumerate(self.template_time):
				periodt = period[0]+period[1]*(self.template_time[i]-period[2])/365.25
				self.template_phase.append((((self.template_time[i]-self.hjd0+self.ps*periodt)%periodt)/periodt))
				
		except:
			self.template_phase = np.divide(np.mod(np.add(np.subtract(self.template_time,self.hjd0),self.ps*period),period),period)

		self.template_phase = np.array(self.template_phase)
	 

	def phasing(self,period=None,ps=0.):
		self.phase = []
		if period==None:
			period=self.period

		#try:
		if True:
			
			for i,t in enumerate(self.time):
				direction = -1.
				if self.time[i] < self.hjd0:
					direction = 1.

				hjd0 = self.hjd0
				while direction*(hjd0 - self.time[i]) < 0.:
					print(t-period[2])
					periodt = period[0]+period[1]*(hjd0-period[2])/365.25
					hjd0 = hjd0 + direction*periodt
					print(hjd0)
				

				periodt = period[0]+period[1]*(self.time[i]-period[2])/365.25
				self.phase.append((((self.time[i]-hjd0+self.ps*periodt)%periodt)/periodt))
				
		#except:
		#	self.phase = np.divide(np.mod(np.add(np.subtract(self.time,self.hjd0),self.ps*period),period),period)

		self.phase = np.array(self.phase)

	def singlephase(self,time):
		return ((self.ps*self.period+(time-self.hjd0))%self.period/self.period)

	

	def fit(self,phase=None,flux=None,k=None,n=None,sort=True):
		if k != None:
			self.k = k
		if n != None:
			self.n = n
	
		if self.tomean:
			if len(self.mag_orig) > 1:
				self.mag = self.mag_orig
				self.flux=self.toflux(self.mag)

		try:

			if phase == None:

				phase = self.phase
				flux = self.flux

		except:
			pass
	
		sort=True
		if sort:
			sortedphase=np.argsort(phase)
			p = []
			f = []
			for pt in sortedphase:
				#if phase[pt] not in p:
				p.append(phase[pt])
				f.append(flux[pt])

		else:
			p = phase
			f = flux
		
		
		
		akima = Akimaspline(p,f,k=self.k,n=self.n)
		if self.tomean:
			f = f-self.mean(akima)
			try:
				if self.mag_orig == None:
					self.mag_orig = self.mag
					self.mag = self.mag-self.mean(akima)
					self.flux_orig = self.toflux(self.mag_orig)
					self.flux = self.toflux(self.mag)
				else:
					self.mag = self.mag-self.mean(akima)
					self.flux = self.toflux(self.mag)
					
			except:
				pass

		akima = Akimaspline(p,f,k=self.k,n=self.n)
		self.akima = akima
		#print self.akima
		self.x = np.arange(0,1.001,0.001)
		#print self.x
		self.y = self.tomag(akima.__call__(self.x))
		#print self.y
		return self.x,self.y,akima

	def mean(self,akima):
		return self.integral(akima,0,1)-self.shift
	
	def integral(self,akima,p0,p1):
		integral = akima.antiderivative()
		return integral(float(p1))-integral(float(p0))

	def rmpt(self,pt):
		self.phase_del.append(self.phase[pt])
		self.mag_del.append(self.mag[pt])
		phasenew = []
		magnew = []
		for i,p in enumerate(self.phase):
			if i != pt:
				phasenew.append(self.phase[i])
				magnew.append(self.mag[i])

		self.phase = phasenew
		self.mag = magnew

	def addpt(self,pt):
		#phase
		#for i,pt in enumerate(self.phase_del):
			#if i !
		self.phase_del.append(self.phase[pt])
		self.mag_del.append(self.mag[pt])
		phasenew = []
		magnew = []
		for i,p in enumerate(self.phase):
			if i != pt:
				phasenew.append(self.phase[i])
				magnew.append(self.mag[i])

		self.phase = phasenew
		self.mag = magnew


	def clean(self,nsigma):
		
		newtime = []
		newphase = []
		newmag = []
		newmagerr = []
		newflux = []
		newfluxerr = []
		newsources = []
		for i,pt in enumerate(self.time):
			#print np.abs(self.akima.__call__(self.phase[i]) - self.flux[i])/self.flux[i]
			if np.abs(self.akima.__call__(self.phase[i]) - self.flux[i])/self.flux[i] <= 2.:
				newtime.append(self.time[i])
				newphase.append(self.phase[i])
				newmag.append(self.mag[i])
				newmagerr.append(self.magerr[i])
				newflux.append(self.flux[i])
				newfluxerr.append(self.fluxerr[i])
				if len(self.sources) > 0:
					newsources.append(self.sources[i])

		self.time = np.array(newtime)
		self.phase = np.array(newphase)
		self.mag = np.array(newmag)
		self.magerr = np.array(newmagerr)
		self.flux = np.array(newflux)
		self.fluxerr = np.array(newfluxerr)
		if len(self.sources) > 1:
			self.sources = newsources
	

	def mc(self,n,phases,fraction = 1.0):
		self.mean_mc = param('mean')
		#testpt = param('test')
		for i,pt in enumerate(self.x):
			self.phasepts_plot.append(param(str(pt)))

		for i,pt in enumerate(phases):
			self.phasepts.append(param(str(pt)))

		for i in range(n):
			
			#subset = sample(np.arange(0,len(self.flux),1),int(fraction*len(self.flux)))
			subset = np.arange(0,len(self.flux),1)
			fluxnew_orig=[]
			fluxerr_new = []
			phasenew = []
			for el in subset:
				fluxnew_orig.append(self.flux[el])
				fluxerr_new.append(self.fluxerr[el])
				phasenew.append(self.phase[el])

			fluxnew = []
			
			for i,pt in enumerate(subset):
				
				fluxnew.append(random(fluxnew_orig[i],fluxerr_new[i]))
				
			
			x,y,akima = self.fit(phasenew,fluxnew,self.k,self.n)
			self.mean_mc.values.append(self.tomag(self.mean(akima)))
			for i,pt in enumerate(self.x):
				self.phasepts_plot[i].values.append(y[i])
			
			for i,pt in enumerate(phases):
				self.phasepts[i].values.append(akima.__call__(pt))
			

			
		self.y = []
		self.yerr = []	
		for i,pt in enumerate(self.x):
			self.phasepts_plot[i].srednia()
			self.y.append(self.phasepts_plot[i].mean)
			self.yerr.append(self.phasepts_plot[i].std)
			self.mean_mc.srednia()

			#if self.findhjd:
			#	self.find_hjd0()
			#	pl = open('hjd0.dat','a')
			#	pl.write(self.name+'\t'+str(self.hjd0)+'\n')
			#	pl.close()

		self.fluxtocalc=[]
		self.fluxerrtocalc=[]
		for i,pt in enumerate(phases):
			self.phasepts[i].srednia()
			self.fluxtocalc.append(self.phasepts[i].mean)
			self.fluxerrtocalc.append(self.phasepts[i].std)
			
		#print len
	
	def find_hjd0(self):
		
		maxt = np.max(self.time)
		ph = self.singlephase(maxt)
		#print maxt,ph
		c = 0
		for i,x in enumerate(self.x):
			if self.line[i] >= self.mean_mc.mean and c == 0:
				continue
			elif self.line[i] < self.mean_mc.mean and c == 0:
				c=1
				continue
			elif self.line[i] < self.mean_mc.mean and c == 1:
				c=1
				continue
			elif self.line[i] >= self.mean_mc.mean and c == 1:
				phase0 = x
				break
		#print phase0
		#print self.hjd0
		self.hjd0 = maxt+(phase0-ph)*self.period
		#print self.hjd0	
		self.x= np.mod(np.subtract(self.x,phase0),1.)
		self.phase = np.mod(np.subtract(self.phase,phase0),1.)
		ind = np.argsort(self.x)
		#print ind,ind.tolist()
		#self.x = np.array([self.x[i] for i in ind.tolist()])
		#self.line = np.array([self.line[i] for i in ind.tolist()])

		
		
def main(lst):
	timecol = 1
	magcol = 7
	errcol = 8
	scols = None
	
	if '-h' in lst:
		print('Usage: akimafit.py file_with_data[timecol,magcol,errcol] object_name -p <float> -k <int> -n <int> -mc <int> -f <float> -s <str> -o <str> -hjd0 <float> -flux\n\n')
		print('\t-p - period\n')
		print('\t-k - step between two knots (default 0.02)\n')
		print('\t-n - width of the range to calculate the mean value for given knot (default 0.05)\n')
		print('\t-mc - number of Monte Carlo simulations (default 2000)\n')
		print('\t-f - fraction of points used in a single MC simulation (default 1.0)\n')
		print('\t-s - name of the source of data\n')
		print('\t-o - output file name for mean value\n')
		print('\t-hjd0 - heliocentric julian date used for the zero phase\n')
		print('\t-flux - if used, magnitudes are converted into flux before fitting\n')
		print('\n\npwielgor@camk.edu.pl, 06.2021')
		sys.exit()


	if '-p' in lst:
			period = float(lst[lst.index('-p')+1])
	else:
		try:
			datafile = os.popen('cat /home/piotr/Dokumenty/science/MW_targets/analysis/GaiaEDR3/paper/MWT2CEP_data.csv').read().split('\n')[1:-1]
			names = []
			periods = []
			for linia in datafile:
			
				if '#' not in linia:
					names.append(linia.split(',')[0])
					periods.append(float(linia.split(',')[3]))

			if ob in names:
				ind = names.index(ob)
				period = periods[ind]

			else:
				print('Object '+ob+' not found in /home/piotr/Dokumenty/science/MW_targets/analysis/GaiaEDR3/paper/MWT2CEP_data.csv')
				sys.exit()
		except:
		
			print('Specify period with -p option and run again\n')
			sys.exit()
		

	plik = lst[1]
	if '[' in plik:
		cols = plik.split('[')[1].replace(']','',1).split(',')
		timecol = int(cols[0])
		magcol = int(cols[1])
		errcol = int(cols[2])
		try:
			scols = int(cols[3])
		except:
			scols = None
		plik = plik.split('[')[0]
		
	ob = lst[2]
	if '-k' in lst:
		ind = lst.index('-k')
		k = float(lst[ind+1])

	else:
		k=0.02

	if '-n' in lst:
		ind = lst.index('-n')
		n = float(lst[ind+1])
	else:
		n=0.05

	if '-mc' in lst:
		ind = lst.index('-mc')
		mc = int(lst[ind+1])
	else:
		mc = 2000

	if '-f' in lst:
		ind = lst.index('-f')
		frac = float(lst[ind+1])
	else:
		frac = 1.

	
	if '-s' in lst:
		ind = lst.index('-s')
		source = lst[ind+1]
	else:
		source = ''

	if '-o' in lst:
		ind = lst.index('-o')
		outfile = open(lst[ind+1],'a')
	else:
		outfile = None
	
	if '-2mass' in lst:
		ind = lst.index('-2mass')
		try:
			hjd = float(lst[ind+1])
			mag = float(lst[ind+2])
		except:
			hjd = None
			mag = None	
	else:
		hjd=None
		mag = None

	if '-hjd0' in lst:
		try:
			ind = lst.index('-hjd0')
			hjd0 = float(lst[ind+1])
		except:
			hjd0 = None
	else:
		hjd0 = None

	if '-flux' in lst:
		flux = True
	else:
		flux = False

	print('Akima spline is fitted with k='+str(k)+', n='+str(n)+'\n'+str(mc)+' Monte Carlo iterations with random subsample of '+str(frac)+ ' of full dataset\n')
	
	st = signal(plik,timecol,magcol,errcol,period,k,n,name=ob,sc=scols,hjd0=hjd0,flux=flux)
	
	if '-ax' in lst:
		ind = lst.index('-ax')
		'''lst[ind+1].plot(st.x,st.y,color='red')
		lst[ind+1].plot(st.x+1,st.y,color='red')
		lst[ind+1].plot(st.phase,st.mag,'.',color='black')
		lst[ind+1].plot(st.phase+1,st.mag,'.',color='black')
		lst[ind+1].hlines(st.meanval,0,2,label=str('{:.3f}'.format(st.meanval)))
		lst[ind+1].gca().invert_yaxis()
		lst[ind+1].legend()
		lst[ind+1].show()'''
		st.mc(mc,fraction = frac,ax = lst[ind+1],source=source,mass2=[hjd,mag])
		
			
	
	else:
		plt.plot(st.x,st.y,color='red')
		plt.plot(st.x+1,st.y,color='red')
		plt.plot(st.phase,st.mag,'.',color='black')
		plt.plot(st.phase+1,st.mag,'.',color='black')
		plt.hlines(st.meanval,0,2,label=str('{:.3f}'.format(st.meanval)))
		plt.gca().invert_yaxis()
		plt.legend()
		plt.show()
		st.mc(mc,[],fraction = frac)#ax=None)
	
	if outfile != None:
		outfile.write(ob+','+str(st.mean_mc.mean)+','+str(st.mean_mc.std)+'\n')
		outfile.close()
	else:
		print(ob+','+str(st.mean_mc.mean)+','+str(st.mean_mc.std)+'\n')

if __name__ == '__main__':
	main(sys.argv)
	
