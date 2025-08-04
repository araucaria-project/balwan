#!/usr/bin/env python
import os,sys
import numpy as np
import scipy
from pwlib import metnk,metnk_uncer
from bw_lib import *
from abl_fit import param,random
from bw_lib import Akimaspline

def bisect_my(x,xerr,y,yerr):
	a1,b1,sigmaa1,sigmab1,sigma1 = metnk_uncer(x,y,yerr)
	a2,b2,sigmaa2,sigmab2,sigma2 = metnk_uncer(y,x,xerr)
	return (a1+1./a2)/2., b1,sigmaa1,sigmab1,sigma1

class balwan():
	def __init__(self,name,v,k,rv,ebv,ebverr=0,R_V=3.136,R_K=0.363,plx=None,plxerr=None,irsb=[3.953,-0.1336],dirsb=0,mc=0,p0=0.,p1=1.,method = 1):
		self.mean_rad = 0.
		self.name = name
		self.kmkpc=206264.81*149597870.7*1000.
		self.rv = rv
		self.v = v
		self.k = k
		self.ebv = ebv
		self.ebverr = ebverr
		self.plx = plx
		self.plxerr = plxerr
		self.irsbterms=irsb
		self.dirsb = dirsb 
		self.method = method
		
		self.p0 = p0
		self.p1 = p1
		self.phase = []
		self.phase_del = []
		for i,pt in enumerate(self.k.phase):
			if pt >= self.p0 and pt < self.p1:
				self.phase.append(pt)
			else:
				self.phase_del.append(pt) 
		
	
		self.R_V=R_V
		self.R_K=R_K

		self.phase_plot = np.arange(0,1.001,0.001)
		self.v_kphase = []
		self.del_v_kphase = []
		self.vk = []
		self.dr = []	
		self.dr_del = []
		self.fi = []
		self.fi_del = []
		self.vk_err = []
		self.dr_err = []
		self.dr_err_del = []
		self.fi_err = []	
		self.fi_err_del = []
		self.v_plot = []
		self.k_plot = []
		self.rv_plot = []
		self.vk_plot = []

		self.dr_plot = []
		self.fi_plot = []

		self.runbw(mc=mc)

	def integral(self,akima,p0,p1):
		integral = akima.antiderivative()
		return integral(float(p1))-integral(float(p0))

	def intrv(self):
		self.dr_plot = []
		self.dr = []
		self.dr_del = []
		for i,pt in enumerate(self.phase_plot):
			self.dr_plot.append(-1.*self.rv.integral(self.rv.akima,0,pt)*self.v.period[0]*86400.)
			#self.dr_err.plot.append()
			
		for i,pt in enumerate(self.k.phase):
			if pt >= self.p0 and pt < self.p1:
				self.dr.append(-1.*self.rv.integral(self.rv.akima,0,pt)*self.v.period[0]*86400.)
			else:
				self.dr_del.append(-1.*self.rv.integral(self.rv.akima,0,pt)*self.v.period[0]*86400.)
			#self.dr_err.append

		self.mean_rad = 0.
		for i,pt in enumerate(self.phase_plot):
			if i > 0:
				self.mean_rad = self.mean_rad + (self.phase_plot[i]-self.phase_plot[i-1])*(self.dr_plot[i]+self.dr_plot[i-1])/2.

		#print self.phase_plot.tolist()
		
		akima = Akimaspline(self.phase_plot[:-1],self.dr_plot[:-1],k=0.02,n=0.02)
		self.dr_akima = akima
		#self.mean_rad = self.integral(akima,0,1)
		self.dr = np.array(self.dr)-self.mean_rad
		self.dr_del = np.array(self.dr_del)-self.mean_rad
		self.dr_plot = np.array(self.dr_plot)-self.mean_rad
		

	def calc_vk(self):
		self.v_plot = []
		self.vk_plot = []
		self.v_kphase= []
		self.vk = []
		if self.v.template_file == '':
			for i,pt in enumerate(self.phase_plot):
				self.v_plot.append(self.v.tomag(self.v.akima.__call__(pt)))
				self.vk_plot.append(self.v.tomag(self.v.akima.__call__(pt))-self.k.tomag(self.k.akima.__call__(pt)))
			for i,pt in enumerate(self.k.phase):
				self.v_kphase.append(self.v.tomag(self.v.akima.__call__(pt)))
				self.vk.append(self.v.tomag(self.v.akima.__call__(pt))-self.k.mag[i])
				#self.vk_err.append'''

	
		else:
			for i,pt in enumerate(self.phase_plot):
				
				self.v_plot.append(self.v.f_shift(pt,self.v.templ_mag_shift,self.v.templ_phase_shift,self.v.templ_ampl))
				self.vk_plot.append(self.v.f_shift(pt,self.v.templ_mag_shift,self.v.templ_phase_shift,self.v.templ_ampl)-self.k.tomag(self.k.akima.__call__(pt)))
			for i,pt in enumerate(self.k.phase):
				self.v_kphase.append(self.v.f_shift(pt,self.v.templ_mag_shift,self.v.templ_phase_shift,self.v.templ_ampl))
				self.vk.append(self.v.f_shift(pt,self.v.templ_mag_shift,self.v.templ_phase_shift,self.v.templ_ampl)-self.k.mag[i])
			

	def irsb(self,vk,ebv):
		vk = np.array(vk)
		irsb = np.zeros(len(vk))
		
		for j,pt in enumerate(irsb):
			for i,el in enumerate(self.irsbterms):
				irsb[j] = irsb[j]+el*np.power(vk[j]-(self.R_V-self.R_K)*ebv,float(i))
				if j == 0:
					irsb[j] = irsb[j]+self.dirsb
		return irsb
		#return np.add(np.multiply(self.irsba,np.subtract(vk,(self.R_V-self.R_K)*ebv)),self.irsbb)

	def calc_fi(self,ebv):
		self.fi = []
		self.fi_del = []
		self.fi_plot = []
		self.fi_plot2 = []
		fv = self.irsb(self.vk,self.ebv)
		#this is for radians 
		fi = np.multiply(np.power(10.,np.multiply(-2.,np.add(np.multiply(0.1,np.subtract(self.v_kphase,self.R_V*ebv)),np.subtract(fv,4.2207)))),np.pi/(180.*3600.*1000.))
		for i,pt in enumerate(self.k.phase):
			if pt >= self.p0 and pt < self.p1:
				self.fi.append(fi[i])
			else:
				self.fi_del.append(fi[i]) 
		#this is for miliarcseconds
		#self.fi = np.power(10.,np.multiply(-2.,np.add(np.multiply(0.1,np.subtract(self.v_kphase,self.R_V*ebv)),np.subtract(fv,4.2207))))
		#fv2 = np.add(np.multiply(self.irsba,np.subtract(self.vk_plot,(self.R_V-self.R_K)*self.ebv)),self.irsbb)
		#self.fi_plot = np.multiply(np.power(10.,np.multiply(-2.,np.add(np.multiply(0.1,np.subtract(self.v_plot,self.R_V*self.ebv)),np.subtract(fv2,4.2207)))),np.pi/(180.*3600.*1000.))
		akima = Akimaspline(self.k.phase,fi,0.1,0.1)
		self.fi_akima = akima
		for i,pt in enumerate(self.phase_plot):
			self.fi_plot.append(akima.__call__(pt))
			
		self.fi = np.array(self.fi)
		self.fi_del = np.array(self.fi_del)
		self.fi_plot=np.array(self.fi_plot)

	
		
	def fit_dr2fi(self):
		pfs = np.arange(1.,2.3,0.01)
		
		fi0s = np.arange(np.min(self.fi),np.max(self.fi),(np.max(self.fi)-np.min(self.fi))/500.)
		#print len(fi0s)
		chi2_min = 10000000000.
		pf_final = None
		fi0_final = None
		for pf in pfs:
			for fi0 in fi0s:
				chi2 = 0.
				for i,pt in enumerate(self.k.phase):
					chi2 = chi2 + np.power(self.fi[i] - fi0 -(self.dr[i])*pf*2.*self.plx/self.kmkpc,2.)

				chi2 = np.power(chi2/len(self.k.phase),0.5)

				if chi2 <= chi2_min:
					pf_final = pf
					fi0_final = fi0
					chi2_min = chi2
					#print chi2,pf,fi0
		#print np.min(self.fi),np.max(self.fi),fi0_final,np.min(fi0s),np.max(fi0s)
		if pf_final == None:
			pf_final = 0
			fi0_final = 0

		pf_final_err = chi2_min/np.power(len(self.k.phase),0.5)


		return pf_final,fi0_final,pf_final_err,pf_final_err,chi2_min

	def fit_dr2fi_dst(self):
		dsts = np.arange(1.,2.0,0.01)
		
		fi0s = np.arange(np.min(self.fi),np.max(self.fi),(np.max(self.fi)-np.min(self.fi))/500.)
		#print len(fi0s)
		chi2_min = 10000000000.
		pf_final = None
		fi0_final = None
		for pf in pfs:
			for fi0 in fi0s:
				chi2 = 0.
				for i,pt in enumerate(self.k.phase):
					chi2 = chi2 + np.power(self.fi[i] - fi0 -(self.dr[i])*pf*2.*self.plx/self.kmkpc,2.)

				chi2 = chi2/len(self.k.phase)

				if chi2 <= chi2_min:
					pf_final = pf
					fi0_final = fi0
					chi2_min = chi2
					#print chi2,pf,fi0

		pf_final_err = chi2_min/np.power(len(self.k.phase),0.5)
		#print np.min(self.fi),np.max(self.fi),fi0_final,np.min(fi0s),np.max(fi0s)
		if pf_final == None:
			pf_final = 0
			fi0_final = 0


		return pf_final,fi0_final,pf_final_err,pf_final_err,chi2_min


	def bawan(self,plx=None,fi = [],dr = [],dr_err=[],fi_err=[]):
		
		#print kmkpc
		if plx != None:
			if self.method == 0:
				p,fi0,sigmap,sigmafi0,sigma = self.fit_dr2fi()
			elif self.method == 1:
				p,fi0,sigmap,sigmafi0,sigma=metnk(np.multiply(dr,2.*plx/self.kmkpc),fi)
			else:
				p,fi0,sigmap,sigmafi0,sigma=lin_bisect(np.multiply(dr,2.*plx/self.kmkpc),fi)

			#if len(dr_err) > 0:
			#	p,fi0,sigmap,sigmafi0,sigma=bisect_my(np.multiply(dr,2.*plx/self.kmkpc),np.multiply(dr_err,2.*plx/self.kmkpc),fi,fi_err)
			#else:
			#	p,fi0,sigmap,sigmafi0,sigma=metnk(np.multiply(dr,2.*plx/self.kmkpc),fi)
			sigma = 0.
			for i,pt in enumerate(fi):
				sigma = sigma+np.power(fi[i]-p*np.multiply(dr[i],2.*plx/self.kmkpc)-fi0,2.)

			sigma = np.power(sigma/float(len(fi)),0.5)
			#print p,fi0
			#print plx
			#print 'p=',0.5*(p*self.kmkpc/plx),'r0=',fi0*0.5*(self.kmkpc/plx)
			#print 'p=',p,'r0=',fi0*0.5*(self.kmkpc/plx)
		'''elif p != None:#p is the distance here
			if self.method == 0:
				plx,fi0,sigmaplx,sigmafi0,sigma = self.fit_dr2fi_dst()
			elif self.method == 1:
				plx,fi0,sigmaplx,sigmafi0,sigma=metnk(np.multiply(dr,2.*p/self.kmkpc),fi)
			else:
				plx,fi0,sigmaplx,sigmafi0,sigma=lin_bisect(np.multiply(dr,2.*p/self.kmkpc),fi)

			p = plx
			sigmap = sigmaplx'''
		print(p)


		return p,fi0,sigmap,sigmafi0,sigma

	def runbw(self,mc=0):
		p = param('$p-factor$')
		fi = param('$\Theta$')
		r = param('$R[R_{\odot}]$')
		vk_err_param = []
		fi_err_param = []
		dr_err_param = []
		fi_err_param_del = []
		dr_err_param_del = []
		self.vk_err = []
		self.fi_err = []
		self.dr_err = []
		self.fi_err_del = []
		self.dr_err_del = []
		for i,pt in enumerate(self.phase):
			vk_err_param.append(param(str(i)))
			fi_err_param.append(param(str(i)))
			dr_err_param.append(param(str(i)))

		for i,pt in enumerate(self.phase_del):
			fi_err_param_del.append(param(str(i)))
			dr_err_param_del.append(param(str(i)))

		if mc > 0:
			k_mag_orig = []
			v_mag_orig = []
			rv_mag_orig = []
			for i,pt in enumerate(self.k.phase):
				k_mag_orig.append(self.k.mag[i])
			
			for i,pt in enumerate(self.v.phase):
				v_mag_orig.append(self.v.mag[i])
			for i,pt in enumerate(self.rv.phase):			
				rv_mag_orig.append(self.rv.mag[i])
			
			for j in range(mc):
				#print 'tutu'

				'''self.vk = []
				self.fi = []
				self.dr = []'''
				
				
				if self.plxerr > 0:
					plx_mc = random(self.plx,self.plxerr)
				else:
					plx_mc = self.plx

				if self.ebverr > 0:
					ebv_mc = random(self.ebv,self.ebverr)
				else:
					ebv_mc = self.ebv

				
				for i,pt in enumerate(self.k.phase):
					self.k.mag[i]=random(k_mag_orig[i],self.k.magerr[i])

				for i,pt in enumerate(self.v.phase):
					self.v.mag[i]=random(v_mag_orig[i],self.v.magerr[i])
		
				for i,pt in enumerate(self.rv.phase):
					self.rv.mag[i] = random(rv_mag_orig[i],self.rv.magerr[i])	
					

				self.v.flux=self.v.toflux(self.v.mag)
				self.v.fit()
				
				self.rv.mag_orig = self.rv.mag
				self.rv.fit()
				self.k.flux=self.k.toflux(self.k.mag)
				self.k.fit()
				self.intrv()
				self.calc_vk()
				self.calc_fi(ebv=ebv_mc)
				for i,pt in enumerate(self.phase):
					vk_err_param[i].values.append(self.vk[i])
					fi_err_param[i].values.append(self.fi[i])
					dr_err_param[i].values.append(self.dr[i])


				for i,pt in enumerate(self.phase_del):
					fi_err_param_del[i].values.append(self.fi_del[i])
					dr_err_param_del[i].values.append(self.dr_del[i])

				self.p,self.fi0,self.sigmap,self.sigmafi0,self.sigma=self.bawan(plx_mc,self.fi,self.dr)
				if self.p < 2.:
					p.values.append(self.p)
					fi.values.append(self.fi0)
					r.values.append(self.fi0*0.5*(self.kmkpc/self.plx)/695700)
					out = open('mc_out.dat','a')
					out.write(str(self.p)+'\t'+str(self.fi0*0.5*(self.kmkpc/self.plx))+'\n')


			for i,pt in enumerate(self.k.phase):
				self.k.mag[i]=k_mag_orig[i]
				vk_err_param[i].srednia()
				self.vk_err.append(vk_err_param[i].std)

			for i,pt in enumerate(self.phase):
				fi_err_param[i].srednia()
				self.fi_err.append(fi_err_param[i].std)
				dr_err_param[i].srednia()
				self.dr_err.append(dr_err_param[i].std)
				#if i ==0:
				#	dr_err_param[i].hist(10)

			for i,pt in enumerate(self.phase_del):
				fi_err_param_del[i].srednia()
				self.fi_err_del.append(fi_err_param_del[i].std)
				dr_err_param_del[i].srednia()
				self.dr_err_del.append(dr_err_param_del[i].std)
			
			for i,pt in enumerate(self.v.phase):
				self.v.mag[i]=v_mag_orig[i]

			for i,pt in enumerate(self.rv.phase):
				self.rv.mag[i]=rv_mag_orig[i]

			self.x_err = np.multiply(self.dr_err,2.*self.plx/self.kmkpc)

			#print self.k.mag

			self.v.flux=self.v.toflux(self.v.mag)
			self.v.fit()
				
			self.rv.fit()
			self.k.flux=self.k.toflux(self.k.mag)
			self.k.fit()
			self.intrv()
			self.calc_vk()
			self.calc_fi(self.ebv)

			p.srednia()
			fi.srednia()
			r.srednia()
			p.hist(50)
			
			#fi.hist(50)
			self.p = p.mean
			self.sigmap = p.std
			self.fi0 = fi.mean
			self.sigmafi0 = fi.std
			
			#p.hist(15)
			#print('Not implemented yet')

			self.fi=np.array(self.fi)
			self.fi_err=np.array(self.fi_err)
			self.fi_del=np.array(self.fi_del)
			self.fi_err_del=np.array(self.fi_err_del)
		else:
			self.intrv()
			self.calc_vk()
			self.calc_fi(ebv=self.ebv)
			
			p,self.fi0,sigmap,self.sigmafi0,self.sigma=self.bawan(self.plx,self.fi,self.dr)
			
			self.p = float(p)
			self.sigmap = float(sigmap)
			#self.p,self.fi0
		print(self.p)
		for i,pt in enumerate(self.phase_plot):
			self.fi_plot2.append(np.multiply(self.dr_plot[i],2.*self.plx*self.p/(self.kmkpc))+self.fi0)

		self.fi_plot2=np.array(self.fi_plot2)	
			
	def text(self):
		#try:
		if True:
			t='object: '+self.name+'\n\nprojection factor ='+str("{:.3f}".format(self.p))+'+-'+str("{:.3f}".format(self.sigmap))+'\n'+'R0='+str(int(self.fi0*0.5*(self.kmkpc/self.plx)))+'+-'+str(int(self.sigmafi0*0.5*(self.kmkpc/self.plx)))+' km = '+str("{:.2f}".format(self.fi0*0.5*(self.kmkpc/self.plx)/695700))+'+-'+str("{:.2f}".format(self.sigmafi0*0.5*(self.kmkpc/self.plx)/695700))+'Rsun\n\nSigma='+str((180.*3600000./np.pi)*self.sigma)+'\n'+str(self.fi0*180*3600000/np.pi)
		#except:
		#	t='object: '+self.name+'\n\nprojection factor ='+str("{:.3f}".format(self.p))+'+-'+str("{:.3f}".format(self.sigmap))+'\n'+'R0='+str(int(self.fi0*0.5*(self.kmkpc/self.plx)))+'+-'+str(int(self.sigmafi0*0.5*(self.kmkpc/self.plx)))+' km = '+str("{:.2f}".format(self.fi0*0.5*(self.kmkpc/self.plx)/695700))+'+-'+str("{:.2f}".format(self.sigmafi0*0.5*(self.kmkpc/self.plx)/695700))+'Rsun\n\nSigma='+str((180.*3600000./np.pi)*self.sigma)+'\n'+str(self.fi0*180*3600000/np.pi)

		#print t
		return t
