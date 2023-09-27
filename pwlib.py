#!/usr/bin/env python
import sys, os
import numpy as np
import math

def metnk(x,y):
	nump=0
	sx=0.
	sy=0.
	sxx=0.
	syy=0.
	sxy=0.
	for i,iks in enumerate(x):
		nump = nump + 1
		sx=sx+x[i]
		sy=sy+y[i]
		sxx=sxx+x[i]*x[i]
		syy=syy+y[i]*y[i]
		sxy=sxy+x[i]*y[i]
		
	s=nump
	delta=s*sxx-sx*sx
	a=(s*sxy-sx*sy)/delta
	b=(sxx*sy-sx*sxy)/delta
	
	suma=0
	for i,iks in enumerate(x):
		suma=suma+(y[i]-b-a*x[i])*(y[i]-b-a*x[i])

	sigmay=math.sqrt(suma/(float(nump)-2.))
	sigmaa=math.sqrt(s*(sigmay*sigmay/delta))
	sigmab=math.sqrt(sigmay*sigmay*sxx/delta)
#	sigmaa=pow(s*suma/(delta*(s-1.)),0.5);
#	sigmab=pow(sigmaa*sigmaa*sxx/s,0.5);
#	sigmay=pow(suma/((double)nump-1.),0.5);
	return a,b,sigmaa,sigmab,sigmay

def metnk_uncer(x,y,yerr):
	nump=0
	sx=0.
	sy=0.
	sxx=0.
	syy=0.
	sxy=0.
	for i,iks in enumerate(x):
		nump = nump + 1/(yerr[i]*yerr[i])
		sx=sx+x[i]/(yerr[i]*yerr[i])
		sy=sy+y[i]/(yerr[i]*yerr[i])
		sxx=sxx+x[i]*x[i]/(yerr[i]*yerr[i])
		syy=syy+y[i]*y[i]/(yerr[i]*yerr[i])
		sxy=sxy+x[i]*y[i]/(yerr[i]*yerr[i])
		
	s=nump
	delta=s*sxx-sx*sx
	a=(s*sxy-sx*sy)/delta
	b=(sxx*sy-sx*sxy)/delta
	
	suma=0
	for i,iks in enumerate(x):
		suma=suma+(y[i]-b-a*x[i])*(y[i]-b-a*x[i])

	sigmay=math.sqrt(suma/(float(len(x))))
	sigmaa=math.sqrt(s/delta)
	sigmab=math.sqrt(sxx/delta)
	#sigmaa=math.sqrt(s*(sigmay*sigmay/delta))
	#sigmab=math.sqrt(sigmay*sigmay*sxx/delta)
#	sigmaa=pow(s*suma/(delta*(s-1.)),0.5);
#	sigmab=pow(sigmaa*sigmaa*sxx/s,0.5);
#	sigmay=pow(suma/((double)nump-1.),0.5);
	return a,b,sigmaa,sigmab,sigmay

def readfile(plik):
	dane = os.popen('cat ' + plik).read().split('\n')[:-1]
	data = []
	nrow = len(dane)
	ncol = 0
	for i,linia in enumerate(dane):
		line_data = linia.split()
		for j,item in enumerate(line_data):
			if i == 0:
				ncol = ncol+1
				data.append([item])

			else:
				data[j].append(item)

	return data,ncol,nrow

def parse_file(plik):
	parsed_file = []
	dane = os.popen('cat ' + plik).read().split('---END---')[:-1]
	for band in dane:
		lines = band.split('\n')[:-1]
		for i,line in enumerate(lines):

			if 'Data name' in line:
				name = line.split('Data name: ')[1]

			if 'from file' in line:
				file_name = line.split('from file: ')[1]

			if 'ID column' in line:
				idcol = line.split('ID column: ')[1]

			if 'RA column' in line:
				racol = line.split('RA column: ')[1]

			if 'DEC column' in line:
				decol = line.split('DEC column: ')[1]

			if 'Data color' in line:
				color = line.split('Data color: ')[1]

			if 'x input data' in line:
				funx = lines[i+1].split('function ')[1]
				xcol = lines[i+2].split('column ')[1]

			if 'y input data' in line:
				funy = lines[i+1].split('function ')[1]
				ycol = lines[i+2].split('column ')[1]
				yval1 = lines[i+3].split(' * ')[0].replace('\t','',1).replace('+','',1)
				ycol1 = lines[i+3].split('column ')[1]
				yval2 = lines[i+4].split(' * ')[0].replace('\t','',1).replace('+','',1)
				ycol2 = lines[i+4].split('column ')[1]
				yval3 = lines[i+5].replace('\t','',1).replace('+','',1)

			if 'x-shift' in line:
				xshift = line.split('x-shift= ')[1]

			if 'y-shift' in line:
				yshift = line.split('y-shift= ')[1]

			if 'min-x-value' in line:
				xmin =  line.split('min-x-value= ')[1]

			if 'max-x-value' in line:
				xmax =  line.split('max-x-value= ')[1]

			if 'min-y-value' in line:
				ymin =  line.split('min-y-value= ')[1]

			if 'max-y-value' in line:
				ymax =  line.split('max-y-value= ')[1]

			if 'clipping number' in line:
				nclip =  line.split('clipping number: ')[1]

			if 'Criterion' in line:
				try:
					critcol = line.split()[2]
					critcond = line.split()[3]
					critval = line.split()[4]
				except:
					critcol = ''
					critcond = ''
					critval = ''

			
		parsed_file.append([name,file_name,idcol,xcol,ycol,yval1,ycol1,yval2,ycol2,yval3,racol,decol,color,funx,funy,xshift,yshift,xmin,xmax,ymin,ymax,nclip,critcol,critcond,critval])
	return parsed_file
