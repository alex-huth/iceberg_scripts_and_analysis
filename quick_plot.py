#!/usr/bin/env python

#First import the netcdf4 library
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import numpy as np  # http://code.google.com/p/netcdf4-python/
import matplotlib
import math
import os
matplotlib.use("GTKAgg")
from pylab import *
#import matplotlib.pyplot as plt
import pdb
import netCDF4 as nc
import sys
import argparse


def load_data_from_file(filename,field,dim_num,layer_number):
	with nc.Dataset(filename) as file:
		if dim_num==2:
			data = file.variables[field][:,:]
		if dim_num==3:
			data = file.variables[field][:,:,:]
			data=np.squeeze(np.mean(data,axis=0))  #Mean over first variable
		if dim_num==4:
			data = file.variables[field][:,:,:,:]
			data=np.squeeze(np.mean(data,axis=0))  #Mean over first variable
			#data=np.squeeze((data[0,:,:,:]))  #Mean over first variable
			data=np.squeeze(data[layer_number,:,:])#Layer choice over second 
	print data.shape
	return data

def plot_data_field(data,field,vmin=None,vmax=None): 
	print 'Starting to plot...'	
	if vmin==None:
		vmin=np.min(data)
	else:
		vmin=int(vmin)

	if vmax==None:
		vmax=np.max(data)
	else:
		vmax=int(vmax)
	cmap='jet'
	print vmin	
	print vmax	

	cNorm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
	#cNorm = mpl.colors.Normalize(vmin=600, vmax=850)
	plt.pcolor(data,norm=cNorm,cmap=cmap)
	plt.colorbar()
	plt.grid(True)
	plt.title(field)

def plot_operated_fields(data1,data2,field,operation): 
	
	print 'Starting to plot...'	
	vmin=min(np.min(data1),np.min(data2))
	vmax=max(np.max(data1),np.max(data2))
	cmap='jet'
	cNorm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

	#Data from file 1:
	plt.subplot(3,1,1)
	plt.pcolor(data1,norm=cNorm,cmap=cmap)
	plt.colorbar()
	plt.title('File 1: ' + field)
	#Data from file 2:
	plt.subplot(3,1,2)
	plt.pcolor(data2,norm=cNorm,cmap=cmap)
	plt.colorbar()
	plt.title('File 2: ' + field)


	#Data from difference:
	data3=data1-data2 #subtraction is the default
	if operation=='divide':
		data3=data1/data2
	if operation=='relative':
		data3=(data1-data2)/data1
	vmax=np.max(abs(data3))
	vmin=-vmax
	cmap='bwr'
	cNorm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
	if operation=='divide':
		cNorm = mpl.colors.Normalize(vmin=0, vmax=2)
	if operation=='relative':
		cNorm = mpl.colors.Normalize(vmin=-1, vmax=1)
	plt.subplot(3,1,3)
	plt.pcolor(data3,norm=cNorm,cmap=cmap)
	plt.colorbar()
	plt.title('Difference')

####################################################################################################################################################
##########################################################  Main Program   #########################################################################
####################################################################################################################################################

def main():
	parser = argparse.ArgumentParser()
        parser.add_argument('--file1', default=None, help='The input data file1 in NetCDF format.')
        parser.add_argument('--file2', default=None, help='The input data file2 in NetCDF format.')
        parser.add_argument('--field', default=None, help='Feild to plot')
        parser.add_argument('--operation', default=None, help='Operation betweeen fields')
        parser.add_argument('--layer_number', default=0, help='Operation betweeen fields')
        parser.add_argument('--dim_num', default=2, help='Number of dimensions of data')
        parser.add_argument('--vmax', default=None, help='Maximum for plotting')
        parser.add_argument('--vmin', default=None, help='Minimum for plotting')
	args = parser.parse_args()

	#Converting input
	field=args.field
	layer_number=int(args.layer_number)
	dim_num=int(args.dim_num)
	operation=args.operation
	filename1=args.file1
	print 'File 1 = ' , filename1

	#filename3='../../ocean_only/Alistair_ISOMIP/rho/ISOMIP_IC.nc'
	#filename1='/lustre/f1/unswept/Alon.Stern/MOM6-examples_Alon/ocean_only/Alistair_ISOMIP/rho/ISOMIP.nc'

	#if filename1==filename3:
	#	print 'BOOM'
	#else:
	#	print 'NO'
	plt.figure(figsize=(15,10))
	#fig = plt.figure(1)

	#Loading data from file1
	data1=load_data_from_file(filename1,field,dim_num,layer_number)

	if args.file2!=None:
		filename2=args.file2
		print 'File 2 = ' ,filename2
		#Loading data from file1
		data2=load_data_from_file(filename2,field,dim_num,layer_number)

		if operation=='subtract':
			print 'Subtracting fields'
			#data=data-data2
			#operation='subtract'

	if args.file2==None:
		plot_data_field(data1,field,args.vmin,args.vmax)	
	else:
		plot_operated_fields(data1,data2,field,operation)	


	#Plotting flags
	#fig.set_size_inches(9,4.5)
	plt.show()
	print 'Script complete'



if __name__ == '__main__':
	main()
	#sys.exit(main())














