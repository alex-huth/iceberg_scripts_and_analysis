#!/usr/bin/env python

###############################################################################################
from netCDF4 import Dataset
import numpy as np
import os
import numpy
import netCDF4 as nc
#import matplotlib
#import matplotlib.pyplot as plt
#from pylab import *
#import pdb
#import argparse
#This combination allows you to access the varibale imputted from the command line.
#import sys

#Clear screen
def main():

	os.system('clear')
	input_geometry_filename='input_files/Isomip_ocean_geometry.nc'

	#Flags
	use_flat_iceshelf=False

	if use_flat_iceshelf==True:
		new_filename='output_files/topog_flat.nc'
	else:
		new_filename='output_files/topog.nc'

	#f=Dataset(input_geometry_filename,'r')
	#g=Dataset(new_filename,'w') # w if for creating a file


	with nc.Dataset(input_geometry_filename) as file:
		Depth = file.variables['D'][:,:]

	M= Depth.shape
	ny=M[0]
	nx=M[1]
	print nx,ny
	
	#Creating the topog file
	g=Dataset(new_filename,'w') # w if for creating a file

	yt=g.createDimension('yt',ny)
	xt=g.createDimension('xt',nx)

	depth_h=g.createVariable('depth','f4',('yt','xt'))
	if use_flat_iceshelf==True:
		g.variables['depth'][:]=Depth*0
	else:
		g.variables['depth'][:]=Depth
	
	
	g.sync()
	g.close()





	print 'Script complete'



if __name__ == '__main__':
        main()
        #sys.exit(main())
