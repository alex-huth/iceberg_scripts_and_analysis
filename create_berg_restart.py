#!/usr/bin/env python

###############################################################################################
from netCDF4 import Dataset
import numpy as np
import os
#import matplotlib
#import matplotlib.pyplot as plt
#from pylab import *
#import pdb
#import argparse
#This combination allows you to access the varibale imputted from the command line.
#import sys

#Clear screen
os.system('clear')
empty_iceberg_path = 'iceberg.res.nc'
calving_time=nc.Dataset(empty_iceberg_path).variables['i']


#Creating the calving file
fnam='iceberg.res.nc'
f=nc.Dataset(fnam,'w',format='NETCDF3_CLASSIC')
i=f.createDimension('i', None)
#yt=f.createDimension('yt',ny)
#xt=f.createDimension('xt',nx)

#lon=f.createVariable('lon','f4',('time','yt','xt'))
#lon=f.createVariable('lon','f4',('time','yt','xt'))
#time=f.createVariable('time','f8',('time'))
#yt=f.createVariable('yt','f8',('yt'))
#xt=f.createVariable('xt','f8',('xt'))

i_list=[1]
f.variables['i'][:]=i_list


f.sync()
f.close()





print 'Script complete'
