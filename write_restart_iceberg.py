#!/usr/bin/env python

#First import the netcdf4 library
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import numpy as np  # http://code.google.com/p/netcdf4-python/
import matplotlib
#import matplotlib.pyplot as plt
import math
import os
matplotlib.use("GTKAgg")
import matplotlib.pyplot as plt
from pylab import *
import pdb



def Create_iceberg_restart_file(f,g,Number_of_bergs, Number_of_bonds,lon,lat,thickness,width,mass,mass_scaling,iceberg_num):
	
	# To copy the global attributes of the netCDF file  

	for attname in f.ncattrs():
		    setattr(g,attname,getattr(f,attname))


	# To copy the dimension of the netCDF file
	for dimname,dim in f.dimensions.iteritems():
		# if you want to make changes in the dimensions of the new file
		# you should add your own conditions here before the creation of the dimension.
		#g.createDimension(dimname,len(dim))
		g.createDimension(dimname,Number_of_bergs)

	# To copy the variables of the netCDF file

	for varname,ncvar in f.variables.iteritems():
		# if you want to make changes in the variables of the new file
		# you should add your own conditions here before the creation of the variable.
		var = g.createVariable(varname,ncvar.dtype,ncvar.dimensions)
		#Proceed to copy the variable attributes
		for attname in ncvar.ncattrs():  
			setattr(var,attname,getattr(ncvar,attname))
		#Finally copy the variable data to the new created variable
		var[:] = ncvar[0]
		
		if varname=='i':
			var[:]=Number_of_bergs
		
		if varname=='iceberg_num':
			for j in range(Number_of_bergs):
				#var[j]=j+1
				var[j]=iceberg_num[j]

		if varname=='uvel' or varname=='vvel' or varname=='uvel_old' or varname=='vvel_old' or varname=='axn' or varname=='ayn'\
		or varname=='bxn' or varname=='byn' or  varname=='halo_berg' or varname=='heat_density' or varname=='lon_old' or varname=='lat_old' \
		or varname=='mass_of_bits' or varname=='start_mass' or  varname=='start_day' or varname=='start_year' or varname=='start_lon' \
		or varname=='start_lat' or varname=='start_mass' or  varname=='start_day' or varname=='start_year' or varname=='start_lon' or varname=='lat_old' :
			var[:]=0

		if varname=='mass_scaling':
			var[:]=mass_scaling

		if varname=='thickness':
			var[:]=thickness
		
		if varname=='mass':
			var[:]=mass

		if varname=='width'  or varname=='length':
			var[:]=width

		if varname=='lon':
			for j in range(Number_of_bergs):
				var[j]=lon[j]

		if varname=='lat':
			for j in range(Number_of_bergs):
				var[j]=lat[j]


	f.close()
	g.close()


def Create_bond_restart_file(q,h, Number_of_bonds,first_berg_num,first_berg_ine,first_berg_jne,other_berg_ine,other_berg_jne,iceberg_num,other_berg_num):
	#Creating the bond restart file


	# To copy the global attributes of the netCDF file  

	for attname in h.ncattrs():
		    setattr(q,attname,getattr(h,attname))


	# To copy the dimension of the netCDF file
	for dimname,dim in h.dimensions.iteritems():
		# if you want to make changes in the dimensions of the new file
		# you should add your own conditions here before the creation of the dimension.
		#g.createDimension(dimname,len(dim))
		q.createDimension(dimname,Number_of_bonds)

	# To copy the variables of the netCDF file

	for varname,ncvar in h.variables.iteritems():
		# if you want to make changes in the variables of the new file
		# you should add your own conditions here before the creation of the variable.
		var = q.createVariable(varname,ncvar.dtype,ncvar.dimensions)
		#Proceed to copy the variable attributes
		for attname in ncvar.ncattrs():  
			setattr(var,attname,getattr(ncvar,attname))
		#Finally copy the variable data to the new created variable
		#var[:] = ncvar[0]
		var[:] = 0.

		if varname=='i':
			var[:]=Number_of_bonds

		if varname=='first_berg_num':
			for j in range(Number_of_bonds):
				var[j]=first_berg_num[j]

		if varname=='first_berg_ine':
			for j in range(Number_of_bonds):
				var[j]=first_berg_ine[j]

		if varname=='first_berg_jne':
			for j in range(Number_of_bonds):
				var[j]=first_berg_jne[j]

		if varname=='other_berg_num':
			for j in range(Number_of_bonds):
				var[j]=other_berg_num[j]

		if varname=='other_berg_ine':
			for j in range(Number_of_bonds):
				var[j]=other_berg_ine[j]

		if varname=='other_berg_jne':
			for j in range(Number_of_bonds):
				var[j]=other_berg_jne[j]

	h.close()
	q.close()



def plotting_iceberg_positions_and_bonds(lat,lon,first_berg_lat,first_berg_lon,other_berg_lat,other_berg_lon,Number_of_bergs,Number_of_bonds,R_earth,Radius,IA_scaling):

	Radius=Radius*IA_scaling
	circ_ind=np.linspace(0,2*pi,100);
	for k in range(Number_of_bergs):
		dR_lat=(Radius/R_earth)*(180/np.pi)
		dR_lon=(Radius/R_earth)*(180/np.pi)*(1/np.cos(lat[k]*np.pi/180))
		plt.plot(lon[k], lat[k],'bo-',linewidth=5)
		plt.plot(lon[k]+(dR_lon*cos(circ_ind)),lat[k]+(dR_lat*sin(circ_ind)),'b');

	for k in range(Number_of_bonds):
		x_bond=[]
		y_bond=[]
		x_bond.append(first_berg_lon[k])
		x_bond.append(other_berg_lon[k])
		y_bond.append(first_berg_lat[k])
		y_bond.append(other_berg_lat[k])
		plt.plot(x_bond, y_bond,'r',linewidth=5)


	#plt.plot(lon, lat,'bo-')
	plt.xlabel('longitude (deg)')
	plt.ylabel('latitude (deg)')
	plt.title('Iceberg initial positions')
	plt.grid(True)
	plt.show()




def Define_iceberg_positions(N,M ,lon_init,lat_init,Radius,R_earth):
	dx_berg=[]  #x distance in cartesian of berg from lon_init
	dy_berg=[]  #y distance in cartesian of berg from lat_init
	dlon_berg=[] #x distance in lon of berg from lon_init
	dlat_berg=[] #y distance in lat of berg from lat_init
	lon=[] #Longitude of iceberg
	lat=[] #Latitude of iceberg
	iceberg_num=[] #ID of iceberg

	berg_count=0
	for j in range(M):
		x_start=(j%2)*Radius
		y_start=np.sqrt(3)*Radius*j

		for i in range(N):
			berg_count=berg_count+1
			iceberg_num.append(berg_count)
			dx_berg.append(x_start+(2*i*Radius))
			dy_berg.append(y_start)

	Number_of_bergs=berg_count

	#Defining lon lat positions:
	#dlat_berg=(180/np.pi)*(1/R_earth)*dy_berg
	#dlon_berg=(180/np.pi)*(1/R_earth)*(1/cos(lat_init*np.pi/180))*dx_berg
	for i in range(Number_of_bergs):
		#Finding latittude
		dlat_dy=(180/np.pi)*(1/R_earth)
		dlat_berg.append(dlat_dy*dy_berg[i])
		lat.append(lat_init+dlat_berg[i])

		#Finding longitude
		dlon_dx=(180/np.pi)*(1/R_earth)*(1/np.cos(lat[i]*np.pi/180)) #Note that this depends on the latitude of the iceberg. Could approx this with lat_init.
		dlon_berg.append(dlon_dx*dx_berg[i])
		lon.append(lon_init+dlon_berg[i])

	return (Number_of_bergs,lon,lat,iceberg_num,dx_berg,dy_berg)

def Define_iceberg_bonds(Number_of_bergs,iceberg_num,lat,lon,dx_berg, dy_berg,Radius):
	
	#Defining Bonds:
	Bond=np.zeros((Number_of_bergs, Number_of_bergs))
	first_berg_num=[]  # Initializing bond list first berg
	first_berg_ine=[]  # Initializing bond list first berg
	first_berg_jne=[]  # Initializing bond list first berg
	first_berg_lat=[]  # Initializing bond list first berg
	first_berg_lon=[]  # Initializing bond list first berg
	other_berg_num=[]  # Initializing bond list other berg
	other_berg_ine=[]  # Initializing bond list other berg
	other_berg_jne=[]  # Initializing bond list other berg
	other_berg_lat=[]  # Initializing bond list other berg
	other_berg_lon=[]  # Initializing bond list other berg
	bond_count=0
	for i in range(Number_of_bergs):
		for j in range(Number_of_bergs):
			if i!=j:
				R_dist=np.sqrt(((dx_berg[i]-dx_berg[j])**2) + ((dy_berg[i]-dy_berg[j])**2))
				if R_dist < (2.01*Radius):
					bond_count=bond_count+1
					Bond[i,j]=1 # This is not really used.
					first_berg_num.append(iceberg_num[i])
					first_berg_ine.append(999)
					first_berg_jne.append(999)
					other_berg_num.append(iceberg_num[j])
					other_berg_ine.append(999)
					other_berg_jne.append(999)
					#Also get the lat lon of each for plotting
					first_berg_lat.append(lat[i])
					first_berg_lon.append(lon[i])
					other_berg_lat.append(lat[j])
					other_berg_lon.append(lon[j])
	Number_of_bonds=bond_count

	return [ Number_of_bonds, first_berg_num,first_berg_ine,first_berg_jne,first_berg_lat,first_berg_lon, other_berg_num,other_berg_ine, other_berg_jne,other_berg_lat,other_berg_lon]

####################################################################################################################################################
##########################################################  Main Program   #########################################################################
####################################################################################################################################################


def main():
	# Read en existing NetCDF file and create a new one
	# f is going to be the existing NetCDF file from where we want to import data
	# and g is going to be the new file.

	f=Dataset('input_files/icebergs.res.nc','r') # r is for read only
	g=Dataset('output_files/New_icebergs.res.nc','w') # w if for creating a file
					      # if the file exists it the file will be deleted to write on it

	h=Dataset('input_files/bonds_iceberg.res.nc','r') # r is for read only
	q=Dataset('output_files/New_bonds_iceberg.res.nc','w') # w if for creating a file



	#Parameters
	thickness=100.
	Radius=3*1000.
	rho_ice=850.
	mass_scaling=1.
	R_earth=6360.*1000.

	#Derived quanities
	width=np.sqrt(np.pi*(Radius**2))
	mass=thickness*rho_ice*np.pi*(Radius**2)

	#Let interactive radius be different from radius for testing the model:
	IA_scaling=1.#(1./2.)
	Radius=Radius/IA_scaling

	#Here we create the lons and lats for a tabular iceberg
	N= 5  # Number of rows in iceberg
	M= 5   # Number of columns in iceberg
	lon_init=-32.9  #longitude of bottom left corner of iceberg
	lat_init=-70.  #latitude  of bottom left corner of iceberg

	#Define the positions of the icebergs
	(Number_of_bergs,lon,lat,iceberg_num,dx_berg, dy_berg)= Define_iceberg_positions(N,M ,lon_init,lat_init,Radius,R_earth)

	#Define the positions of the iceberg bonds
	(Number_of_bonds, first_berg_num,first_berg_ine,first_berg_jne,first_berg_lat,first_berg_lon, other_berg_num,other_berg_ine, other_berg_jne,other_berg_lat,other_berg_lon)=\
			Define_iceberg_bonds(Number_of_bergs,iceberg_num,lat,lon,dx_berg, dy_berg,Radius)



	#####################################################################################
	#Creating iceberg restart file
	Create_iceberg_restart_file(f,g,Number_of_bergs, Number_of_bonds,lon,lat,thickness,width,mass,mass_scaling,iceberg_num)
	
	#Creating bond restart file
	Create_bond_restart_file(q,h, Number_of_bonds,first_berg_num,first_berg_ine,first_berg_jne,other_berg_ine,other_berg_jne,iceberg_num,other_berg_num)

	##################################################################################
	
	plotting_iceberg_positions_and_bonds(lat,lon,first_berg_lat,first_berg_lon,other_berg_lat,other_berg_lon,Number_of_bergs,Number_of_bonds,R_earth,Radius,IA_scaling)
	# Plotting the positions and bonds of the newly formed formation


	print 'Script complete'



if __name__ == '__main__':
	main()
	#sys.exit(main())














