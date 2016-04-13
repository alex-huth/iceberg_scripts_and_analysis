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


def Create_iceberg_restart_file(f,g,Number_of_bergs, Number_of_bonds,lon,lat,thickness,width,mass,mass_scaling,iceberg_num):
	
	print 'Writing iceberg restart files...'
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
			for j in range(Number_of_bergs):
				var[j]=thickness[j]
		
		if varname=='mass':
			for j in range(Number_of_bergs):
				var[j]=mass[j]

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


	print 'Writing bond restart files...'
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



def Define_iceberg_thickness_and_mass(Number_of_bergs,dx_berg,dy_berg,h_ice_vec,xi,yi,rho_ice,Radius):
	thickness=[]
	mass=[]
	for berg_count in range(Number_of_bergs):
		x_val=dx_berg[berg_count] ; y_val=dy_berg[berg_count]
		#y_ind=(abs(y-y_val)).argmin()
		#x_ind=(abs(x-x_val)).argmin()
		R_ind=(abs(xi-x_val)+(abs(yi-y_val))).argmin()
		
		#Interpolate thickness from 4 corners - possibly do this later, but I worry about when you are between a shelf and non shelf piece.
		#dir=np.sign(x_val-x[x_ind])
		#
		#Th=h_ice[y_ind,x_ind]
		Th=h_ice_vec[R_ind]
		thickness.append(Th)
		mass.append(Th*rho_ice*np.pi*(Radius**2))
	return [thickness, mass]

def Create_distance_and_ice_mask_vector(x,y,ice_mask,h_ice,input_is_cartesian,R_earth,lat_ref):
	N=len(x)  ; M = len(y)
	x_min=np.min(x)  ;y_min=np.min(y)
	count=-1
	#R=np.zeros([N*M,1])
	xi=np.zeros([N*M,1])
	yi=np.zeros([N*M,1])
	ice_mask_vec=np.zeros([N*M,1])
	h_ice_vec=np.zeros([N*M,1])
	x_temp=np.zeros([N*M,1])
	y_temp=np.zeros([N*M,1])

	for i in range(N):
		for j in range(M):
			count=count+1
			if input_is_cartesian==True:
				#R[count]=sqrt(((x[i]-x_min)**2) + ((y[j]-y_min)**2))
				xi[count]=x[i]-x_min
				yi[count]=y[j]-y_min
			else:
				#x and y are given in lon and lat
				dlon=x[i]-x_min
				dlat=y[j]-y_min
				dy_dlat=(np.pi/180)*(R_earth)			
				dx_dlon=(np.pi/180)*(R_earth)*(np.cos(lat_ref*np.pi/180))
				dx=dlon*dx_dlon
				dy=dlat*dy_dlat
				#R[count]=sqrt((dx**2)+(dy**2))
				xi[count]=dx
				yi[count]=dy
			ice_mask_vec[count,0]=ice_mask[j,i]
			h_ice_vec[count,0]=h_ice[j,i]
			x_temp[count,0]=x[i]
			y_temp[count,0]=y[j]
			

	ice_mask_vec[np.where(np.isnan(ice_mask_vec))]=0
	#cNorm = mpl.colors.Normalize(vmin=-1., vmax=5.)
	#plt.pcolor(x,y,ice_mask,norm=cNorm,cmap='jet')
	#plt.scatter(x_temp, y_temp,c=ice_mask_vec[:,0],cmap='jet',norm=cNorm)
	plt.show()

	return [xi,yi, ice_mask_vec,h_ice_vec]


def check_if_it_is_in_domain(x_val,y_val,x_min,x_max,y_min,y_max,input_is_cartesian,R_earth,lat_ref):
	point_is_in_domain=True
	if input_is_cartesian==False:
		lat_ref=y_max 
		dlat_dy=(180/np.pi)*(1/R_earth)
		dlon_dx=(180/np.pi)*(1/R_earth)*(1/np.cos(lat_ref*np.pi/180))
		x_val=(x_val*dlon_dx)
		y_val=(y_val*dlat_dy)
	
	if (x_val > (x_max-x_min)) or (x_val < 0) or  (y_val > (y_max-y_min)) or (y_val < 0):
		point_is_in_domain=False
	
	return point_is_in_domain

def check_if_it_is_ice(x_val,y_val,xi,yi,ice_mask_vec,input_is_cartesian):
	#R_val=np.sqrt(((x_val)**2) +((y_val)**2))
	#R_ind=(abs(R-R_val)).argmin()
	R_ind= (abs(x_val-xi) + abs(yi-y_val)).argmin()
	if ice_mask_vec[R_ind]==1:
		it_is_ice=True
	else:
		it_is_ice=False
	return it_is_ice


def Create_icebergs(N,M ,lon_init,lat_init,Radius,R_earth, x, y,ice_mask,h_ice,Convert_to_lat_lon,rho_ice,input_is_cartesian):
	dx_berg=[]  #x distance in cartesian of berg from lon_init
	dy_berg=[]  #y distance in cartesian of berg from lat_init
	dlon_berg=[] #x distance in lon of berg from lon_init
	#dlat_berg=[] #y distance in lat of berg from lat_init
	lon=[] #Longitude of iceberg
	lat=[] #Latitude of iceberg
	iceberg_num=[] #ID of iceberg
	lat_ref=np.max(y)



	if input_is_cartesian==False:
		lat_init=np.min(y)
		lon_init=np.min(x)

	#Create a vector of distances from a minimum point.
	(xi,yi,ice_mask_vec,h_ice_vec)=Create_distance_and_ice_mask_vector(x,y,ice_mask,h_ice,input_is_cartesian,R_earth,lat_ref)

	X_min=np.min(x); X_max=np.max(x)  #In lat lon or cartesian
	Y_min=np.min(y); Y_max=np.max(y)  #In lat lon or cartesian

	#N=2*int(ceil((R_max)/(2*Radius))+2)
	
	N=int(ceil((np.max(xi))/(Radius/2))+2)
	M=int(ceil((np.max(yi))/(Radius/2))+2)

	berg_count=0
	#for j in range(M):
	for j in range(N):
		x_start=((j%2)*Radius)
		y_start=(np.sqrt(3)*Radius*j)
		for i in range(M):
			x_val=x_start+(2*i*Radius)  ; y_val=y_start
			if check_if_it_is_in_domain(x_val,y_val,X_min,X_max,Y_min,Y_max,input_is_cartesian,R_earth,lat_ref):
			#if True:
				#R_val=np.sqrt(((x_val-x0)**2) +((y_val-y0)**2))
				if check_if_it_is_ice(x_val,y_val,xi,yi,ice_mask_vec,input_is_cartesian):
				#if True:
					berg_count=berg_count+1
					iceberg_num.append(berg_count)
					dx_berg.append(x_val)
					dy_berg.append(y_val)
		

	#print np.max(dx_berg),np.max(dy_berg)
	#plt.scatter(dx_berg, dy_berg)
	#plt.show()

	Number_of_bergs=berg_count
	print 'Number of bergs = ', Number_of_bergs
	
	#Defining the thickness of the icebergs
	(thickness, mass)=Define_iceberg_thickness_and_mass(Number_of_bergs,dx_berg,dy_berg,h_ice_vec,xi,yi,rho_ice,Radius)
	width=np.sqrt(np.pi*(Radius**2))
	

	if Convert_to_lat_lon==True:
		#Defining lon lat positions:
		#dlat_berg=(180/np.pi)*(1/R_earth)*dy_berg
		#dlon_berg=(180/np.pi)*(1/R_earth)*(1/cos(lat_init*np.pi/180))*dx_berg
		for i in range(Number_of_bergs):
			#Finding latittude
			dlat_dy=(180/np.pi)*(1/R_earth)
			lat.append(lat_init+(dlat_dy*dy_berg[i]))

			#Finding longitude
			dlon_dx=(180/np.pi)*(1/R_earth)*(1/np.cos(lat_ref*np.pi/180)) #Note that this depends on the latitude of the iceberg. Could approx this with lat_init.
			lon.append(lon_init+(dlon_dx*dx_berg[i] ))

			#dlon_berg.append(dlon_dx*dx_berg[i])
			#lon.append(lon_init+dlon_berg[i])
	else:
		lon=dx_berg+np.min(x)
		lat=dy_berg+np.min(y)

	return (Number_of_bergs,lon,lat,iceberg_num,dx_berg,dy_berg,thickness, mass,width)

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

def load_ISOMIP_ice_geometry(filename):
	with nc.Dataset(filename) as file:
		ocean_mask = file.variables['openOceanMask'][:,:]
		upperSurface = file.variables['upperSurface'][:,:]
		lowerSurface = file.variables['lowerSurface'][:,:]
		x = file.variables['x'][:]
		y = file.variables['y'][:]
		
		#x, y= np.meshgrid(x,y)

		ice_mask=1-ocean_mask #one if it ice, zero if it is ocean
		h_ice=upperSurface-lowerSurface #The ice thickness
	
	return [x,y,ice_mask,h_ice]


def load_Weddel_ice_geometry(filename):
	with nc.Dataset(filename) as file:
		h_ice = file.variables['thickness'][:,:]
		ice_mask = file.variables['icemask_grounded_and_shelves'][:,:]
		x = file.variables['lon'][:]
		y = file.variables['lat'][:]


		ice_mask[np.where(ice_mask<-1)]=0.
		#ice_mask=1-ocean_mask #one if it ice, zero if it is ocean
	
	return [x,y,ice_mask,h_ice]


def plotting_iceberg_positions_and_bonds(lat,lon,first_berg_lat,first_berg_lon,other_berg_lat,other_berg_lon,Number_of_bergs,Number_of_bonds,R_earth,Radius,IA_scaling,Convert_to_lat_lon,\
		plot_bonds,plot_circles,h_ice,ice_mask,x,y,plot_ice_mask,plot_ice_thickness,thickness,plot_icebergs_positions):
	print 'Starting to plot...'	
	Radius=Radius*IA_scaling
	circ_ind=np.linspace(0,2*pi,100);

	if plot_ice_mask==True:
		cNorm = mpl.colors.Normalize(vmin=0., vmax=1.)
		plt.pcolor(x,y,ice_mask,norm=cNorm)
	elif plot_ice_thickness==True:
		cNorm = mpl.colors.Normalize(vmin=0., vmax=1000.)
		plt.pcolor(x,y,h_ice,cmap='jet',norm=cNorm)


	if plot_icebergs_positions==True:
		#plt.scatter(lon, lat,color='yellow')
		plt.scatter(lon, lat,c=thickness)

	#plt.plot(lon, lat,'bo-',linewidth=5)
	if plot_circles==True:
		for k in range(Number_of_bergs):
			if Convert_to_lat_lon==True:
				dR_lat=(Radius/R_earth)*(180/np.pi)
				dR_lon=(Radius/R_earth)*(180/np.pi)*(1/np.cos(lat[k]*np.pi/180))
				plt.plot(lon[k]+(dR_lon*cos(circ_ind)),lat[k]+(dR_lat*sin(circ_ind)),'b');
			else:
				plt.plot(lon[k]+(Radius*cos(circ_ind)),lat[k]+(Radius*sin(circ_ind)),'b');


	if plot_bonds==True:
		for k in range(Number_of_bonds):
			x_bond=[]
			y_bond=[]
			x_bond.append(first_berg_lon[k])
			x_bond.append(other_berg_lon[k])
			y_bond.append(first_berg_lat[k])
			y_bond.append(other_berg_lat[k])
			plt.plot(x_bond, y_bond,'r',linewidth=5)


	#plt.plot(lon, lat,'bo-')
	if Convert_to_lat_lon==True:
		plt.xlabel('longitude (deg)')
		plt.ylabel('latitude (deg)')
	else:
		plt.xlabel('longitude (m)')
		plt.ylabel('latitude (m)')
	plt.title('Iceberg initial positions')
	plt.grid(True)
	plt.show()




####################################################################################################################################################
##########################################################  Main Program   #########################################################################
####################################################################################################################################################

def main():

	#Flags
	save_restart_files=True
	Convert_to_lat_lon=False
	input_is_cartesian=True

	#Plotting flags
	plot_bonds=False
	plot_circles=False
	plot_ice_mask=False
	plot_ice_thickness=True
	plot_icebergs_positions=True

	#Which experiment
	#Ice_geometry_source='ISOMIP'     ; Convert_to_lat_lon=False       ; input_is_cartesian=True  
	Ice_geometry_source='Weddell'     ; Convert_to_lat_lon=True        ; input_is_cartesian=False


	#Input and output files	
	f=Dataset('input_files/icebergs.res.nc','r') # r is for read only
	g=Dataset('output_files/shape_icebergs.res.nc','w') # w if for creating a file
	h=Dataset('input_files/bonds_iceberg.res.nc','r') # r is for read only
	q=Dataset('output_files/shape_bonds_iceberg.res.nc','w') # w if for creating a file

	ISOMIP_ice_geometry_filename='input_files/Isomip_ice_geometry.nc'
	Weddell_ice_geometry_filename='input_files/Bedmap2_gridded_subset_Weddell_Sea_Region.nc'


	#Parameters
	thickness=100.
	#Radius=2*1000
	Radius=20*1000
	rho_ice=850.
	mass_scaling=1.
	R_earth=6360.*1000.

	#Let interactive radius be different from radius for testing the model:
	IA_scaling=1.#(1./2.)
	Radius=Radius/IA_scaling

	#Here we create the lons and lats for a tabular iceberg
	N= 5  # Number of rows in iceberg
	M= 4   # Number of columns in iceberg
	lon_init=-32.9  #longitude of bottom left corner of iceberg
	lat_init=-70.  #latitude  of bottom left corner of iceberg

	#####################################################################################
	#####################################################################################

	if  Ice_geometry_source=='ISOMIP':
		(x,y,ice_mask,h_ice)=load_ISOMIP_ice_geometry(ISOMIP_ice_geometry_filename)
	
	if  Ice_geometry_source=='Weddell':
		(x,y,ice_mask,h_ice)=load_Weddel_ice_geometry(Weddell_ice_geometry_filename)


	#Define the positions,thickness, mass,  of the icebergs
	(Number_of_bergs,lon,lat,iceberg_num,dx_berg, dy_berg,thickness,mass,width)= Create_icebergs(N,M ,lon_init,lat_init,Radius,R_earth\
			,x,y,ice_mask,h_ice,Convert_to_lat_lon,rho_ice,input_is_cartesian)

	#Define the positions of the iceberg bonds
	(Number_of_bonds, first_berg_num,first_berg_ine,first_berg_jne,first_berg_lat,first_berg_lon, other_berg_num,other_berg_ine, other_berg_jne,other_berg_lat,other_berg_lon)=\
			Define_iceberg_bonds(Number_of_bergs,iceberg_num,lat,lon,dx_berg, dy_berg,Radius)


	if save_restart_files==True:	
		#Creating iceberg restart file
		Create_iceberg_restart_file(f,g,Number_of_bergs, Number_of_bonds,lon,lat,thickness,width,mass,mass_scaling,iceberg_num)
		
		#Creating bond restart file
		Create_bond_restart_file(q,h, Number_of_bonds,first_berg_num,first_berg_ine,first_berg_jne,other_berg_ine,other_berg_jne,iceberg_num,other_berg_num)

	
	# Plotting the positions and bonds of the newly formed formation
	plotting_iceberg_positions_and_bonds(lat,lon,first_berg_lat,first_berg_lon,other_berg_lat,other_berg_lon,Number_of_bergs,Number_of_bonds,R_earth,Radius,IA_scaling,Convert_to_lat_lon,\
			plot_bonds,plot_circles,h_ice,ice_mask,x,y,plot_ice_mask,plot_ice_thickness,thickness,plot_icebergs_positions)

#
	print 'Script complete'



if __name__ == '__main__':
	main()
	#sys.exit(main())














