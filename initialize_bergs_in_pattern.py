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
from hexagon_area import Divide_hexagon_into_4_quadrants_old
from hexagon_area import Hexagon_into_quadrants_using_triangles

def Create_iceberg_restart_file(Number_of_bergs, lon,lat,thickness,width,mass,mass_scaling,iceberg_num,Ice_geometry_source):
	
	print 'Writing iceberg restart files...'
	# To copy the global attributes of the netCDF file  

	#Input and output files	
	f=Dataset('input_files/icebergs.res.nc','r') # r is for read only
	g=Dataset('output_files/' + Ice_geometry_source + '_icebergs.res.nc','w') # w if for creating a file

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

	
def Create_bond_restart_file(Number_of_bonds,first_berg_num,first_berg_ine,first_berg_jne,other_berg_ine,other_berg_jne,iceberg_num,other_berg_num,Ice_geometry_source):
	#Creating the bond restart file

	print 'Writing bond restart files...'
	# To copy the global attributes of the netCDF file  

	#Input and output files	
	h=Dataset('input_files/bonds_iceberg.res.nc','r') # r is for read only
	q=Dataset('output_files/' + Ice_geometry_source + '_bonds_iceberg.res.nc','w') # w if for creating a file

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



def Define_iceberg_thickness_and_mass(Number_of_bergs,dx_berg,dy_berg,h_ice_vec,xi,yi,rho_ice,Radius,x_ind_vec,y_ind_vec,h_ice,x,y,\
		set_all_thicknesses_to_one,element_area,Interpolate_from_four_corners):
	thickness=[]
	mass=[]
	for berg_count in range(Number_of_bergs):
		x_val=dx_berg[berg_count] ; y_val=dy_berg[berg_count]
		#y_ind=(abs(y-y_val)).argmin()
		#x_ind=(abs(x-x_val)).argmin()
		R_ind=(abs(xi-x_val)+(abs(yi-y_val))).argmin()

		#Interpolate thickness from 4 corners - possibly do this later, but I worry about when you are between a shelf and non shelf piece.
		if Interpolate_from_four_corners==True:
			ny,nx=h_ice.shape
			i1=x_ind_vec[R_ind][0]
			j1=y_ind_vec[R_ind][0]
			if x_val>xi[R_ind]:
				i2=i1-1
			else:
				i2=i1+1
			if y_val>yi[R_ind]:
				j2=j1-1
			else:
				j2=j1+1
			
			if i2>=0 and i2 <nx and j2>=0 and j2<ny:
				A=(h_ice[j1,i1]*(abs(y_val-y[j1]))+(h_ice[j2,i1]*(abs(y_val-y[j2]))))/( abs(y_val-y[j1]) + abs(y_val-y[j2])  )
				B=(h_ice[j1,i2]*(abs(y_val-y[j1]))+(h_ice[j2,i2]*(abs(y_val-y[j2]))))/( abs(y_val-y[j1]) + abs(y_val-y[j2])  )
				Th=(A*(abs(x_val-x[i1]))+(B*(abs(x_val-x[i2]))))/( abs(x_val-x[i1]) + abs(x_val-x[i2])  )
				#Th=h_ice[j2,i2]
			else:
				Th=h_ice_vec[R_ind][0]

		else:
			Th=h_ice_vec[R_ind][0]
		if set_all_thicknesses_to_one==True:
			Th=1.
		thickness.append(Th)
		mass.append(Th*rho_ice*element_area) 
		#if element_type=='square': 
		#	mass.append(Th*rho_ice*((2*Radius)**2)) # For square elements
		#else:
		#	mass.append(Th*rho_ice*np.pi*(Radius**2))  #For circular elements
	if set_all_thicknesses_to_one==True:
		print 'Thickness set to 1'
	return [thickness, mass]

def Create_distance_and_ice_mask_vector(x,y,ice_mask,h_ice,input_is_cartesian,R_earth,lat_ref,adjust_lat_ref):
	N=len(x)  ; M = len(y)
	x_min=np.min(x)  ;y_min=np.min(y)
	count=-1
	#R=np.zeros([N*M,1])
	xi=np.zeros([N*M,1])
	yi=np.zeros([N*M,1])
	ice_mask_vec=np.zeros([N*M,1])
	h_ice_vec=np.zeros([N*M,1])
	x_ind_vec=np.zeros([N*M,1])
	y_ind_vec=np.zeros([N*M,1])

	print N,M
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
				if adjust_lat_ref==True:
					lat_ref=y[j]

				dx_dlon=(np.pi/180)*(R_earth)*(np.cos(lat_ref*np.pi/180))
				dx=dlon*dx_dlon
				dy=dlat*dy_dlat
				#R[count]=sqrt((dx**2)+(dy**2))
				xi[count]=dx
				yi[count]=dy
			ice_mask_vec[count,0]=ice_mask[j,i]
			h_ice_vec[count,0]=h_ice[j,i]
			x_ind_vec[count,0]=i
			y_ind_vec[count,0]=j
			

	ice_mask_vec[np.where(np.isnan(ice_mask_vec))]=0
	#cNorm = mpl.colors.Normalize(vmin=-1., vmax=5.)
	#plt.pcolor(x,y,ice_mask,norm=cNorm,cmap='jet')
	#plt.scatter(x_temp, y_temp,c=ice_mask_vec[:,0],cmap='jet',norm=cNorm)
	plt.show()

	return [xi,yi, ice_mask_vec,h_ice_vec,x_ind_vec,y_ind_vec]


def check_if_it_is_in_domain(x_val,y_val,x_min,x_max,y_min,y_max,input_is_cartesian,R_earth,lat_init,adjust_lat_ref,dx,dy):
	point_is_in_domain=True
	if input_is_cartesian==False:
		dlat_dy=(180/np.pi)*(1/R_earth)
		if adjust_lat_ref==True:
			lat_ref=lat_init+(dlat_dy*y_val)
		else:
			lat_ref=y_max 
		dlon_dx=(180/np.pi)*(1/R_earth)*(1/np.cos(lat_ref*np.pi/180))
		x_val=(x_val*dlon_dx)
		y_val=(y_val*dlat_dy)
	
	if (x_val >= (x_max-x_min+(dx))) or (x_val<= 0) or  (y_val >= (y_max-y_min+(dy))) or (y_val <= 0):
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

def calculate_element_area(element_type,Radius):
	if element_type=='square':
		element_area=(2*Radius)**2
	elif element_type=='hexagon':
		element_area=(3.*np.sqrt(3.)/2.)*((4./3.)*(Radius)**2) #Area of hexagon around circle (used for packing)  
		#Another derivation uses innner hexagon with two more triangles added, which is a 1/6 of the hexagon area each (two since there are 6, shared by 3 each)
		#element_area=(4./3.)*H_i, where H_i=(3.*np.sqrt(3.)/2.)*((Radius)**2)  is the area of the inner hexagon (with sides equal to the radius)

	return element_area

def Create_icebergs(lon_init,lat_init,Radius,R_earth, x, y,ice_mask,h_ice,Convert_to_lat_lon,rho_ice,input_is_cartesian,\
		element_type,scale_the_grid_to_lat_lon,lat_ref,adjust_lat_ref,set_all_thicknesses_to_one,Interpolate_from_four_corners):
	print 'Starting to create icebergs...'
	dx_berg=[]  #x distance in cartesian of berg from lon_init
	dy_berg=[]  #y distance in cartesian of berg from lat_init
	dlon_berg=[] #x distance in lon of berg from lon_init
	#dlat_berg=[] #y distance in lat of berg from lat_init
	lon=[] #Longitude of iceberg
	lat=[] #Latitude of iceberg
	iceberg_num=[] #ID of iceberg

	if input_is_cartesian==False:
		lat_init=np.min(y)
		lon_init=np.min(x)

	#Create a vector of distances from a minimum point.
	(xi,yi,ice_mask_vec,h_ice_vec,x_ind_vec,y_ind_vec)=Create_distance_and_ice_mask_vector(x,y,ice_mask,h_ice,input_is_cartesian,R_earth,lat_ref,adjust_lat_ref)
	
	#dx=x[1]-x[0]  ;dy=y[1]-y[0]
	#X_min=np.min(x)-dx; X_max=np.max(x)+dx  #In lat lon or cartesian
	#Y_min=np.min(y)-dy; Y_max=np.max(y)+dy  #In lat lon or cartesian
	X_min=np.min(x); X_max=np.max(x)  #In lat lon or cartesian
	Y_min=np.min(y); Y_max=np.max(y)  #In lat lon or cartesian

	#N=2*int(ceil((R_max)/(2*Radius))+2)
	
	N=int(ceil((np.max(xi))/(Radius/2))+2) +2
	M=int(ceil((np.max(yi))/(Radius/2))+2) +2
	#N=3
	#M=10

 	dx=x[1]-x[0]	
 	dy=y[1]-y[0]	

	berg_count=0
	#for j in range(N):
	for i in range(M):
		if element_type=='square':
			x_start=Radius+(2*Radius*i)
			y_start=(Radius)
			#y_start=(Radius)+(2*Radius*j)
		else:
			#x_start=(Radius)+(((j%2)*Radius))
			#y_start=(Radius)+(np.sqrt(3)*Radius*j)
			y_start=(Radius)+(((i%2)*Radius))
			x_start=(Radius)+(np.sqrt(3)*Radius*i)

		for j in range(N):
		#for i in range(M):
			#x_val=x_start+(2*i*Radius)  ; y_val=y_start
			y_val=y_start+(2*j*Radius)  ; x_val=x_start
			if check_if_it_is_in_domain(x_val,y_val,X_min,X_max,Y_min,Y_max,input_is_cartesian,R_earth,lat_init,adjust_lat_ref,dx,dy):
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
	print 'Icebergs created. Number of bergs = ', Number_of_bergs
	
	
	#width=np.sqrt(np.pi*(Radius**2))
	#width=2*Radius
	element_area=calculate_element_area(element_type,Radius)
	width=np.sqrt(element_area)
	
	#Defining the thickness of the icebergs
	(thickness, mass)=Define_iceberg_thickness_and_mass(Number_of_bergs,dx_berg,dy_berg,h_ice_vec,xi,yi,rho_ice,Radius,x_ind_vec,y_ind_vec, h_ice,x,y,\
			set_all_thicknesses_to_one,element_area,Interpolate_from_four_corners)	
	

	if Convert_to_lat_lon==True:
		#Defining lon lat positions:
		#dlat_berg=(180/np.pi)*(1/R_earth)*dy_berg
		#dlon_berg=(180/np.pi)*(1/R_earth)*(1/cos(lat_init*np.pi/180))*dx_berg
		for i in range(Number_of_bergs):
			#Finding latittude
			dlat_dy=(180/np.pi)*(1/R_earth)

			lat.append(lat_init+(dlat_dy*dy_berg[i]))

			#Finding longitude
			if adjust_lat_ref==True:
				lat_ref=lat_init+(dlat_dy*dy_berg[i])
				#print lat_ref
			dlon_dx=(180/np.pi)*(1/R_earth)*(1/np.cos(lat_ref*np.pi/180)) #Note that this depends on the latitude of the iceberg. Could approx this with lat_init.
			lon.append(lon_init+(dlon_dx*dx_berg[i] ))

			#dlon_berg.append(dlon_dx*dx_berg[i])
			#lon.append(lon_init+dlon_berg[i])
	else:
		if scale_the_grid_to_lat_lon==True:
			Scale_up=1./2000.
			Radius=Radius*Scale_up
			dx_berg = [(i*Scale_up) for i in dx_berg] ; dy_berg = [i*Scale_up for i in dy_berg] 
			x=x*Scale_up ; y=y*Scale_up
			dx=dx*Scale_up ;dy=dy*Scale_up
		x=(x-np.min(x))+(dx/2) ; y=(y-np.min(y))+(dy/2)
		lon=dx_berg  ; lat=dy_berg
	
	#print np.unique(lat)
	#print y
	return (Number_of_bergs,lon,lat,iceberg_num,dx_berg,dy_berg,thickness, mass,width,x,y,Radius)



def Create_calving_event(lat1,lon1,lat2,lon2):
	Calve_lat =20  
	Calve_lon=160
	R1=np.sqrt((lon1-Calve_lon)**2+ (lat1-Calve_lat)**2)
	R2=np.sqrt((lon2-Calve_lon)**2+ (lat2-Calve_lat)**2)
	R_calve=15
	bond_broken=False
	#if ((R1 < R_calve)*(R2<R_calve) )  < 0.5 :
	#	bond_broken=False
	#if ((R1 > R_calve)*(R2>R_calve) )  < 0.5 :
	#	bond_broken=False
	if ((R1 < R_calve)*(R2>R_calve) )  > 0.5 :
		bond_broken=True
	if ((R2 < R_calve)*(R1>R_calve) )  > 0.5 :
		bond_broken=True
	#if (R2 < R_calve)  < 0.5 :
	#	bond_broken=True
	return bond_broken

def Define_iceberg_bonds(Number_of_bergs,iceberg_num,lat,lon,dx_berg, dy_berg,Radius,break_some_bonds):
	print 'Starting to create bonds...'
	#Defining Bonds:
	Bond=np.zeros((Number_of_bergs, Number_of_bergs))
	bond_broken=False
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
		#for j in range(Number_of_bergs):
		for j in range(i):
			if i!=j:
				R_dist=np.sqrt(((dx_berg[i]-dx_berg[j])**2) + ((dy_berg[i]-dy_berg[j])**2))
				if break_some_bonds==True:
					bond_broken=Create_calving_event(lat[i],lon[i],lat[j],lon[j])
				if R_dist < (2.01*Radius) and (bond_broken==False):
					bond_count=bond_count+2
					#Connect bond in the first direction
					first_berg_num.append(iceberg_num[i]);	first_berg_ine.append(999); first_berg_jne.append(999);first_berg_lat.append(lat[i]); first_berg_lon.append(lon[i])
					other_berg_num.append(iceberg_num[j]); 	other_berg_ine.append(999); other_berg_jne.append(999);other_berg_lat.append(lat[j]); other_berg_lon.append(lon[j])
					#Connect bond in the other direction
					first_berg_num.append(iceberg_num[j]);	first_berg_ine.append(999); first_berg_jne.append(999);first_berg_lat.append(lat[j]); first_berg_lon.append(lon[j])
					other_berg_num.append(iceberg_num[i]); 	other_berg_ine.append(999); other_berg_jne.append(999);other_berg_lat.append(lat[i]); other_berg_lon.append(lon[i])
	Number_of_bonds=bond_count

	return [ Number_of_bonds, first_berg_num,first_berg_ine,first_berg_jne,first_berg_lat,first_berg_lon, other_berg_num,other_berg_ine, other_berg_jne,other_berg_lat,other_berg_lon]

def load_ISOMIP_ice_geometry(filename,buffer_number):
	with nc.Dataset(filename) as file:
		ocean_mask = file.variables['openOceanMask'][:,:]
		upperSurface = file.variables['upperSurface'][:,:]
		lowerSurface = file.variables['lowerSurface'][:,:]
		x = file.variables['x'][:]
		y = file.variables['y'][:]
	print np.max(x)	-np.min(x)
	print np.max(y)	-np.min(y)
	#x, y= np.meshgrid(x,y)

	ice_mask=1-ocean_mask #one if it ice, zero if it is ocean
	h_ice=upperSurface-lowerSurface #The ice thickness
	M=ice_mask.shape
	print M
	#Setting the boundaries to non-ice
	A=np.arange(0,buffer_number)
	B=np.arange(M[0]-buffer_number,M[0])
	C=np.arange(M[1]-buffer_number,M[1])
	print A, B, C
	#if buffer_number>0:
	ice_mask[A,:]=0; ice_mask[B,:]=0
	ice_mask[:,A]=0; ice_mask[:,C]=0
	#print x
	#print y
	
	return [x,y,ice_mask,h_ice]


def load_ISOMIP_reduced_ice_geometry(filename,buffer_number):
	with nc.Dataset(filename) as file:
		h_ice = file.variables['thick'][:,:]
		ice_mask=h_ice>0.
	M=h_ice.shape
	print M
	y=np.linspace(1000,79000,M[0],endpoint=True)
	x=np.linspace(321000,799000,M[1],endpoint=True)
	#print x
	#print y
	#Setting the boundaries to non-ice
	A=np.arange(0,buffer_number)
	B=np.arange(M[0]-buffer_number,M[0])
	C=np.arange(M[1]-buffer_number,M[1])
	print A, B, C
	#if buffer_number>0:
	ice_mask[A,:]=0; ice_mask[B,:]=0
	ice_mask[:,A]=0; ice_mask[:,C]=0
	
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

def Select_just_one_berg(lon,lat,thickness,width,mass,iceberg_num,chosen_berg_num):
	print 'You have chosen to choose just one icebergs!!!'
	for k in range(len(lat)):
		if iceberg_num[k]==chosen_berg_num:
			berg_ind=k
	lon_temp=lon[berg_ind]; lon=[] ; lon.append(lon_temp)
	lat_temp=lat[berg_ind]; lat=[] ; lat.append(lat_temp)
	thickness_temp=thickness[berg_ind]; thickness=[] ; thickness.append(thickness_temp)
	mass_temp=mass[berg_ind]; mass=[] ; mass.append(mass_temp)
	#width_temp=width[berg_ind]; width=[] ; width.append(width_temp)
	iceberg_num_temp=iceberg_num[berg_ind]; iceberg_num=[] ; iceberg_num.append(iceberg_num_temp)
	Number_of_bergs=1
	return [Number_of_bergs,lon,lat,thickness,width,mass,iceberg_num]

def plotting_iceberg_positions(lat,lon,Number_of_bergs,R_earth,Radius,IA_scaling,Convert_to_lat_lon, \
		plot_circles,h_ice,ice_mask,x,y,plot_ice_mask,plot_ice_thickness,thickness,plot_icebergs_positions):
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
		cNorm = mpl.colors.Normalize(vmin=0., vmax=1000.)
		plt.scatter(lon, lat,c=thickness,norm=cNorm,cmap='jet')

	#plt.plot(lon, lat,'bo-',linewidth=5)
	if plot_circles==True:
		for k in range(Number_of_bergs):
			if Convert_to_lat_lon==True:
				dR_lat=(Radius/R_earth)*(180/np.pi)
				dR_lon=(Radius/R_earth)*(180/np.pi)*(1/np.cos(lat[k]*np.pi/180))
				plt.plot(lon[k]+(dR_lon*cos(circ_ind)),lat[k]+(dR_lat*sin(circ_ind)),'b');
			else:
				plt.plot(lon[k]+(Radius*cos(circ_ind)),lat[k]+(Radius*sin(circ_ind)),'b');



	#plt.plot(lon, lat,'bo-')
	if Convert_to_lat_lon==True:
		plt.xlabel('longitude (deg)')
		plt.ylabel('latitude (deg)')
	else:
		plt.xlabel('longitude (m)')
		plt.ylabel('latitude (m)')
	plt.title('Iceberg initial positions')
	plt.grid(True)


def plotting_iceberg_bonds(first_berg_lat,first_berg_lon,other_berg_lat,other_berg_lon,Number_of_bonds):

	for k in range(Number_of_bonds):
		x_bond=[]
		y_bond=[]
		x_bond.append(first_berg_lon[k])
		x_bond.append(other_berg_lon[k])
		y_bond.append(first_berg_lat[k])
		y_bond.append(other_berg_lat[k])
		plt.plot(x_bond, y_bond,'r',linewidth=5)


def spread_mass_to_ocean(i,j,mass_on_ocean,x,y,Area,Mass,element_type,grid_area):
		#Note that the x,y coming into this routine are the position within a cell (from 0 to 1), with 0.5,0.5 being in the center of the cell.

		#Initialize weights for each cell	
		yDxL=0.  ; yDxC=0. ; yDxR=0. ; yCxL=0. ; yCxR=0. 
		yUxL=0.  ; yUxC=0. ; yUxR=0. ; yCxC=1.

		if element_type=='square':
		#if True:
			L = min(( (np.sqrt(Area)  / np.sqrt(grid_area))),1) ;  #Non dimensionalize element length by grid area. (This gives the non-dim length of the square)
			xL=min(0.5, max(0., 0.5-(x/L)))
			xR=min(0.5, max(0., (x/L)+(0.5-(1/L) )))
			xC=max(0., 1.-(xL+xR))
			yD=min(0.5, max(0., 0.5-(y/L)))
			yU=min(0.5, max(0., (y/L)+(0.5-(1/L) )))
			yC=max(0., 1.-(yD+yU))

			yDxL=yD*xL#*grd%msk[i-1,j-1]
			yDxC=yD*xC#*grd%msk[i  ,j-1]
			yDxR=yD*xR#*grd%msk[i+1,j-1]
			yCxL=yC*xL#*grd%msk[i-1,j  ]
			yCxR=yC*xR#*grd%msk[i+1,j  ]
			yUxL=yU*xL#*grd%msk[i-1,j+1]
			yUxC=yU*xC#*grd%msk[i  ,j+1]
			yUxR=yU*xR#*grd%msk(i+1,j+1]
			yCxC=1.-( ((yDxL+yUxR)+(yDxR+yUxL)) + ((yCxL+yCxR)+(yDxC+yUxC)) )

		if element_type=='hexagon':
		#if False:#element_type=='hexagon':
			H = min(( (np.sqrt(Area/(2*sqrt(3)))  / np.sqrt(grid_area))),1) ;  #Non dimensionalize element length by grid area. (This gives the non-dim Apothen of the hexagon)
			S=(2/np.sqrt(3))*H #Side of the hexagon
			if S>0.5:
				print 'Elements must be smaller than a whole gridcell', 'i.e.: S= ' , S , '>=0.5'
				halt

			#Subtracting the position of the nearest corner from x,y
			origin_x=1 ; origin_y=1
			if x<0.5:
				origin_x=0
			if y<0.5:
				origin_y=0

			x0=(x-origin_x) #Position of the hexagon center, relative to origin at the nearest vertex
			y0=(y-origin_y)
			
			#(Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)= Divide_hexagon_into_4_quadrants_old(x0,y0,H)
			(Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)= Hexagon_into_quadrants_using_triangles(x0,y0,H,0.)
			if min(min(Area_Q1,Area_Q2),min(Area_Q3, Area_Q4)) <0:
				print 'Yolo'
				print x0,y0,H
				#print min(min(Area_Q1,Area_Q2),min(Area_Q3, Area_Q4))
				print Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4

			Area_Q1=Area_Q1/Area_hex
			Area_Q2=Area_Q2/Area_hex
			Area_Q3=Area_Q3/Area_hex
			Area_Q4=Area_Q4/Area_hex

			#Now, you decide which quadrant belongs to which mass on ocean cell.
			#Top right vertex
			if x>=0.5 and y>= 0.5:
				yUxR=Area_Q1
				yUxC=Area_Q2
				yCxC=Area_Q3
				yCxR=Area_Q4
			
			#Top left vertex
			if x<0.5 and y>= 0.5:
				yUxC=Area_Q1
				yUxL=Area_Q2
				yCxL=Area_Q3
				yCxC=Area_Q4
			
			#Bottom left vertex
			if x<0.5 and y< 0.5:
				yCxC=Area_Q1
				yCxL=Area_Q2
				yDxL=Area_Q3
				yDxC=Area_Q4

			#Bottom right vertex
			if x>=0.5 and y< 0.5:
				yCxR=Area_Q1
				yCxC=Area_Q2
				yDxC=Area_Q3
				yDxR=Area_Q4

			#if Sector<4 and Sector>-1:
			#	print Sector
				
					
		#Check that this is true
		if abs(yCxC-(1.-( ((yDxL+yUxR)+(yDxR+yUxL)) + ((yCxL+yCxR)+(yDxC+yUxC)) )))>0.001:
			print 'All the mass is not being used!!!'
			#print W1 , W2 , W3 , W4 , W5 , W6
			#print Area_Upper, Area_Lower, Area_right, Area_left
			print 'Areas: ',Area_hex,Area_hex*Area_Q1, Area_hex*Area_Q2, Area_hex*Area_Q3, Area_hex*Area_Q4
			print 'x0=',x0, 'y0=',y0, 'H=', H
			print 'Total area= ',(Area_Q1+Area_Q2+Area_Q3+Area_Q4)#, Sector
			
			#if W1==False and W2==True:
			#	print W1, W2
			#	print 'Stop the party!'
			#	print W1, W2, W3, W4, W5, W5
			#	print x0, y0
			#	print H, S
			#	halt
			#print 'Total added= ',  yCxC, (1.-( ((yDxL+yUxR)+(yDxR+yUxL)) + ((yCxL+yCxR)+(yDxC+yUxC)) )) , abs(yCxC-(1.-( ((yDxL+yUxR)+(yDxR+yUxL)) + ((yCxL+yCxR)+(yDxC+yUxC)) )))


		mass_on_ocean[i,j,1]=mass_on_ocean[i,j,1]+yDxL*Mass
		mass_on_ocean[i,j,2]=mass_on_ocean[i,j,2]+yDxC*Mass
		mass_on_ocean[i,j,3]=mass_on_ocean[i,j,3]+yDxR*Mass
		mass_on_ocean[i,j,4]=mass_on_ocean[i,j,4]+yCxL*Mass
		mass_on_ocean[i,j,5]=mass_on_ocean[i,j,5]+yCxC*Mass
		mass_on_ocean[i,j,6]=mass_on_ocean[i,j,6]+yCxR*Mass
		mass_on_ocean[i,j,7]=mass_on_ocean[i,j,7]+yUxL*Mass
		mass_on_ocean[i,j,8]=mass_on_ocean[i,j,8]+yUxC*Mass
		mass_on_ocean[i,j,9]=mass_on_ocean[i,j,9]+yUxR*Mass



		return mass_on_ocean

def regrid_iceberg_thickness(lat,lon,Number_of_bergs,thickness,mass,h_ice,x,y,rho_ice,element_type):
	Nx=len(x)  ; Ny=len(y)
	dx=x[1]-x[0]    ;dy=y[1]-y[0]
	New_mass=h_ice*0 #Initializing regrided field to zero
	Area = [(mass[i]/(rho_ice*thickness[i])) for i in range(Number_of_bergs)] ;
	grid_area=dx*dy #Assuming a regular grid
	Orig_mass=h_ice*rho_ice*grid_area
	
	print grid_area, dx, dy

	mass_on_ocean=np.zeros([Nx,Ny,10])   #Setting up matrix to spread mass to ocean.  Note that I have used 10 points  so that I can ignore 0 and match with python numbering
	#mass_on_ocean=np.zeros([Nx,Ny])   #Setting up matrix to spread mass to ocean.
	#Note: You may want to import the grid that the ocean model actually sees.
	#print np.unique(lat)
	#print np.unique(lon)
	#print np.unique(lat)/dy
	#print np.unique(lon)/dy
	#return
	for berg_count in range(Number_of_bergs):  
		x_val=lon[berg_count]  
		y_val=lat[berg_count]
		Area_val=Area[berg_count]
		mass_val=mass[berg_count]

		#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		#These two lines are wrong!!!
		#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		#j_val=(abs(y-y_val)).argmin()
		#i_val=(abs(x-x_val)).argmin()
		j_val=floor(y_val/dy)
		i_val=floor(x_val/dx)
		#print i_val,j_val
		
		x_cell=(x_val-x[i_val])/sqrt(grid_area)+0.5
		y_cell=(y_val-y[j_val])/sqrt(grid_area)+0.5
		#print x_cell, y_cell
		#print x_val, y_val
		#print i_val,j_val


		mass_on_ocean=spread_mass_to_ocean(i_val,j_val,mass_on_ocean,x_cell,y_cell,Area_val,mass_val,element_type,grid_area)

		#T_gr[j,i]=T_gr[j,i]+h_ice[j,i]
	
	#Adding mass_onto_ocean  (skipping sides for now)
	for i in range(1,Nx-1):
		for j in range(1,Ny-1):
			New_mass[j,i]=mass_on_ocean[i,j,5]  \
					+  ( ( (mass_on_ocean[i-1,j-1,9]+mass_on_ocean[i+1,j+1,1])   \
					+  (mass_on_ocean[i+1,j-1,7]+mass_on_ocean[i-1,j+1,3]) ) \
					+   ( (mass_on_ocean[i-1,j  ,6]+mass_on_ocean[i+1,j  ,4])   \
					+  (mass_on_ocean[i  ,j-1,8]+mass_on_ocean[i  ,j+1,2]) ) )
	
	
	
	#print 'Warning: the interpolation scheme is not the same as the one used in the model'	
	#plot_data=(New_mass-Orig_mass)/(rho_ice*grid_area)
	#plot_data=(New_mass-Orig_mass)/grid_area
	
	#plot_data=(Orig_mass-New_mass)/Orig_mass
	#plot_data=h_ice
	plot_data=(New_mass)/(rho_ice*grid_area)
	#plot_data=(Orig_mass)/(rho_ice*grid_area)

	#Zero'ing out the sides since we are not yet doing the interpolation properly here.
	plot_data[0,:]=0. ;	plot_data[:,0]=0. ;
	plot_data[:,Nx-1]=0. ;	plot_data[Ny-1,:]=0. ;

	vmax=np.max(abs(plot_data))
	#cNorm = mpl.colors.Normalize(vmin=-1., vmax=1.)
	#cNorm = mpl.colors.Normalize(vmin=-vmax, vmax=vmax)
	#cNorm = mpl.colors.Normalize(vmin=0, vmax=vmax)
	cNorm = mpl.colors.Normalize(vmin=0.98, vmax=1.02)
        #plt.pcolor(x,y,h_ice-T_gr,cmap='bwr',norm=cNorm)
        plt.pcolor(x,y,plot_data,cmap='bwr',norm=cNorm)
        #plt.pcolor(x,y,plot_data,cmap='jet',norm=cNorm)
	plt.title('Error between bergs and shelf')
	plt.colorbar()
	#plt.show()



####################################################################################################################################################
##########################################################  Main Program   #########################################################################
####################################################################################################################################################

def main():

	#Flags
	save_restart_files=True
	Convert_to_lat_lon=False
	input_is_cartesian=True

	#Plotting flags
	plotting_bergs_over_grid=False
	plot_bonds=False
	plot_circles=False
	plot_ice_mask=False
	plot_ice_thickness=True
	plot_icebergs_positions=True
	Convert_axes_to_lat_lon=False
	only_choose_one_berg=False  ; chosen_berg_num=1
	Create_icebergs_bonds=False
	scale_the_grid_to_lat_lon=True  #Remember to change this back when providing fields for Isomip
	break_some_bonds=False
	adjust_lat_ref=True
	set_all_thicknesses_to_one=True
	Interpolate_from_four_corners=False
	regrid_icebergs_onto_grid=True
	Switch_regridding_element_type=False

	#element_type='square' #'hexagonal'
	element_type='hexagon'

	#Which experiment
	Ice_geometry_source='ISOMIP'   ; Convert_to_lat_lon=False       ; input_is_cartesian=True ; 
	#Ice_geometry_source='Weddell'     ; Convert_to_lat_lon=True        ; input_is_cartesian=False

	ISOMIP_reduced=True  #Reduced uses 2X2 grid, not reduced uses 1X1 grid

	ISOMIP_ice_geometry_filename='input_files/Isomip_ice_geometry.nc'
	ISOMIP_reduced_ice_geometry_filename='input_files/isomip_ice_shelf1.nc'
	Weddell_ice_geometry_filename='input_files/Bedmap2_gridded_subset_Weddell_Sea_Region.nc'


	#Parameters
	thickness=100.
	#Radius=0.25*1000
	#Radius=sqrt(3)/2.*1000
	#Radius=1.*1000
	Radius=0.5*1000.  #Hexagon only valid for S<half gridcell.  (about 0.85 using 2km)
	print 'Radius = ', Radius
	rho_ice=850.
	mass_scaling=1.
	R_earth=6360.*1000.
	buffer_number=0   #Number of buffer points on the sides of the domain

	#Let interactive radius be different from radius for testing the model:
	IA_scaling=1.#(1./2.)
	Radius=Radius/IA_scaling

	#Here we create the lons and lats for a tabular iceberg
	#N= 5  # Number of rows in iceberg
	#M= 4   # Number of columns in iceberg

	print 'Element type= ', element_type, ';  Switched_regridding= ', Switch_regridding_element_type

	#####################################################################################
	#####################################################################################

	if  Ice_geometry_source=='ISOMIP':
		if ISOMIP_reduced==False:
			lon_init=0  ; lat_init=-70.  #latitude  of bottom left corner of iceberg
			(x,y,ice_mask,h_ice)=load_ISOMIP_ice_geometry(ISOMIP_ice_geometry_filename,buffer_number)
		else:
			lon_init=0  ; lat_init=-70.  #latitude  of bottom left corner of iceberg
			(x,y,ice_mask,h_ice)=load_ISOMIP_reduced_ice_geometry(ISOMIP_reduced_ice_geometry_filename,buffer_number)
	
	if  Ice_geometry_source=='Weddell':
		lon_init=-32.9  ; lat_init=-70.  #latitude  of bottom left corner of iceberg
		(x,y,ice_mask,h_ice)=load_Weddel_ice_geometry(Weddell_ice_geometry_filename)
	
	lat_ref=np.max(y)
	if (Convert_to_lat_lon==True)  and (input_is_cartesian==True):
		Convert_axes_to_lat_lon=True
		lat_ref=lat_init
		x=x-np.min(x)
		y=y-np.min(y)
		x=lon_init+((x/R_earth)*(180./np.pi)*(1/np.cos(lat_ref*np.pi/180.)))
		y=lat_init+((y/R_earth)*(180./np.pi))
		input_is_cartesian=False


	#Define the positions,thickness, mass,  of the icebergs
	(Number_of_bergs,lon,lat,iceberg_num,dx_berg, dy_berg,thickness,mass,width,x,y,Radius)= Create_icebergs(lon_init,lat_init,Radius,R_earth\
			,x,y,ice_mask,h_ice,Convert_to_lat_lon,rho_ice,input_is_cartesian,element_type,scale_the_grid_to_lat_lon,lat_ref,adjust_lat_ref,\
			set_all_thicknesses_to_one,Interpolate_from_four_corners)

	#Define the positions of the iceberg bonds
	if Create_icebergs_bonds==True:
		(Number_of_bonds, first_berg_num,first_berg_ine,first_berg_jne,first_berg_lat,first_berg_lon, other_berg_num,other_berg_ine, other_berg_jne,other_berg_lat,other_berg_lon)=\
				Define_iceberg_bonds(Number_of_bergs,iceberg_num,lat,lon,dx_berg, dy_berg,Radius,break_some_bonds)

	if only_choose_one_berg==True:
		(Number_of_bergs,lon,lat,thickness,width,mass,iceberg_num)= Select_just_one_berg(lon,lat,thickness,width,mass,iceberg_num,chosen_berg_num)

	if save_restart_files==True:	
		#Creating iceberg restart file
		Create_iceberg_restart_file(Number_of_bergs, lon,lat,thickness,width,mass,mass_scaling,iceberg_num,Ice_geometry_source)
		
		#Creating bond restart file
		if Create_icebergs_bonds==True:
			Create_bond_restart_file(Number_of_bonds,first_berg_num,first_berg_ine,first_berg_jne,other_berg_ine,other_berg_jne,iceberg_num,other_berg_num,Ice_geometry_source)
			

	if regrid_icebergs_onto_grid==True:
		if scale_the_grid_to_lat_lon==True:
			print 'Regridding should be run with scale_the_grid_to_lat_lon off'
		else:
			if Switch_regridding_element_type==True:
				if element_type=='square':
					element_type='hexagon'
				else:
					element_type='square'
			regrid_iceberg_thickness(lat,lon,Number_of_bergs,thickness,mass,h_ice,x,y,rho_ice,element_type)


	# Plotting the positions and bonds of the newly formed formation
	if plotting_bergs_over_grid==True:
		plotting_iceberg_positions(lat,lon,Number_of_bergs,R_earth,Radius,IA_scaling,Convert_to_lat_lon,\
				plot_circles,h_ice,ice_mask,x,y,plot_ice_mask,plot_ice_thickness,thickness,plot_icebergs_positions)
		if (plot_bonds==True) and (Create_icebergs_bonds):
			plotting_iceberg_bonds(first_berg_lat,first_berg_lon,other_berg_lat,other_berg_lon,Number_of_bonds)


	#x0=-0.31580
	#y0=0.2
	#H=0.3
	#(Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)= Divide_hexagon_into_4_quadrants_old(x0,y0,H)
	#print x0,y0,H
	#print Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
	
	plt.show()
	print 'Script complete'



if __name__ == '__main__':
	main()
	#sys.exit(main())














