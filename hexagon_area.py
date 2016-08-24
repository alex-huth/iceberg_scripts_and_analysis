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



#######################

def Divide_hexagon_into_4_quadrants_old(x0,y0,H):
	S=(2/sqrt(3))*H
	Area_hex=3*(np.sqrt(3)/2)*(S**2)  #Area of the hexagon (should be equal to Area/grid_area, since it is not dim)  - check this
	#print 'y0,H',y0, H	
	#Defining boundaries of hexagon, and if statements to see which side of the boundary you are on
	W1= False ; W2=False ;W3=False ;W4=False ;W5=False ;W6=False ;H1=False ;V1=False ;V2=False
	W1=(-y0<-sqrt(3)*(-x0) + ((sqrt(3)*(S))));#upper right
	W2=(-y0<(H));#Top
	W3=(-y0<sqrt(3)*(-x0) + ((sqrt(3)*(S))));#Upper left
	W4=(-y0<sqrt(3)*(-x0) + (-(sqrt(3)*(S))));#Lower right
	W5=(-y0>-H);#Bottom
	W6=(-y0<-sqrt(3)*(-x0) + (-(sqrt(3)*(S))));#Lower left
		

	T1=(-x0<S) #Right
	T2=(-y0<(H));#Top
	T3= (-x0>-S) #Left
	T4=(-y0>-H);#Bottom

	H1=(y0<0);
	V1=(x0<-(S/2));
	V2=(x0<(S/2));
	  
	#Deciding if the origin is within the hexagon
	#print W1 , W2 , W3 , W4 , W5 , W6
	#if In_hex:
	#	print In_hex
		
	#Calculating the area of the top and bottom half of the hexagon, 2 Cases for the majority above and below the y0=0 line
	#(and two more for the hexagon totally above and below the y0=0 line)
	if abs(y0)<H:
		Trapesium=((sqrt(3)*H)-(abs(y0)/sqrt(3)))*(H-abs(y0));
		if y0>=0:
			Area_Lower=Trapesium;
			Area_Upper=Area_hex-Trapesium;
		else:
			Area_Upper=Trapesium;
			Area_Lower=Area_hex-Trapesium;
	else:
		if y0>=0:
			Area_Lower=0.;
			Area_Upper=Area_hex;
		else: 
			Area_Lower=Area_hex;
			Area_Upper=0.;

	#Calcularing Left and Right area of the hexagon, about the x0=0 line, 3 cases:
	#(and two more for when the hexagon is totally to the left or right of the x0=0 line)
	if abs(x0)<S:
		if abs(x0)<S/2:
			Rectangle=(abs(x0)*2*H);
			Big_side =(Area_hex/2) +Rectangle;
			Small_side=Area_hex-Big_side;
		else:
			Triangle=(sqrt(3))*((S-abs(x0))**2);
			Small_side=Triangle;
			Big_side=Area_hex-Small_side;
		if x0>=0.:
			Area_right=Big_side;
			Area_left=Small_side;
		else:
			Area_right=Small_side;
			Area_left=Big_side;
	else:
		if x0>=0.:
			Area_right=Area_hex;
			Area_left=0.;
		else:
			Area_right=0.;
			Area_left=Area_hex;

	In_hex= W1 & W2 & W3 & W4 & W5 & W6;
	In_hex_box=T1 & T2 & T3 & T4
	Area_Q1=0.; Area_Q2=0. ; Area_Q3=0.; Area_Q4=0.;
	Sector=0
	#if In_hex==False: #Then the hexagon is completely contained in the middle cell
	if In_hex_box==False: #Then the hexagon is completely contained in the middle cell
		Sector=-1
		#mass_on_ocean[i,j,5]=mass_on_ocean[i,j,5]+Mass
		if min(Area_Upper,Area_Lower)==0.:
			Sector=-2
			if Area_Upper==0.:
				Area_Q3=Area_left;
				Area_Q4=Area_right;
			if Area_Lower==0.:
				Area_Q1=Area_right;
				Area_Q2=Area_left;

		elif min(Area_right,Area_left)==0.:
			Sector=-3
			if Area_right==0.:
				Area_Q2=Area_Upper;
				Area_Q3=Area_Lower;
		
			if Area_left==0.:
				Area_Q1=Area_Upper;
				Area_Q4=Area_Lower;
		

		#yCxC=1.
		#print 'out of hex'

	else:
		#Determine which sector within the hexagon you are in. (sectors 1 to 6 go counter clockwise starting with top right)
		if (H1==True):  #Bottom half
			if V1:
				#if W1==False:
				if ((y0+(sqrt(3)*(x0+S)))<=0.):
					Sector=1;
				else:
					Sector=2;
			elif (V1==False) & (V2==True):
				Sector=3;
			else:  
				#if (W3==True):
				if ((y0-(sqrt(3)*(x0-S)))>=0.):
					Sector=4;
				else:
					Sector=5;
		else:  #Bottom half
			if V1:
				#if W6==False:
				if ((y0 -(sqrt(3)*(x0+S)))>=0.):
					Sector=10;
				else:
					Sector=9;
			elif (V1==False) & (V2==True):
				Sector=8;
			else:  
				#if (W4==True):
				if ((y0+(sqrt(3)*(x0-S)))<=0.):
					Sector=7;
				else:
					Sector=6;

		#print Sector

		#If the hexagon is in Sector 1,3,4 or 6, then the intersetion of the hexagon and the corresponding sector forms a baby triangle
		#If the hexagon is in Sector 2,5 then the intersetion of the hexagon and the corresponding sector forms a baby trapesoid
		if Sector==2 or Sector==4  or Sector==7 or Sector==9:
			Baby_triangle=(1/(2*sqrt(3)))*((-abs(y0)+(sqrt(3)*(S-abs(x0))))**2);
		else:
			#Baby_trap= (H-abs(y0)) * ((-H-abs(y0)+(2*sqrt(3)*(S-abs(x0))))/(2*sqrt(3)));  
			Baby_trap=(H-abs(y0))*((S-abs(x0) - ((H+abs(y0))/(2*sqrt(3))))) ;

		#Finally, we assign the correct areas in each quadrant (Q1,Q2,Q3,Q4), depending on which sector you are in.
		C1=0.;C2=0.;C3=0.;C4=0.;

		#Corner cases
		if Sector==2:
			Area_Q1=Baby_triangle;
			Area_Q2=Area_Upper-Area_Q1
			Area_Q3=Area_left-Area_Q2
			Area_Q4=Area_right-Area_Q1

		if Sector==4:
			Area_Q2=Baby_triangle;
			Area_Q1=Area_Upper-Area_Q2
			Area_Q3=Area_left-Area_Q2
			Area_Q4=Area_right-Area_Q1

		if Sector==7:
			Area_Q3=Baby_triangle;
			Area_Q2=Area_left-Area_Q3
			Area_Q1=Area_Upper-Area_Q2
			Area_Q4=Area_right-Area_Q1

		if Sector==9:
			Area_Q4=Baby_triangle;
			Area_Q1=Area_right-Area_Q4
			Area_Q2=Area_Upper-Area_Q1
			Area_Q3=Area_left-Area_Q2
	
		#Center cases
		if Sector==3:
			if x0<=0.:
				Area_Q1=Baby_trap;
				Area_Q2=Area_Upper-Area_Q1;
				Area_Q3=Area_left-Area_Q2;
				Area_Q4=Area_right-Area_Q1;
			else:
				Area_Q2=Baby_trap;
				Area_Q1=Area_Upper-Area_Q2;
				Area_Q3=Area_left-Area_Q2;
				Area_Q4=Area_right-Area_Q1;
			
		if Sector==8:
			if x0<=0.:
				Area_Q4=Baby_trap;
				Area_Q3=Area_Lower-Area_Q4;
				Area_Q1=Area_right-Area_Q4;
				Area_Q2=Area_Upper-Area_Q1;
			else:
				Area_Q3=Baby_trap;
				Area_Q4=Area_Lower-Area_Q3;
				Area_Q1=Area_right-Area_Q4;
				Area_Q2=Area_Upper-Area_Q1;

		#Outside triangle cases:
		if Sector==1:
			Area_Q1=0.;
			Area_Q2=Area_Upper;
			Area_Q4=Area_right;
			Area_Q3=Area_left-Area_Q2;
		
		if Sector==5:
			Area_Q2=0.;
			Area_Q1=Area_Upper;
			Area_Q3=Area_left;
			Area_Q4=Area_Lower-Area_Q3;

		if Sector==6:
			Area_Q3=0.;
			Area_Q2=Area_left;
			Area_Q4=Area_Lower;
			Area_Q1=Area_right-Area_Q4;
		
		if Sector==10:
			Area_Q4=0.;
			Area_Q3=Area_Lower;
			Area_Q1=Area_right;
			Area_Q2=Area_left-Area_Q3;

		#print x0,y0,Sector
	return [Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4]



def Hexagon_into_quadrants_using_triangles(x0,y0,H,theta):
		 
	#Length of side of Hexagon
	S=(2/sqrt(3))*H;  
	 
	#Finding positions of corners
	C1x=S  +x0        ; C1y=0.+y0;  #Corner 1 (right)
	C2x=H/sqrt(3) +x0 ; C2y=H+y0;  #Corner 2 (top right)
	C3x=-H/sqrt(3)+x0 ; C3y=H+y0;  #Corner 3 (top left)
	C4x=-S     +x0    ; C4y=0.+y0;  #Corner 4 (left)
	C5x=-H/sqrt(3) +x0; C5y=-H+y0; #Corner 5 (top left)
	C6x=H/sqrt(3) +x0 ; C6y=-H+y0; #Corner 3 (top left)
	 
	 
	#Area of Hexagon is the sum of the triangles
	[T12_Area,T12_Q1,T12_Q2,T12_Q3,T12_Q4]=Triangle_divided_into_four_quadrants(x0,y0,C1x,C1y,C2x,C2y); #T012
	[T23_Area,T23_Q1,T23_Q2,T23_Q3,T23_Q4]=Triangle_divided_into_four_quadrants(x0,y0,C2x,C2y,C3x,C3y); #T023
	[T34_Area,T34_Q1,T34_Q2,T34_Q3,T34_Q4]=Triangle_divided_into_four_quadrants(x0,y0,C3x,C3y,C4x,C4y); #T034
	[T45_Area,T45_Q1,T45_Q2,T45_Q3,T45_Q4]=Triangle_divided_into_four_quadrants(x0,y0,C4x,C4y,C5x,C5y); #T023
	[T56_Area,T56_Q1,T56_Q2,T56_Q3,T56_Q4]=Triangle_divided_into_four_quadrants(x0,y0,C5x,C5y,C6x,C6y); #T023
	[T61_Area,T61_Q1,T61_Q2,T61_Q3,T61_Q4]=Triangle_divided_into_four_quadrants(x0,y0,C6x,C6y,C1x,C1y); #T023

	#Summing up
	Area_hex=T12_Area+T23_Area+T34_Area+T45_Area+T56_Area+T61_Area;
	Area_Q1=T12_Q1+T23_Q1+T34_Q1+T45_Q1+T56_Q1+T61_Q1;
	Area_Q2=T12_Q2+T23_Q2+T34_Q2+T45_Q2+T56_Q2+T61_Q2;
	Area_Q3=T12_Q3+T23_Q3+T34_Q3+T45_Q3+T56_Q3+T61_Q3;
	Area_Q4=T12_Q4+T23_Q4+T34_Q4+T45_Q4+T56_Q4+T61_Q4;
	#Area_Q4=Area_hex-(Area_Q1+Area_Q2+Area_Q3)

	Area_Q1=max(Area_Q1,0.);
	Area_Q2=max(Area_Q2,0.);
	Area_Q3=max(Area_Q3,0.);
	Area_Q4=max(Area_Q4,0.);


	Error=Area_hex-(Area_Q1+Area_Q2+Area_Q3+Area_Q4)
	if (abs(Error)>0.01):
		print  'diamonds, hex error, H,x0,y0, Error', H, x0 , y0, Error
        	print 'diamonds, hex error, Areas',Area_hex, (Area_Q1+Area_Q2 + Area_Q3+Area_Q4), Area_Q1,  Area_Q2 , Area_Q3,  Area_Q4
	      	print 'Hexagon error is large!!'



	#Adjust Areas so that the error is zero by subtracting the error from the largest sector.
	if  (((Area_Q1>=Area_Q2) and (Area_Q1>=Area_Q3)) and (Area_Q1>=Area_Q4)):
		print 'fix1',Error
		Area_Q1=Area_Q1+Error
	elif  (((Area_Q2>=Area_Q1) and (Area_Q2>=Area_Q3)) and (Area_Q2>=Area_Q4)):
		print 'fix2', Error
		Area_Q2=Area_Q2+Error
	elif  (((Area_Q3>=Area_Q1) and (Area_Q3>=Area_Q2)) and (Area_Q3>=Area_Q4)):
		print 'fix3',Error
		Area_Q3=Area_Q3+Error
   	elif  (((Area_Q4>=Area_Q1) and (Area_Q4>=Area_Q2)) and (Area_Q4>=Area_Q3)):
		print 'fix4',Error
		Area_Q4=Area_Q4+Error
	else:
		print 'There is some thing wrong with this hexagon. Something very wrong'

	#Error=Area_hex-(Area_Q1+Area_Q2+Area_Q3+Area_Q4)
	Error=Area_hex-(Area_Q1+Area_Q2+Area_Q3+Area_Q4)
	if ((abs(Error)>0.01)):
		print 'The hexagon error is still too large!!!', Error



	return [Area_hex ,Area_Q1, Area_Q2, Area_Q3, Area_Q4]

def Triangle_divided_into_four_quadrants(Ax,Ay,Bx,By,Cx,Cy):
	Area_triangle=Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy);

	#Round of numbers before proceeding further.
        #print Cy
	#sig_fig=12; #Significan figures
	#Ax=roundoff(Ax,sig_fig)
	#Ay=roundoff(Ay,sig_fig)
	#Bx=roundoff(Bx,sig_fig)
	#By=roundoff(By,sig_fig)
	#Cx=roundoff(Cx,sig_fig)
	#Cy=roundoff(Cy,sig_fig)
        #print Cy
	
	#Calculating area across axes
	[Area_Upper ,Area_Lower]=divding_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,'x');
	[Area_Right ,Area_Left]=divding_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,'y');
	 
	#Decide if the origin is in the triangle
	if point_in_triangle(Ax,Ay,Bx,By,Cx,Cy,0.,0.):  #Then you have to divide area 4 ways.
		#Find a line in the triangle that cuts both axes in/on the trianlge
	    	[px, py]=intercept_of_a_line(Ax,Ay,Bx,By,'x'); #x_intercept
	    	[qx, qy]=intercept_of_a_line(Ax,Ay,Bx,By,'y'); #y_intercept
	    	if (point_in_interval(Ax,Ay,Bx,By,px,py) & point_in_interval(Ax,Ay,Bx,By,qx,qy))==False:
			[px, py]=intercept_of_a_line(Ax,Ay,Cx,Cy,'x'); #x_intercept
			[qx, qy]=intercept_of_a_line(Ax,Ay,Cx,Cy,'y'); #y_intercept
			if (point_in_interval(Ax,Ay,Cx,Cy,px,py) & point_in_interval(Ax,Ay,Cx,Cy,qx,qy))==False:
				[px, py]=intercept_of_a_line(Bx,By,Cx,Cy,'x'); #x_intercept
				[qx, qy]=intercept_of_a_line(Bx,By,Cx,Cy,'y'); #y_intercept
				if (point_in_interval(Bx,By,Cx,Cy,px,py) & point_in_interval(Bx,By,Cx,Cy,qx,qy))==False:
					'Houston, we have a problem'
					#halt

		#Assigning quadrants. Key_quadrant is the quadrant with the baby triangle in it.
		Area_key_quadrant=Area_of_triangle(px,py,qx,qy,0.,0.);
		if px>=0. and qy>=0.: #First quadrant
			Key_quadrant=1;
		elif px<0. and qy>=0.:  #Second quadrant
			Key_quadrant=2;
		elif px<0. and qy<0.:
			Key_quadrant=3;  #Third quadrant
		else:
			Key_quadrant=4;  #Forth quadrant
	    
	else:  #Then at least one quadrant is empty, and this can be used to find the areas in the other quadrant.  Assigning quadrants. Key_quadrant is the empty quadrant.
		print 'Mother...'
		print 'Ax, Ay',Ax,Ay
		print 'Bx, By',Bx,By
		print 'Cx, Cy',Cx,Cy
		Area_key_quadrant=0;
		if (((Ax>0. and Ay>0.) or (Bx>0. and By>0.) or (Cx>0. and Cy>0.))==False)  and (Area_Upper+Area_Right<=Area_triangle):
			#No points land in this quadrant and triangle does not cross the quadrant
			Key_quadrant=1;
		elif  (((Ax<0. and Ay>=0) or (Bx<0. and By>=0.) or (Cx<0. and Cy>=0.))==False)  and (Area_Upper+Area_Left<=Area_triangle):
			Key_quadrant=2;
		elif (((Ax<0. and Ay<0.) or (Bx<0. and By<0.) or (Cx<0. and Cy<0.))==False) & (Area_Lower+Area_Left<=Area_triangle):
			Key_quadrant=3;
		else:
			Key_quadrant=4;
	
		print Key_quadrant
	#Assign values to quadrants
	if Key_quadrant==1:
		Area_Q1=Area_key_quadrant;
		Area_Q2=Area_Upper-Area_Q1;
		Area_Q4=Area_Right-Area_Q1;
		#Area_Q3=Area_Left-Area_Q2;
		Area_Q3=Area_triangle-Area_Q1-Area_Q2-Area_Q4;
	elif Key_quadrant==2:
		Area_Q2=Area_key_quadrant;
		Area_Q1=Area_Upper-Area_Q2;
		Area_Q4=Area_Right-Area_Q1;
		#Area_Q3=Area_Left-Area_Q2;
		Area_Q3=Area_triangle-Area_Q1-Area_Q2-Area_Q4;
	elif Key_quadrant==3:
		Area_Q3=Area_key_quadrant;
		Area_Q2=Area_Left-Area_Q3;
		Area_Q1=Area_Upper-Area_Q2;
		#Area_Q4=Area_Right-Area_Q1;
		Area_Q4=Area_triangle-Area_Q1-Area_Q2-Area_Q3;
	elif Key_quadrant==4:
		Area_Q4=Area_key_quadrant;
		Area_Q1=Area_Right-Area_Q4;
		Area_Q2=Area_Upper-Area_Q1;
		#Area_Q3=Area_Left-Area_Q2;
		Area_Q3=Area_triangle-Area_Q1-Area_Q2-Area_Q4;
	else:
	    print 'Help, I need somebody, help!'
	    halt
	
	Area_Q1=max(Area_Q1,0.); 
	Area_Q2=max(Area_Q2,0.);
	Area_Q3=max(Area_Q3,0.);
	Area_Q4=max(Area_Q4,0.); 

	Error=abs(Area_Q1+Area_Q2+Area_Q3+Area_Q4-Area_triangle)
	if Error>0.01:
		print 'The triangles are not accurate enough. This is a problem!'
		print 'Triangle corners: ' ,Ax,Ay,Bx,By,Cx,Cy
		print 'Error',Error
		return

	return [Area_triangle, Area_Q1, Area_Q2 ,Area_Q3 ,Area_Q4]

def divding_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,axes1):
	if axes1=='x':  #Use the y-coordinates for if statements to see which side of the line you are on
		A0=Ay; B0=By;  C0=Cy;
	if axes1=='y':  #Use the y-coordinates for if statements to see which side of the line you are on
		A0=Ax; B0=Bx;  C0=Cx;
 
	A_triangle=Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy);
	if B0*C0>0.: #then B and C are on the same side  (and non-zero)
		if A0*B0>=0.: #then all three on the the same side (if it equals zero, then A0=0 and the otehrs are not)
			if (A0>0.)  or  (A0==0. and  B0>0.):
				Area_positive= A_triangle;
				Area_negative= 0.;
			else:
				Area_positive= 0.;
				Area_negative= A_triangle;
		else:  #A is on the opposite side to B and C
			[Area_positive, Area_negative]=Area_of_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,axes1);

	elif B0*C0<0.: #then B and C are on the opposite sides
		if A0*B0>=0.:  #then C is all alone
			[Area_positive, Area_negative]=Area_of_triangle_across_axes(Cx,Cy,Bx,By,Ax,Ay,axes1);
		else: #then B is all alone
			[Area_positive, Area_negative]=Area_of_triangle_across_axes(Bx,By,Cx,Cy,Ax,Ay,axes1);

	else:  #This is the case when either B or C is equal to zero (or both), A0 could be zero too.
		if (A0==0. and B0==0. and C0==0.):
			Area_positive= 0.;
			Area_negative= 0.;
		elif (A0*B0<0.)  or  (A0*C0<0.):    #A, B are on opposite sides, and C is zero.  OR  A, C are on opposite sides, and B is zero.
			[Area_positive, Area_negative]=Area_of_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,axes1);
		elif (A0*B0>0.) or (A0*C0>0.) or (abs(A0)>0. and (B0==0.) and (C0==0.)):
		        if (A0>0.):
		        	Area_positive= A_triangle;
				Area_negative= 0.;
        		else:
				Area_positive= 0.;
				Area_negative= A_triangle;
    
    		elif A0==0.:  #(Also, one of B,C is zero too)
			if B0>0. or C0>0.:
				Area_positive= A_triangle;
				Area_negative= 0.;
			elif B0<0. or C0<0.:
				Area_positive= 0.;
				Area_negative= A_triangle;
		        else:
				print 'You should not get here1'
				halt
    		else:
			print 'You should not get here2'
			halt
	
	return [Area_positive, Area_negative]

def check_if_point_is_on_the_line(Ax,Ay,Bx,By,qx,qy):
	tol=0.00000000000000;
	dxc = qx - Ax;
	dyc = qy - Ay;
 
	dxl = Bx - Ax;
	dyl = By - Ay;
 
	cross = dxc * dyl - dyc * dxl;
	if abs(cross)<=tol:
		point_is_on_line=True
	else:
		point_is_on_line=False

	return point_is_on_line

def intercept_of_a_line(Ax,Ay,Bx,By,axes1):
	No_intercept_val=100000000000.; #Huge value used to make sure that the intercept is outside the triange in the parralel case. 
	#No_intercept_val=np.NaN;
	 
	if axes1=='x': #x intercept
		if (Ay==By)==False:
			x0=Ax -(((Ax-Bx)/(Ay-By))*Ay);
			y0=0.; 
		else:
			x0=No_intercept_val;
			y0=No_intercept_val;
	
	if axes1=='y': #y intercept
		if (Ax==Bx)==False:
			x0=0.; 
			y0=-(((Ay-By)/(Ax-Bx))*Ax)+Ay;  
		else:
			x0=No_intercept_val;
			y0=No_intercept_val;
	
	return [x0, y0]


def Area_of_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,axis1):#In this fuction, A is the point on one side of the axis, and B,C are on the opposite sides
	A_triangle2=Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy);
	[pABx, pABy]=intercept_of_a_line(Ax,Ay,Bx,By,axis1);
	[pACx, pACy]=intercept_of_a_line(Ax,Ay,Cx,Cy,axis1);
	if axis1=='x':
		A0=Ay; #Value used for if statements (deciding up/down vs left/right)
	if axis1=='y':
		A0=Ax; #Value used for if statements (deciding up/down vs left/right)
	
	A_half_triangle=Area_of_triangle(Ax,Ay,pABx,pABy,pACx,pACy);
	if (A0>=0.):
		Area_positive= A_half_triangle;
		Area_negative= A_triangle2-A_half_triangle;

	else:
		Area_positive= A_triangle2-A_half_triangle;
		Area_negative= A_half_triangle;

	return [Area_positive, Area_negative]

def Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy):
	Area    =   abs(    0.5*((Ax*(By-Cy))+(Bx*(Cy-Ay))+(Cx*(Ay-By))) );
	return Area

def point_in_interval(Ax,Ay,Bx,By,px,py):
	point_is_in_interval=False
	if ((px < max(Ax,Bx)) and (px > min(Ax,Bx))):
		if ((py < max(Ay,By)) and (py > min(Ay,By))):
			point_is_in_interval=True
	return point_is_in_interval



def point_in_triangle_old(Ax,Ay,Bx,By,Cx,Cy,qx,qy):
	if (Ax==qx and Ay==qy) or (Bx==qx and By==qy) or (Cx==qx and Cy==qy):  #Exclude the pathelogical case
	        	point_is_in_triangle = 0.;
			print 'Mark1'
	else:
		if (check_if_point_is_on_the_line(Ax,Ay,Bx,By,qx,qy) or (check_if_point_is_on_the_line(Ax,Ay,Cx,Cy,qx,qy)) or (check_if_point_is_on_the_line(Bx,By,Cx,Cy,qx,qy))):
			point_is_in_triangle = 0;
			print 'Mark2'
		else:
			print 'Mark3'

			#Compute vectors
			v0x=Cx-Ax;
			v1x=Bx-Ax;
			v2x=qx-Ax;
			v0y=Cy-Ay;
			v1y=By-Ay;
			v2y=qy-Ay;

			#%Compute dot products
			dot00 = (v0x*v0x)+(v0y*v0y);
			dot01 = (v0x*v1x)+(v0y*v1y);
			dot02 = (v0x*v2x)+(v0y*v2y);
			dot11 = (v1x*v1x)+(v1y*v1y);
			dot12 = (v1x*v2x)+(v1y*v2y);

			#Compute barycentric coordinates
			invDenom= 1 / ((dot00 * dot11) - (dot01*dot01));
			u=((dot11*dot02)-(dot01*dot12))*invDenom;
			v=((dot00*dot12)-(dot01*dot02))*invDenom;
			print u+v-1

			point_is_in_triangle = (((u)>=0) & ((v)>=0) & ((u+v)<(1)));
			print point_is_in_triangle 
	return point_is_in_triangle



def point_in_triangle(Ax,Ay,Bx,By,Cx,Cy,qx,qy):
	point_is_in_triangle = 0;
	if (Ax==qx and Ay==qy) or (Bx==qx and By==qy) or (Cx==qx and Cy==qy):  #Exclude the pathelogical case
	        	point_is_in_triangle = 0.;
			print 'Mark1'
	else:
		if (check_if_point_is_on_the_line(Ax,Ay,Bx,By,qx,qy) or (check_if_point_is_on_the_line(Ax,Ay,Cx,Cy,qx,qy)) or (check_if_point_is_on_the_line(Bx,By,Cx,Cy,qx,qy))):
			point_is_in_triangle = 0;
			print 'Mark2'
		else:
			print 'Mark3'
			l0=(qx-Ax)*(By-Ay)-(qy-Ay)*(Bx-Ax)
			l1=(qx-Bx)*(Cy-By)-(qy-By)*(Cx-Bx)
			l2=(qx-Cx)*(Ay-Cy)-(qy-Cy)*(Ax-Cx)

			p0=np.sign( l0);
			p1=np.sign( l1);
			p2=np.sign( l2);
			if (l0 == 0.):
				p0=0.
			if (l1==0.):
				p1=0.
			if (l2==0.):
				p2=0.

			if ( (abs(p0)+abs(p2))+(abs(p1)) == abs((p0+p2)+(p1)) ):
				point_is_in_triangle = 1;
	print point_is_in_triangle	
	return point_is_in_triangle


def roundoff(x,sig_fig):
	x_roundoff=round(x*(10**(sig_fig)))/(10**(sig_fig))
	#x_roundoff=(FLOAT (INT(x * (10.**sig_fig) + 0.5)) / (10.**sig_fig))
	return x_roundoff



def square_spreading_calculation(x,y,L):
	xL=min(0.5, max(0., 0.5-(x/L)))
        xR=min(0.5, max(0., (x/L)+(0.5-(1/L) )))
        xC=max(0., 1.-(xL+xR))
        yD=min(0.5, max(0., 0.5-(y/L)))
        yU=min(0.5, max(0., (y/L)+(0.5-(1/L) )))
        yC=max(0., 1.-(yD+yU))
	print 'x0=',x, ' y0=',y,' L=', L
	print 'xL,xC,xR', xL, xC, xR
	print 'yU,yC,yD', yU, yC, yD


####################################################################################################################################################
##########################################################  Main Program   #########################################################################
####################################################################################################################################################

def main():

	#test_case='hexagon'
	test_case='triangle'
	#test_case='square'



	if test_case=='hexagon':
		H=0.4
		x0=0.4
		y0=0.4
		(Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)= Divide_hexagon_into_4_quadrants_old(x0,y0,H)
		print x0,y0,H
		print 'First version', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
		(Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)= Hexagon_into_quadrants_using_triangles(x0,y0,H,theta=0)
		print 'Triangle version',Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
		
		print 'Analysis'	
		print 'sum of quadrants: ', (Area_Q1+Area_Q2+Area_Q3+Area_Q4)
		print 'Scaled areas: '
		print (Area_Q1+Area_Q2+Area_Q3+Area_Q4)/Area_hex, Area_Q1/Area_hex, Area_Q2/Area_hex, Area_Q3/Area_hex, Area_Q4/Area_hex
		

	elif test_case=='square':
	 
		x0=0.6999999999999
		y0=0.3
		L=0.6
		square_spreading_calculation(x0,y0,L)

	elif test_case=='triangle':
		#x0=4.537920628519032E-002
		#y0=-0.300000000000000
		#C2x=0.218584287042078
		#C2y=2.220446049250313E-016
		#C3x=-0.127825874471697
		#C3y=2.220446049250313E-016

		Ax=-9.591318871326848E-002
		Ay=-0.399999999999999

		Bx=0.135026918962581
		By= 5.551115123125783E-017

		Cx=-0.326853296389118
		Cy=5.551115123125783E-017

		[T23_Area,T23_Q1,T23_Q2,T23_Q3,T23_Q4]=Triangle_divided_into_four_quadrants(Ax,Ay,Bx,By,Cx,Cy); #T023

		print 'Triangle 2: ',Ax,Ay,Bx,By,Cx,Cy
		print 'Q1, Q2,Q3,Q4',T23_Q1,T23_Q2,T23_Q3,T23_Q4
		print 'Full triangle area',T23_Area
		print 'Error',(T23_Q1+T23_Q2+T23_Q3+T23_Q4-T23_Area)
	print 'Script complete'


if __name__ == '__main__':
	main()
	#sys.exit(main())














