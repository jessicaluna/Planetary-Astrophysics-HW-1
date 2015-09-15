from pylab import*
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import csv

def rvcalculator(a,e,i,omega,m1,m2,T,to,t,distance):
    G= 6.673E-8
    M=(2*pi)*(t-to)/T
    rad_as=206264.806 # arcseconds in 1 radian
    #setting initial values for trancendental equation
    Eo=M+e*sin(M)
    g=Eo-e*sin(Eo)-M
    gprime=1-e*cos(Eo)
    Bigomega=0.0
    #solving trancendental equation
    for j in arange(1,20):
        E=Eo-g/gprime
        Eo=E
        g=Eo-e*sin(Eo)-M
        gprime=1-e*cos(Eo)
    
    
    E=Eo

    ##finding r,f
    r=a*(1.0-e*cos(E))

    f=2.0*arctan(sqrt((1+e)/(1-e))*tan(E/2.0))

    a1=m2*a/(m1+m2)
    a2=m1*a/(m1+m2)
    ##CALCULATE THE RELATIVE ASTROMETRY
    X_rel=r*(cos(Bigomega)*cos(omega+f)-sin(Bigomega)*sin(omega+f)*cos(i))
    Y_rel=r*(sin(Bigomega)*cos(omega+f)+cos(Bigomega)*sin(omega+f)*cos(i))
    Z_rel=r*sin(omega+f)+sin(i)

    ##finding projected separation of planet relative to star
    rho_rel=sqrt(X_rel**2+Y_rel**2)#/distance

    #calculating position angle of planet relative to star
    if (X_rel)<0 and (Y_rel>=0):
        theta_rel=arctan((X_rel/Y_rel))+pi/2.0
      
     
    elif (X_rel)<=0 and (Y_rel)<0:
        theta_rel=arctan((Y_rel/X_rel))+pi
        
   
    elif (X_rel>0) and (Y_rel<=0):
        theta_rel=arctan((X_rel/Y_rel))+ 3.0*pi/4.0
        

    elif (X_rel>=0) and (Y_rel>0):
        theta_rel=arctan((Y_rel/X_rel))
     
       

    ##CALCULATE THE ABSOLUTE ASTROMETRY
    ##FOR PRIMARY (1)
    X_abs_1=(a1/a)*-1.0*X_rel
    Y_abs_1=(a1/a)*-1.0*Y_rel
    Z_abs_1=(a1/a)*-1.0*Z_rel

    rho_abs_1=sqrt(X_abs_1**2+Y_abs_1**2)/distance

    if (X_abs_1)<0 and (Y_abs_1>=0):
            theta_abs_1=arctan((X_abs_1/Y_abs_1))+ pi/2.0
            
    elif (X_abs_1)<=0 and (Y_abs_1)<0:
            theta_abs_1=arctan((Y_abs_1/X_abs_1))+pi
           
    elif (X_abs_1>0) and (Y_abs_1<=0):
            theta_abs_1=arctan((X_abs_1/Y_abs_1))+ 3.0*pi/4.0
           
    elif (X_abs_1>=0) and (Y_abs_1>0):
            theta_abs_1=arctan((Y_abs_1/X_abs_1))
            
    ##calculating RA and DEC due to planetary companion

    RA_p=X_abs_1/distance   ##change in RA in radians
    Dec_p=Y_abs_1/distance  ##change in RA in radians

   



##Calculating RA and DEC due to proper motion
    

    ##FOR SECONDARY (2)
    X_abs_2=(a2/a)*-1.0*X_rel
    Y_abs_2=(a2/a)*-1.0*Y_rel
    Z_abs_2=(a2/a)*-1.0*Z_rel

    rho_abs_2=sqrt(X_abs_2**2+Y_abs_2**2)/distance
    if (X_abs_2)<0 and (Y_abs_2>=0):
            theta_abs_2=arctan((X_abs_2/Y_abs_2)) +pi/2.0      
    elif (X_abs_2)<=0 and (Y_abs_2)<0:
            theta_abs_2=arctan((Y_abs_2/X_abs_2))+pi
    elif (X_abs_2>0) and (Y_abs_2<=0):
            theta_abs_2=arctan((X_abs_2/Y_abs_2))+ 3.0*pi/4.0
    elif (X_abs_2>=0) and (Y_abs_1>0):
            theta_abs_2=arctan((Y_abs_2/X_abs_2))
    


    #calculate RVS of star and planet
    rvstar=(2*pi*(a1)*sin(i)/(T*sqrt(1-e**2)))*(cos(f+omega)+e*cos(omega))
    rvplanet=(2*pi*(a2)*sin(i)/(T*sqrt(1-e**2)))*(cos(f+omega)+e*cos(omega))
    
    return (rvstar/1.0e5,rvplanet/1.0e5,rho_rel,theta_rel,rho_abs_1,theta_abs_1,rho_abs_2,theta_abs_2,X_rel,Y_rel,Z_rel,X_abs_1,Y_abs_1,X_abs_2,Y_abs_2,RA_p,Dec_p)#,RA_parallax, DEC_parallax)

aa,e,i,omega,m1a,m2a,to,T,distance1=loadtxt('HD80606belements.txt',unpack=True, usecols=[0,1,2,3,4,5,6,7,8])


au_cm=1.496e13  #cm in 1au
m_jup=1.898e30 #Mjupiter in grams
m_sun=1.989e33 #Msun in grams
m1=m1a*m_sun
m2=m2a*m_jup


pc_cm=3.08e18   #cm in 1pc
distance=distance1*pc_cm  #distance in cm
a=aa*au_cm   #a in cms


i=radians(i)
omega=radians(omega)
t_sec=86400.0 #conversion from days to seconds
to=to*t_sec    #time of periapse
tstart=2457235.5*t_sec #JD of August 1st in seconds
tend=2457387.50*t_sec  #JD of December 31st in seconds
T=T*t_sec


time=linspace(tstart,tend,10000)
n=len(time)

##setting lists for outputs of function
rvs=[]
rvp=[]
X_rel=[]
Y_rel=[]
Z_rel=[]
X_abs_1=[]
Y_abs_1=[]
X_abs_2=[]
Y_abs_2=[]
rho_rel=[]
theta_rel=[]
rho_abs_1=[]
theta_abs_1=[]
rho_abs_2=[]
theta_abs_2=[]
RA_p=[]
DEC_p=[]


##running code for multiple times
for cc in range (0,n):
    rvs1,rvp1,rho_rel1,theta_rel1,rho_abs_1a,theta_abs_1a,rho_abs_2a,theta_abs_2a,X_rel1,Y_rel1,Z_rel1,X_abs_1a,Y_abs_1a,X_abs_2a,Y_abs_2a,RA_pa,DEC_pa=rvcalculator(a,e,i,omega,m1,m2,T,to,time[cc],distance)


    ##list of RVs of planet and star
    rvs.append(rvs1)
    rvp.append(rvp1)
    X_rel.append(X_rel1)
    Y_rel.append(Y_rel1)
    Z_rel.append(Z_rel1)
    X_abs_1.append(X_abs_1a)
    Y_abs_1.append(Y_abs_1a)
    X_abs_2.append(X_abs_2a)
    Y_abs_2.append(Y_abs_2a)
    ##relative astrometry of projected separation and PA
    rho_rel.append(rho_rel1)
    theta_rel.append(theta_rel1)
    
    ##absolute astrometry of projected separation and PA of star(1) and planet(2)
    rho_abs_1.append(rho_abs_1a)
    theta_abs_1.append(theta_abs_1a)
    rho_abs_2.append(rho_abs_2a)
    theta_abs_2.append(theta_abs_2a)

    ##RA and DEC from just planet
    RA_p.append(RA_pa)
    DEC_p.append(DEC_pa)




    
JD_times=Time(time/t_sec, format='jd', scale='utc')

## Question 2 finding the time of min and max of the RV
extreme_min=min(rvs)
index_extreme_min=rvs.index(extreme_min)
time_extreme_min= time[index_extreme_min]
extreme_min_date=JD_times[index_extreme_min].iso
print('Date of Extreme RV min=', extreme_min_date)

extreme_max=max(rvs)
index_extreme_max=rvs.index(extreme_max)
time_extreme_max= time[index_extreme_max]
extreme_max_date=JD_times[index_extreme_max].iso
print('Date of Extreme RV max=', extreme_max_date)

##Question 3 finding the time of transit
closedist=min(rho_rel)
indexclose=rho_rel.index(closedist)
timeclose=time[indexclose]/t_sec

Z_rela=np.asarray(Z_rel)
rho_rela=np.asarray(rho_rel)
index_posz=np.where(Z_rela>0)
index_posza=index_posz[0]
element_minsep=np.argmin(rho_rela[index_posza])
true_index=index_posza[element_minsep]


transit_time=JD_times[true_index].iso
print('Date of transit=',transit_time)



fname='ra.txt'
fname2='dec.txt'
savetxt(fname, RA_p, fmt='%.18e', delimiter=' ', newline='\n', header='RA radians')
savetxt(fname2, DEC_p, fmt='%.18e', delimiter=' ', newline='\n', header='DEC radians')

figure(1)
plot((time)/t_sec - 2457000,rvs,c='r')
annotate('Oct 25',xy=(time_extreme_min/t_sec-2457000,extreme_min), xytext=(290, -.3),arrowprops=dict(facecolor='black', shrink=0.05))
annotate('Oct 26',xy=(time_extreme_max/t_sec-2457000,extreme_max), xytext=(350, 0.75),arrowprops=dict(facecolor='black', shrink=0.05))
xlabel('2475k + JD time(days) ')
ylabel('Radial Velocity (km/s)')
title('RV vs Time (Aug to Dec)')



figure(2)
plot(time/t_sec-2457000,rho_rel,c='k')
title('Relative Projected Separation Vs time ')
xlabel('2475k + JD time(days) ')
ylabel('Relative Projected Separation (cm)')
annotate('Nov 11, 2:17 am',xy=(time[true_index]/t_sec-2457000,0.0), xytext=(350, 1.0*10**(12)),arrowprops=dict(facecolor='black', shrink=0.05))
annotate('Transit Date:',xy=(time[true_index]/t_sec-2457000,0.0), xytext=(350, 1.3*10**(12)))

show()


