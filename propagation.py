from pylab import*
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import random

############## CALCULATE THE RA AND DEC DUE TO PLANET


def ra_dec(a,e,i,omega,m1,m2,P,to,t,distance):
    

    M=(2*pi)*(t-to)/P

    #setting initial values for trancendental equation
    Eo=M+e*sin(M)
    g=Eo-e*sin(Eo)-M
    gprime=1-e*cos(Eo)
    Bigomega=0.0
        #solving trancendental equation
    for jj in arange(1,20):
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

    ##CALCULATE THE ABSOLUTE ASTROMETRY
        ##FOR PRIMARY (1)
    X_abs_1=(a1/a)*-1.0*X_rel
    Y_abs_1=(a1/a)*-1.0*Y_rel
    Z_abs_1=(a1/a)*-1.0*Z_rel

    delta_RA_p=X_abs_1/distance   ##change in RA in radians
    delta_DEC_p=Y_abs_1/distance  ##change in DEC in radians

    return(delta_RA_p,delta_DEC_p)

aa,e,i,omega,m1a,m2a,to,P,distance1=loadtxt('HD80606belements.txt',unpack=True, usecols=[0,1,2,3,4,5,6,7,8])
au_cm=1.496e13  #cm in 1au
m_jup=1.898e30 #Mjupiter in grams
m_sun=1.989e33 #Msun in grams
t_sec=86400.0 #conversion from days to seconds

m1=m1a*m_sun
m2=m2a*m_jup


pc_cm=3.08e18   #cm in 1pc
distance=distance1*pc_cm  #distance in cm
a=aa*au_cm   #a in cms

P=P*t_sec
i=radians(i)
omega=radians(omega)
to=to*t_sec    #time of periapse
tstart=2457235.5 #august 1st 2015 in JD
tend=2459061.5  #july 31 2015 2020 in JD



n=tend-tstart

t=linspace(tstart, tend,n*10.0)



ts=t*t_sec

delta_RA_p,delta_DEC_p=ra_dec(a,e,i,omega,m1,m2,P,to,ts,distance)


##################################################################

RA=140.65635797 ##RA at 1991.25
DEC=50.60370469  #DEC at 1991.25

deg_mas=3600000.0  #number of milliarc seconds in 1 degree

d_yr=365.25



delta_RA_p=degrees(delta_RA_p)
delta_DEC_p=degrees(delta_DEC_p)


sigma= 15 ##mas

#### Propagating the RA and DEC due to proper motion from when Hipparcos observed

PM_ra= 46.98/(d_yr*deg_mas) ##degree/day
PM_dec=6.92/(d_yr*deg_mas) ##degree/day

# CAlculate the RA and DEC of August 1st 2015



RA_start= RA + PM_ra*(tend-tstart)
DEC_start=DEC +PM_dec*(tend-tstart)



RA_p=RA_start+delta_RA_p
DEC_p=DEC_start+delta_DEC_p


indices=[]
for rr in range(0,100):
    index=int(random.random() * len(t))
    trand=t[int(random.random() * len(t))]
    indices.append(index)

index_rand=np.asarray(indices)
#print(index_rand)
#print(t[index_rand])
trand=t[index_rand]

RA_p_rand=RA_p[index_rand]
DEC_p_rand=DEC_p[index_rand]






## calculating RA and DEC from August 1st 2015-2020 due to proper motion
    
RA_pm=RA_start+PM_ra*t+RA_p
RA_pm_rand=RA_pm[index_rand]

DEC_pm=DEC_start+PM_dec*t+RA_p
DEC_pm_rand=DEC_pm[index_rand]




##calculating RA and DEC due to parallax
def parallax(RA,DEC,t):
    t_J2000=2451544.5  #JD of J2000
    parallax=17.13   #mas
    ra_mas=206264806.247 #mas in 1 radian
    parallax=parallax/ra_mas #converting parallax to radians
    alpha=radians(RA)     ## RA in radians
    delta=radians(DEC)      ## DEC in radians

    #transforming from RA and DEC to eciptic coordinates
    am_deg=60.0 #arc minute in 1 degree
    as_deg= 3600.0 #arc seconds in 1 degree
    T= (t -t_J2000)/(d_yr*100)
    epsilon=(23.0+26.0/am_deg+ 21.45/as_deg)-(46.815/as_deg)*T-(0.0006/as_deg)*T**2+(0.00181/as_deg)*T**3 #solving for obliquity
    epsilon=radians(epsilon)
    beta=arcsin(sin(delta)*cos(epsilon)-cos(delta)*sin(epsilon)*sin(alpha))   #solving for ecliptic latituede (beta)
    lamb=arccos(cos(alpha)*cos(delta)/cos(beta)  ) ##solving for elcliptic longitude(lambda)    

    ##    finding change in ecliptic coordinates due to parallax
    do= 2447161.5  
    dd=t-do
    ee= 0.016714
    gg= -0.04144+0.0172019696*dd
    L= -1.38691+0.0172021240*dd
    lambda_sun=L +2.0*ee*sin(gg)+1.2*ee**2*sin(2.0*gg)
        

    deltalambda=parallax*(sin(lambda_sun-lamb))/cos(beta)
    deltabeta=parallax*cos(lambda_sun-lamb)*sin(beta)

        
    lambda_new=lamb+deltalambda
    beta_new=beta+deltabeta

    #transform from ecliptic to RA and DEC due to parallax
    delta_new=arcsin(sin(beta_new)*cos(epsilon)+cos(beta_new)*sin(epsilon)*sin(lambda_new))
    alpha_new=arccos(cos(lambda_new)*cos(beta_new)/cos(delta_new))



    ##convert RA and DEC from radians to mas
        
    RA_parallax=degrees(alpha_new)
    DEC_parallax=degrees(delta_new)
    return(RA_parallax, DEC_parallax)
    
pos_err= 6*10**(-3) #mas
pm_err= 5*10**(-3) #mas
parallax_err= 6*10**(-3) #mas

RA_parallax,DEC_parallax=parallax(RA_p,DEC_p,t)
RA_parallax_pm,DEC_parallax_pm=parallax(RA_pm,DEC_pm,t)

RA_parallax_rand=RA_parallax[index_rand]
DEC_parallax_rand=DEC_parallax[index_rand]
RA_parallax_pm_rand=RA_parallax_pm[index_rand]
DEC_parallax_pm_rand=DEC_parallax_pm[index_rand]

deg_as=3600.0

#print(RA_p*deg_mas)
#print(DEC_p*deg_mas)
#print((max(DEC_p)-min(DEC_p))*deg_mas)

figure(1)
plt.subplot(311)
title('DEC vs RA')
plot(RA_p*deg_mas,DEC_p*deg_mas)
plot(RA_p_rand*deg_mas,DEC_p_rand*deg_mas,'.')
plt.errorbar(RA_p_rand*deg_mas, DEC_p_rand*deg_mas,xerr=pos_err , yerr=pos_err,fmt='none')

xlabel('RA (mas)')
ylabel('DEC (mas) ')

plt.subplot(312)
title('RA vs time')
plot(t,RA_p*deg_mas)
plot(trand, RA_p_rand*deg_mas,'.')
plt.errorbar(trand,RA_p_rand*deg_mas,yerr=pos_err,fmt='none')
xlabel('time in JD days')
ylabel('RA (mas)')

plt.subplot(313)
title('DEC vs time')
plot(t,DEC_p*deg_mas)
plot(trand,DEC_p_rand*deg_mas,'.')
plt.errorbar(trand,DEC_p_rand*deg_mas,yerr=pos_err,fmt='none')
xlabel('time in JD days')
ylabel('DEC(mas)')


figure(2)
plt.subplot(311)
title('DEC vs RA')
plot(RA_parallax*deg_mas,DEC_parallax*deg_mas)
plot(RA_parallax_rand*deg_mas,DEC_parallax_rand*deg_mas,'.')
plt.errorbar(RA_parallax_rand*deg_mas, DEC_parallax_rand*deg_mas,xerr=parallax_err , yerr=parallax_err,fmt='none')
xlabel('RA (mas)')
ylabel('DEC (mas) ')

plt.subplot(312)
title('RA vs time')
plot(t,RA_parallax*deg_mas)
plot(trand, RA_parallax_rand*deg_mas, '.')
plt.errorbar(trand,RA_parallax_rand*deg_mas, yerr=parallax_err,fmt='none')
xlabel('time in JD days')
ylabel('RA (mas)')

plt.subplot(313)
title('DEC vs time')
plot(t,DEC_parallax*deg_mas)
plot(trand, DEC_parallax_rand*deg_mas,'.')
plt.errorbar(trand,DEC_parallax_rand*deg_mas, yerr=parallax_err,fmt='none')
xlabel('time in JD days')
ylabel('DEC(mas)')

figure(3)
plt.subplot(311)
title('DEC vs RA')
plot(RA_parallax_pm*deg_mas,DEC_parallax_pm*deg_mas)
plot(RA_parallax_pm_rand*deg_mas,DEC_parallax_pm_rand*deg_mas,'.')
plt.errorbar(RA_parallax_pm_rand*deg_mas, DEC_parallax_pm_rand*deg_mas,xerr=pm_err , yerr=pm_err,fmt='.')
xlabel('RA (mas)')
ylabel('DEC (mas) ')

plt.subplot(312)
title('RA vs time')
plot(t,RA_parallax_pm*deg_mas)
plot(trand,RA_parallax_pm_rand*deg_mas,'.')
plt.errorbar(trand,RA_parallax_pm_rand*deg_mas,yerr=pm_err,fmt='none')
xlabel('time in JD days')
ylabel('RA (mas)')

plt.subplot(313)
title('DEC vs time')
plot(t,DEC_parallax_pm*deg_mas)
plot(trand,DEC_parallax_pm_rand*deg_mas,'.')
plt.errorbar(trand,DEC_parallax_pm_rand*deg_mas,yerr=pm_err,fmt='none')
xlabel('time in JD days')
ylabel('DEC(mas)')
show()

