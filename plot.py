#!/usr/bin/env python3

# import csv
import matplotlib.pyplot as plt
import scipy.optimize as sopt
import numpy as np
import pandas as pd


"""`SodsExact(...)` (with `pressureEqn(...)`) copied from
https://github.com/fergu/EulerWeno5/blob/master/Python/eulerweno_LF.py
"""
#def pressureEqn(p2p1,p4p1,g,a1a4):
    #num = (g-1)*a1a4*(p2p1-1);
    #den = np.sqrt(2*g*(2*g+(g+1)*(p2p1-1)))
    #powterm = np.power(1-num/den,-2*g/(g-1))
    #return p2p1*powterm-p4p1;


#def SodsExact(rho4,u4,p4,rho1,u1,p1,tend):
    ## (4) (expansion) (3) (contact surface) (2) (shock) (1)
    #gamma=1.4;
    #x0 = 0.5;
    #gpm = (gamma+1)/(gamma-1);

    #p4p1 = p4/p1; # Left/right pressure ratio
    #a4 = np.sqrt(gamma*p4/rho4); # Left speed of sound
    #a1 = np.sqrt(gamma*p1/rho1); # Right speed of sound
    #a4a1 = a4/a1; # Ratio of the speed of sounds
    #p2p1 = sopt.fsolve(pressureEqn,1.0,args=(p4p1,gamma,1.0/a4a1)) # This gives the pressure ratio of the shock
    #ushock = a1*np.sqrt((gamma+1)/(2*gamma)*(p2p1-1)+1)
    ## Calculate region between shock and contact surface
    #p2 = p2p1*p1;
    #rho2rho1 = (1+gpm*p2p1)/(gpm+p2p1); rho2 = rho2rho1*rho1;
    #u2 = a1/gamma*(p2p1-1)*np.sqrt(((2*gamma)/(gamma+1))/(p2p1+1.0/gpm));
    ## Now calculate the expansion wave
    #p3p4 = p2p1/p4p1; p3 = p3p4*p4;
    #rho3rho4 = np.power(p3p4,1/gamma); rho3 = rho3rho4*rho4;
    #u3 = u2; # Velocity is unchanged across the contact surface
    ## Now need where each of these regions are based on time
    #x1 = x0 - a4*tend; # Location of the left part of the expansion fan
    #a3 = np.sqrt(gamma*p3/rho3)
    #x2 = x0 + (u3-a3)*tend; # Location of the right part of the expansion fan
    #x3 = x0 + u2*tend; # Location of the contact surface
    #x4 = x0 + ushock*tend; # Location of the shock wave
    ## # Now calculate behavior in expansion region
    #xpts = np.linspace(0,1,1000);
    #ypts = np.zeros((3,len(xpts)))
    #for i in np.arange(0,len(xpts)):
        #x = xpts[i]
        #if (x <= x1): # In the left region
            #ypts[0,i] = rho4;
            #ypts[1,i] = rho4*u4;
            #ypts[2,i] = p4/(gamma-1)+rho4*np.power(u4,2.0)/2;
        #elif (x > x1 and x <= x2): # In the expansion fan
            #u = 2/(gamma+1)*(a4+(x-x0)/tend)
            #p = p4*np.power(1-(gamma-1)/2*(u/a4),(2*gamma)/(gamma-1))
            #rho = rho4*np.power(1-(gamma-1)/2*(u/a4),2/(gamma-1))
            #ypts[0,i] = rho;
            #ypts[1,i] = rho*u;
            #ypts[2,i] = p/(gamma-1)+rho*np.power(u,2.0)/2;
        #elif (x > x2 and x <= x3): # Between the expansion and contact surface
            #ypts[0,i] = rho3;
            #ypts[1,i] = rho3*u3;
            #ypts[2,i] = p3/(gamma-1)+rho3*np.power(u3,2.0)/2;
        #elif (x > x3 and x <= x4): # Between the contact surface and the shock
            #ypts[0,i] = rho2;
            #ypts[1,i] = rho2*u2;
            #ypts[2,i] = p2/(gamma-1)+rho2*np.power(u2,2.0)/2;
        #elif (x > x4): # The right region
            #ypts[0,i] = rho1;
            #ypts[1,i] = rho1*u1;
            #ypts[2,i] = p1/(gamma-1)+rho1*np.power(u1,2.0)/2;
	#return [xpts,ypts];

#[xex,yex] = SodsExact(1.,0,1.,0.125,0.,.1, t)


#Rhoex = yex[0,:]
#Uex = yex[1,:]/yex[0,:]
#Pex = (yex[2,:]-yex[1,:]*Uex/2)*(1.4-1)


x = []
rho = []
u = []
p = []
s = []

#with open('res.dat', 'r') as datafile:
    #plotting = csv.reader(datafile, delimiter=' ')

    #for ROWS in plotting:
        #x.append(float(ROWS[0]))
        #rho.append(float(ROWS[1]))
        #u.append(float(ROWS[2]))
        #p.append(float(ROWS[3]))
        #s.append(float(ROWS[4]))
n = 100
filename = f'res_n_{n}.dat'
df = pd.read_csv(filename,
                 skiprows = 3, sep=' ',
				 names = ['x', 'rho', 'u', 'p', 'e'])

with open(filename, 'r') as datafile:
    title = datafile.readline()
    t = float(title[title.find('t=')+2:title.rfind('"')])

# L = 1250
L = 1

x = list(map(float, df['x'].tolist()))
rho = list(map(float, df['rho'].tolist()))
u = list(map(float, df['u'].tolist()))
p = list(map(float, df['p'].tolist()))
e = list(map(float, df['e'].tolist()))



plt.figure(1)
plt.subplot(4, 1, 1)
# plt.subplot(3, 1, 1)
plt.title("t={:0.4f}".format(t))
# plt.plot(x, rho, '.b-')
plt.plot(x, rho, 'b-')
plt.grid(True, ls=':')
# plt.plot(x500, rho500, 'g-')
# plt.plot(x1000, rho1000, 'b-')
# plt.plot(xex,Rhoex,'r--', linewidth=0.5)
plt.ylabel('Density ρ')
# plt.xlim([0, L])
plt.subplot(4, 1, 2)
plt.grid(True, ls=':')
#plt.subplot(3, 1, 2)
# plt.plot(x, u, '.b-')
plt.plot(x, u, 'b-')
# plt.plot(x500, u500, 'g-')
# plt.plot(x1000, u1000, 'b-')
# plt.plot(xex,Uex,'r--', linewidth=0.5)
plt.ylabel('Velocity v')
# plt.xlim([0, L])
plt.subplot(4, 1, 3)
plt.grid(True, ls=':')
# plt.subplot(3, 1, 3)
# plt.plot(x, p, '.b-')
plt.plot(x, p, 'b-')
# plt.plot(x500, p500, 'g-')
# plt.plot(x1000, p1000, 'b-')
# plt.xlim([0, L])
# plt.plot(xex,Pex,'r--', linewidth=0.5)
plt.ylabel('Pressure p')
plt.subplot(4, 1, 4)
plt.grid(True, ls=':')
# plt.plot(x, e, '.b-');
plt.plot(x, e, 'b-');
#plt.plot(x, [1 - dp for dp in s], 'r-');
#plt.xlim([0,L])
#plt.ylabel('Species Fraction')
plt.ylabel('Energy e')


plt.tight_layout()

plt.savefig(filename[:-4] + '.svg', format='svg', transparent=True)
plt.show();
