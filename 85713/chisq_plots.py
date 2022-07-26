#python3 parser of the output from WISE-steep-MCMC-PJDFC
#run this in the directory made for the object

import sys
import os
import matplotlib.pyplot as plt
import numpy
from math import modf,floor
if sys.version[0] != '3':
    print("Please run with python3")
    exit()

inf=open("PJDFC.out")
inf.readline()

D=[]
chi=[]
gamma=[]

for line in inf.readlines():
    dat=line.strip().split()
    if len(dat)<10:
        #capture odd cases where col 6 was runing into col 5.  Only a handful, so just skip 
        continue
    
    D.append(numpy.e**float(dat[5]))
    chi.append(float(dat[8]))
    gamma.append(numpy.e**float(dat[4]))


plt.figure(figsize=[6,6])
plt.subplots_adjust(left=0.15,right=0.95,top=0.95)
plt.loglog(D,chi,color='k',marker='o',ms=0.5,ls='None')
plt.xlabel("diameter (km)",fontsize=13)
plt.ylabel(r"fit $\chi^2$",fontsize=13)
plt.savefig("DvChi.png")
plt.savefig("DvChi.pdf")
plt.close()

plt.figure(figsize=[6,6])
plt.subplots_adjust(left=0.15,right=0.95,top=0.95)
plt.plot(D,chi,color='k',marker='o',ms=0.5,ls='None')
plt.xlabel("diameter (km)",fontsize=13)
plt.ylabel(r"fit $\chi^2$",fontsize=13)
if 'Phaethon' in os.getcwd():
    plt.axvline(x=4.6,ls='dashed',color='k')
    plt.axvline(x=4.8,ls='dotted',color='k')
    plt.axvline(x=4.3,ls='dotted',color='k')
    

    plt.xlim(3.5,5.5)
    plt.ylim(55,80)
elif '2005UD' in os.getcwd():
    plt.axvline(x=1.2,ls='dashed',color='k')
    plt.axvline(x=0.8,ls='dotted',color='k')
    plt.axvline(x=1.6,ls='dotted',color='k')
    plt.xlim(0.5,2.5)
    plt.ylim(7,25)
else:    
    plt.ylim(0.9*min(chi),2*min(chi))
    
    
plt.savefig("DvChi_zoom.png")
plt.savefig("DvChi_zoom.pdf")
plt.close()


plt.figure(figsize=[6,6])
plt.subplots_adjust(left=0.15,right=0.95,top=0.95)
plt.loglog(gamma,chi,color='k',marker='o',ms=0.5,ls='None')
plt.xlabel(r"Thermal inertia ($J~m^{-2}~s^{-0.5}~K^{-1})$)",fontsize=13)
plt.ylabel(r"fit $\chi^2$",fontsize=13)

if 'Apophis' in os.getcwd():
#    plt.axhline(y=6.5,ls='dotted',color='#444444')
    plt.axhline(y=7.5,ls='dotted',color='#444444')
#    plt.axvline(x=600,ls='dotted',color='#444444')
#    plt.axvline(x=5000,ls='dotted',color='#444444')
    plt.xlim(100,30000)
    plt.ylim(5,30)

plt.savefig("GammavChi.png")
plt.savefig("GammavChi.pdf")
plt.close()

###Working on it

plt.scatter(D, gamma, marker='o',color='k')
plt.title('Diameters vs Gamma Plot')
plt.xlabel('Diameter (km)')
plt.ylabel('Gamma')
plt.savefig("diam_gamma_plot.png")