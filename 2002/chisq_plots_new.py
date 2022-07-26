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

output_file = open("PJDFC.out")
output_file.readline()

diameters = []
chis = []
gammas = []

for line in output_file.readlines():
    datum = line.strip().split()
    if len(datum) < 10:
        #capture odd cases where col 6 was runing into col 5.  Only a handful, so just skip 
        continue
    
    diameters.append(numpy.e ** float(datum[5]))
    chis.append(float(datum[8]))
    gammas.append(numpy.e ** float(datum[4]))


plt.figure(figsize=[6, 6])
plt.subplots_adjust(left=0.15, right=0.95, top=0.95)
plt.loglog(diameters, chis, color='k', marker='o', ms=0.5, ls="None")
plt.xlabel("Diameter (km)", fontsize=13)
plt.ylabel(r"fit $\chi^2$", fontsize=13)
plt.savefig("diameter_vs_chi.png")
plt.savefig("diameter_vs_chi.pdf")
plt.close()

plt.figure(figsize=[6, 6])
plt.subplots_adjust(left=0.15, right=0.95, top=0.95)
plt.plot(diameters, chis, color='k', marker='o', ms=0.5, ls='None')
plt.xlabel("Diameter (km)", fontsize=13)
plt.ylabel(r"fit $\chi^2$", fontsize=13)
if "Phaethon" in os.getcwd():
    plt.axvline(x=4.6, ls="dashed", color='k')
    plt.axvline(x=4.8, ls="dotted", color='k')
    plt.axvline(x=4.3, ls="dotted", color='k')
    plt.xlim(3.5, 5.5)
    plt.ylim(55, 80)
elif "2005UD" in os.getcwd():
    plt.axvline(x=1.2, ls="dashed", color='k')
    plt.axvline(x=0.8, ls="dotted", color='k')
    plt.axvline(x=1.6, ls="dotted", color='k')
    plt.xlim(0.5, 2.5)
    plt.ylim(7, 25)
else:    
    plt.ylim(0.9 * min(chis), 2 * min(chis))
    
    
plt.savefig("diameter_vs_chi_zoom.png")
plt.savefig("diameter_vs_chi_zoom.pdf")
plt.close()


plt.figure(figsize=[6, 6])
plt.subplots_adjust(left=0.15, right=0.95, top=0.95)
plt.loglog(gammas, chis, color='k', marker='o', ms=0.5, ls='None')
plt.xlabel(r"Thermal inertia ($J~m^{-2}~s^{-0.5}~K^{-1})$)", fontsize=13)
plt.ylabel(r"fit $\chi^2$", fontsize=13)

if "Apophis" in os.getcwd():
    plt.xlim(100, 30000)
    plt.ylim(5, 30)

plt.savefig("gamma_vs_chi.png")
plt.savefig("gamma_vs_chi.pdf")
plt.close()




exit()

logD=[numpy.log10(d) for d in D]
foo=logD[:]
nD=len(logD)
foo.sort()
sig1low=foo[int(nD*0.16)]
sig1high=foo[int(nD*0.84)]
sig2low=foo[int(nD*0.025)]
sig2high=foo[int(nD*0.975)]

if sig2high-sig2low > 0.2:
    histstep=0.01
elif sig2high-sig2low > 0.02:
    histstep=0.001
else:
    histstep=0.0001

plt.figure(figsize=[6,6])
histlimlow=int(min(logD)*1000)/1000.
histlimhigh=(int(max(logD)*1000)+1)/1000.
plt.hist(logD,bins=numpy.arange(histlimlow,histlimhigh,histstep),histtype='step',color='black')


plt.axvline(sig1low,color='#cc0000',ls='dashed')
plt.axvline(sig1high,color='#cc0000',ls='dashed')
plt.axvline(sig2low,color='#cc0000',ls='dotted')
plt.axvline(sig2high,color='#cc0000',ls='dotted')

plt.xlabel("log diameter (km)")
plt.ylabel("number of Monte Carlo results")
plt.savefig("diam_hist.png")


Dper=[]
period=[]
inf=open("DvsPeriod.dat")
for line in inf.readlines():
    dat=line.rstrip().split()
    Dper.append(float(dat[0]))
    period.append(float(dat[1]))

if max(period)-min(period)==0:
    print("Fixed period used, skipping D vs Period plot")
else:

    plt.figure(figsize=[6,6])
    
    plt.scatter(Dper,period,s=1.5,c='black',marker='o',edgecolors='')
    
    (Hbin,xbin,ybin)=numpy.histogram2d(Dper,period,bins=20)
    xbino=[(xbin[i+1] - xbin[i])/2. + xbin[i] for i in range(len(xbin)-1)]
    ybino=[(ybin[i+1] - ybin[i])/2. + ybin[i] for i in range(len(ybin)-1)]
    Hout=numpy.transpose(Hbin)
    Hmax=numpy.amax(Hout)
    
    foo=list(Hout.flatten())
    foo.sort()
    totH=sum(foo)
    
    sumH=0
    sig1=-1
    sig2=-1
    for i in range(len(foo)):
        sumH+=foo[i]
        if sumH>0.05*totH and sig2<0:
            sig2=i
        if sumH>0.32*totH and sig1<0:
            sig1=i
    
    plt.contour(xbino,ybino,Hout,colors='#cc0000',levels=[foo[sig2],foo[sig1]])
    plt.text(numpy.mean(Dper),0.95*max(period),"Contours contain 68%\nand 95.5% of all points")
    plt.ylim(-2,0)
    plt.xlabel("diameter (km)")
    plt.ylabel("rotation period")
    plt.savefig("DvsPer.png")
    plt.close()
