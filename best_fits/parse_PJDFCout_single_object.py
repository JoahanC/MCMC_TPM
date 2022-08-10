#python3 parser of the output from WISE-steep-MCMC-PJDFC
#run this in the directory made for the object

import os, sys, numpy, argparse, matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from astropy.table import Table
from pathlib import Path

plt.rcParams['axes.facecolor'] = 'whitesmoke'
COLOR = 'gray'
mpl.rcParams['text.color'] = COLOR
mpl.rcParams['axes.labelcolor'] = 'dimgray'
mpl.rcParams['xtick.color'] = COLOR
mpl.rcParams['ytick.color'] = COLOR
plt.rcParams["font.family"] = "monospace"

class Asteroid(object):
    def __init__(self, fname=None):
        #print(fname)
        result = open(fname)
        for line in result:
            x = line.split()
            #print(x)
            if "dia=" in line:
                # print(x)
                diameter_median = float(x[4].split("+")[0])
                diameter_84th_percentile = float(x[5].split("-")[0])
                diameter_16th_percentile = float(x[6].split("%")[0])
                self.diameter_median = diameter_median
                self.diameter_16th_percentile = diameter_16th_percentile
                self.diameter_84th_percentile = diameter_84th_percentile

            if "p_V = " in line:
                # print(x)
                albedo_median = float(x[6])
                albedo_84th_percentile = float(x[8].split("%")[0])/100.
                albedo_16th_percentile = float(x[8].split("%")[0])/100.
                # print(albedo_median, albedo_84th, albedo_16th)
                self.albedo_median = albedo_median
                self.albedo_16th_percentile = albedo_16th_percentile
                self.albedo_84th_percentile = albedo_84th_percentile

            if "theta1=" in line:
                # print(x)
                theta_median = float(x[4].split("+")[0])
                theta_84th_percentile = float(x[5].split("-")[0])
                theta_16th_percentile = float(x[6].split("%")[0])
                self.theta_median = theta_median
                self.theta_16th_percentile = theta_16th_percentile
                self.theta_84th_percentile = theta_84th_percentile

            if "Period [h]" in line:
                # print(x)
                period_median = float(x[3].split("+")[0])
                period_84th_percentile = float(x[4])
                period_16th_percentile = float(x[5].split("%")[0])
                # print(period_median, period_84th, period_16th)
                self.period_median = period_median
                self.period_16th_percentile = period_16th_percentile
                self.period_84th_percentile = period_84th_percentile

            if "sqrt(kappa*rho*C)" in line:
                #print(line)
                #print(float(line.split("+")[1].split("-")[0]))
                gamma_median = float(line.split("=")[1].split("+")[0])
                gamma_84th_percentile = float(line.split("+")[1].split("-")[1].split("%")[0])/100.
                gamma_16th_percentile = float(line.split("+")[1].split("-")[0])/100.
                # print(gamma_median, gamma_84th, gamma_16th)
                self.gamma_median = gamma_median
                self.gamma_16th_percentile = gamma_16th_percentile
                self.gamma_84th_percentile = gamma_84th_percentile

            if "crater fraction" in line:
                # print(x)
                crater_fraction_median = float(x[2].split("+")[0])
                crater_fraction_84th_percentile = float(x[2].split("+")[1].split("-")[0])
                crater_fraction_16th_percentile = float(x[2].split("-")[1])
                # print(crater_fraction_median, crater_fraction_84th, crater_fraction_16th)
                self.crater_fraction_median = crater_fraction_median
                self.crater_fraction_16th_percentile = crater_fraction_16th_percentile
                self.crater_fraction_84th_percentile = crater_fraction_84th_percentile

            if "p_IR/p_V=" in line:
                # print(x)
                ir_fraction_median = float(x[1].split("+")[0])
                ir_fraction_84th_percentile = float(x[1].split("+")[1].split("-")[0])
                ir_fraction_16th_percentile = float(x[1].split("-")[1].split("%")[0])
                # print(ir_fraction_median, ir_fraction_84th, ir_fraction_16th)
                self.ir_fraction_median = ir_fraction_median
                self.ir_fraction_16th_percentile = ir_fraction_16th_percentile
                self.ir_fraction_84th_percentile = ir_fraction_84th_percentile

            if "pole peak at = " in line:
                # print(x)
                pole_peak_ra = float(x[4])
                pole_peak_dec = float(x[5])
                # print(pole_peak_ra, pole_peak_dec)
                self.pole_peak_ra = pole_peak_ra
                self.pole_peak_dec = pole_peak_dec

            if "mean pole at = " in line:
                pole_mean_ra = float(x[4])
                pole_mean_dec = float(x[5])
                pole_bar = float(x[7])
                # print(pole_mean_ra, pole_mean_dec, pole_bar)
                self.pole_mean_ra = pole_mean_ra
                self.pole_mean_dec = pole_mean_dec
                self.pole_bar = pole_bar


def merger(A,B):
    '''
    Function that couples values from two lists together
    :param A: list A
    :param B: list B
    :return: [(Ai, Bi): where i = len(list)]
    '''
    answer = []
    for l in range(len(A)):
        a=A[l]
        b=B[l]
        array = [a,b]
        answer.append(array)
    return answer

def med1sigma(parameter,weights):
    '''
    Function that computes the weighted values, weighted median, & 16th & 84th percentiles
    for a parameter from the TPM.
    :param parameter:
    :param weights:
    :return:
    '''

    arr = merger(parameter, weights)
    arr = sorted(arr, key=lambda x: x[0], reverse=False)

    # appending sorted parameter and weights into separate lists for ease of manipulation
    value = [arr[i][0] for i in range(len(parameter))]
    weights2 = [arr[i][1] for i in range(len(parameter))]

    #print("weights2",numpy.shape(weights2),weights2[0:10])
    test_weighted_points = []
    #test_weighted_points2 += [value[i]]*int(weights2[i])
    for i in range(0, len(weights2)):
        test_weighted_points.append([value[i]]*int(weights2[i])) # the number of points will increase
        #print(value[i], (str(weights2[i])))

    sorted_V1 = sorted([item for sublist in test_weighted_points for item in sublist])

    nJ = len(sorted_V1)
    sig1low=sorted_V1[int(nJ*0.16)]
    sig1high=sorted_V1[int(nJ*0.84)]
    med=sorted_V1[int(nJ*0.5)]
    sig2low=sorted_V1[int(nJ*0.025)]
    sig2high=sorted_V1[int(nJ*0.975)]

    return sorted_V1, med, sig1low, sig1high, sig2low, sig2high

if sys.version[0] != '3':
    print("Please run with python3")
    exit()

# plt.rcParams['axes.facecolor'] = 'whitesmoke'
# COLOR = 'gray'
# matplotlib.rcParams['text.color'] = COLOR
# matplotlib.rcParams['axes.labelcolor'] = 'dimgray'
# matplotlib.rcParams['xtick.color'] = COLOR
# matplotlib.rcParams['ytick.color'] = COLOR
# plt.rcParams["font.family"] = "monospace"
#sns.set_theme(style="ticks")
#sns.set_context("talk")

parser = argparse.ArgumentParser(description="Run the TPM parser on an individual object")
parser.add_argument("--path", dest="path", action="store", help="the directory where the object is located", default=None)
parser.add_argument("--name", dest="name", action="store", help="the name of the object", default=False)
args = parser.parse_args()

path = args.path
name = args.name


possible_symbols = ["o", "s", "v", "*", "^", "D", "h", "X", ">", "<", "P", "8", "x", 4]
possible_colors = ["blue", "deepskyblue", "cyan", "lime", "green", "gold", "orange", "tomato",
                   "red", "pink", "magenta", "purple", "black", "silver"]
objects = {}


# for item in dirs:
#i = int(args.dir)-1
# for item in [dirs[i]]:
dirs = [path]
name = args.name

#os.system("cd "+ path + item)
os.chdir(path)
print("name", name)

#cp, not mv, so that there is a record of the actual outputs
os.system("/bin/cp "+path+"/fort.8 "+path+"/PJDFC_patch.out")
#fort.22 is the postscript format of the SED
os.system("/bin/cp "+path+"/fort.22 "+path+"/SED_data.out")


#read the SED file here and load it into a list of lines
listOfLines = list()
# with open (path+item+"/SED_data.out", "r") as inFile:
with open (path+"/SED_data.out", "r") as inFile:
    for line in inFile:
        listOfLines.append(line.strip())

#get the wavelengths for the SED plots
start = listOfLines.index("/wvl [")
# print(listOfLines[6])

end = listOfLines.index("] def")
# print(listOfLines[end])

wave=[]
for l in listOfLines[start+1:end]:
    wdat=l.strip().split()
    for w in wdat:
        wave.append(float(w))
# print("wave", wave)

#now get the predicted SED fluxes for each epoch
starts = [i for i, e in enumerate(listOfLines) if "/ft" in e]
ends = [i for i, e in enumerate(listOfLines) if "] def doit" in e]
esed = []
for j in range(0, len(starts)):
    # print("epoch %d" % j)
    # print([float(x) for x in listOfLines[starts[j]+1:ends[j]]])
    esed.append([float(x) for x in listOfLines[starts[j]+1:ends[j]]])

# print("END********")
# print(esed)

#now get the **measured** fluxes for each epoch
start = listOfLines.index("] def")
# print("start", start)
starts = [i for i, e in enumerate(listOfLines) if "def" in e]
# print(starts)
ends = [i for i, e in enumerate(listOfLines) if "/ft" in e]
# print(ends)

nepoch = len(starts)-1
obsdatax = []
obsdatayerr = []
obsdatay = []
colors = []
symbols = []

for j in range(0,len(starts)-1):
    # print("epoch", j)
    # print(listOfLines[s])
    someshit = [e for i, e in enumerate(listOfLines[starts[j]:ends[j]]) if "QQ" in e]
    # print(someshit)
    print("number of bands with measurements", len(someshit))
    bands_available_for_this_epoch = []
    errors_available_for_this_epoch = []
    fluxes_available_for_this_epoch = []
    # colors_available_for_this_epoch = []
    # symbols_available_for_this_epoch = []

    for s in someshit:
        (xx,err,yy,foo) = s.split()
        # print(xx, err, yy)
        bands_available_for_this_epoch.append(float(xx))
        errors_available_for_this_epoch.append(float(err))
        fluxes_available_for_this_epoch.append(float(yy))
        # symbols_available_for_this_epoch.append(possible_symbols[j])
        # colors_available_for_this_epoch.append(possible_colors[j])

        # if 3. < float(xx) < 4.:
        #     colors_available_for_this_epoch.append("b")
        # elif 4. < float(xx) < 5.:
        #     colors_available_for_this_epoch.append("c")
        # elif 11. < float(xx) < 12.:
        #     colors_available_for_this_epoch.append("r")
        # else:
        #     colors_available_for_this_epoch.append("m")
        #

    obsdatax.append(bands_available_for_this_epoch)
    obsdatayerr.append(errors_available_for_this_epoch)
    obsdatay.append(fluxes_available_for_this_epoch)
    # colors.append(colors_available_for_this_epoch)
    # symbols.append(symbols_available_for_this_epoch)
    colors.append(possible_colors[j])
    symbols.append(possible_symbols[j])


# print(obsdatax)
# print(obsdatayerr)
# print(obsdatay)
# print(colors)
# print(symbols)
print("there are %d epochs" % nepoch)

plt.figure(figsize=[9, 9])
ax=plt.gca()
for epoch in range(nepoch):

    # print("epoch ", epoch, obsdatax[epoch], colors[epoch], symbols[epoch])
    plt.plot(wave,esed[epoch],color=colors[epoch],ls='solid')

    ax.set_xscale('log')
    ax.set_yscale('log')

    erroroutl = [obsdatay[epoch][n] * (1 - 1 / (10 ** (obsdatayerr[epoch][n] / 2.5))) for
                 n in range(len(obsdatay[epoch]))]
    errorouth = [obsdatay[epoch][n] * (10 ** (obsdatayerr[epoch][n] / 2.5) - 1) for
                 n in range(len(obsdatay[epoch]))]

    the_label = "Epoch " + str(int(epoch) + 1)
    plt.errorbar(obsdatax[epoch], obsdatay[epoch], yerr=[erroroutl, errorouth],
                 color=colors[epoch], marker=symbols[epoch], linestyle="none",
                 label=the_label)


#plt.title(name)
plt.axvspan(3.14, 3.78, facecolor='c', alpha=0.2)
plt.axvspan(4.09, 5.19, facecolor='b', alpha=0.2)
plt.xlabel("log Wavelength (microns)",fontsize=15)
plt.ylabel(r"log $\nu$F$_\nu$ (erg/cm$^{2}$/s)",fontsize=15)
plt.legend(loc="lower right", prop={'size': 10})
plt.tight_layout()
plt.savefig(path+"/plots/bestfit_SED.pdf")
plt.close()


#always rerun, to get analysis parameters
if os.path.isfile("fort.4"):
    print("PJDFC parser output files already exist; using those")
else:

# ostring = "echo "+item+"/PJDFC.out | /Users/amainzer/PycharmProjects/pyamysandbox/"
# ostring += "TPM/mcmc_tpm-master/Code/read-WISE-rc-MCMC-PJDFC > "+name+"_out.txt"
    ostring = "echo "+path+"/PJDFC_patch.out | /home/u21/satpathyakash/MCMC/" #removed "+item" after path, unsure
    ostring += "TPM/read-WISE-rc-MCMC-PJDFC > "+name+"_out.txt"
    os.system(ostring)

os.system("/bin/cp "+path+"/fort.2 "+path+"/Dhist.dat")
os.system("/bin/cp "+path+"/fort.32 "+path+"/Dhist_fine.dat")
os.system("/bin/cp "+path+"/fort.3 "+path+"/DvsPeriod.dat")
os.system("/bin/cp "+path+"/fort.4 "+path+"/DvsAlb.dat")

#load all final results from the output file into a dictionary
#fname = path + "/" + name + "_out.txt"
fname = path + "/PJDFC_patch.out"
#print("fname", fname)
objects[name] = Asteroid(fname)
#print(name, objects[name].diameter_median, objects[name].diameter_16th_percentile)

#Dhist:
#Ned's code spits out a list of numbers (index 1 below)
#entry 1 = 1m
#entry 24 = 10 m
#entry 47 = 100 m
#entry 70 = 1 km
#final entry 123 = 200 km
#Dhist_fine has a factor of 10 more bins, each 1/10th the width
#But, I'm just going to use the diameters from DvsAlb to build my own histograms


# D=[]
# pV=[]
# inf=open(path+"/DvsAlb.dat")
# for line in inf.readlines():
#     dat=line.rstrip().split()
#     D.append(float(dat[0]))
#     pV.append(float(dat[1]))

inf=open(path+"/PJDFC_patch.out")
inf.readline()

# opening the output file with all MC results
inf2=open(path + "../PJDFC.out")
inf2.readline()

D=[]
D_full = []
pV=[]
pV_full = []
pIRpV=[]
gamma=[]
gamma_full = []
fc = []
fc_full = []
weights=[]
weights_full=[]

for line in inf.readlines():
    dat=line.strip().split()
    if len(dat)<10:
        #capture odd cases where col 6 was runing into col 5.  Only a handful, so just skip
        continue
    D.append(numpy.e**(float(dat[5])) * 1000)
    pV.append(numpy.e**(float(dat[2])))
    gamma.append(numpy.e**(float(dat[4])))
    pIRpV.append(numpy.e**(float(dat[7])))
    weights.append(int(dat[9]))
    fc.append(int(dat[6]))

for line2 in inf2.readlines():
    dat2=line2.strip().split()
    if len(dat2)<10:
        #capture odd cases where col 6 was runing into col 5.  Only a handful, so just skip
        continue
    D_full.append(numpy.e**(float(dat2[5])) * 1000)
    pV_full.append(numpy.e**(float(dat2[2])))
    gamma_full.append(numpy.e**(float(dat2[4])))
    #pIRpV.append(numpy.e**(float(dat2[7])))
    weights_full.append(int(dat2[9]))
    fc_full.append(dat2[6])


plt.figure(figsize=[9,9])
ax=plt.gca()
plt.scatter(D,pV,s=1.5,c='black',marker='o',edgecolor='none')
# ax.set_xscale('log')
ax.set_yscale('log')

(Hbin,xbin,ybin)=numpy.histogram2d(D,pV,bins=30)
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

#print(sig2,foo[sig2],sig1,foo[sig1])

#plt.contour(xbino,ybino,Hout,colors='#777777',levels=[0.25*Hmax,0.5*Hmax,0.75*Hmax])
#plt.text(0.2,-0.2,"Contours at 25%,\n50%, 75% of peak")
plt.contour(xbino,ybino,Hout,colors='#cc0000',levels=[foo[sig2],foo[sig1]])
plt.text(max(D)*.7,0.6,"Contours contain 68%\nand 95.5% of all points")
#plt.text(max(D)*0.7,0.5,"Best-fit $p_{v}$ =%5.2f" % objects[name].albedo_median)
#plt.text(max(D)*0.7,0.4,"Best-fit $D$ =%5.2f km" % objects[name].diameter_median)
#plt.ylim(-2,0)
plt.ylim(0.12,1)
plt.title(name)
plt.xlabel("diameter (m)", fontsize=15)
plt.ylabel("albedo", fontsize=15)
plt.tight_layout()
#plt.savefig(path+"/plots/DvsAlb.pdf")
plt.close()

Dper=[]
period=[]
inf=open(path+"/DvsPeriod.dat")
for line in inf.readlines():
    dat=line.rstrip().split()
    Dper.append(float(dat[0]))
    period.append(float(dat[1]))

# print("period", len(period))
# print("Dper", len(Dper))

# print(min(period))

if max(period)-min(period)==0:
    print("Fixed period used, skipping D vs Period plot")
else:

    plt.figure(figsize=[6,6])

    plt.scatter(Dper,period,s=1.5,c='black',marker='o',edgecolor='none')

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
    plt.yscale("log")
    plt.ylim(.01,1e4)
    plt.title(name)
    plt.xlabel("Diameter (m)", fontsize=15)
    plt.ylabel("Rotation Period (h)", fontsize=15)
    plt.tight_layout()
    #plt.savefig(path+"/plots/DvsPer.pdf")
    plt.close()


#### Make diameter histogram ####
weighted_diameter,sigmed,sig1low,sig1high,sig2low,sig2high = med1sigma(D, weights)
logD = [numpy.log10(d) for d in weighted_diameter]
print("len(D), len(weighted_diameter)", len(D), len(weighted_diameter))
print(sigmed,sig1low,sig1high)
sigmed = numpy.log10(sigmed)
sig1low = numpy.log10(sig1low)
sig1high = numpy.log10(sig1high)
sig2low = numpy.log10(sig2low)
sig2high = numpy.log10(sig2high)

weighted_diameter_full,sigmed_full,sig1low_full,sig1high_full,sig2low_full,sig2high_full = med1sigma(D_full,weights_full)
logD_full = [numpy.log10(d_full) for d_full in weighted_diameter_full]
print("len(D_full), len(weighted_diameter_full)", len(D_full), len(weighted_diameter_full))
print(sigmed_full,sig1low_full,sig1high_full)
sigmed_full = numpy.log10(sigmed_full)
sig1low_full = numpy.log10(sig1low_full)
sig1high_full = numpy.log10(sig1high_full)
sig2low_full = numpy.log10(sig2low_full)
sig2high_full = numpy.log10(sig2high_full)

if sig2high-sig2low > 0.2:
    histstep=0.01
elif sig2high-sig2low > 0.02:
    histstep=0.001
else:
    histstep=0.0001

if sig2high_full-sig2low_full > 0.2:
    histstep_full=0.01
elif sig2high_full-sig2low_full > 0.02:
    histstep_full=0.001
else:
    histstep_full=0.0001

plt.figure(figsize=[9,9])
histlimlow=int(min(logD)*1000)/1000.
histlimhigh=(int(max(logD)*1000)+1)/1000.

histlimlow_full=int(min(logD_full)*1000)/1000.
histlimhigh_full=(int(max(logD_full)*1000)+1)/1000.

#plt.hist(logD,bins=numpy.arange(histlimlow,histlimhigh,histstep),
#         histtype='step',color='black',label="Best-fit D = %5.3f km"
#            % 10**sigmed)
plt.hist(logD,bins=numpy.arange(histlimlow,histlimhigh,histstep),
         histtype='step',color='black', edgecolor='#cc0000')
plt.hist(logD_full,bins=numpy.arange(histlimlow_full,histlimhigh_full,histstep_full),
         histtype='step',color='black', lw = 0.25)
print("length of logD:",len(logD), "& length of logD_full:",len(logD_full))

plt.axvline(sig1low,color='#cc0000',ls='dashed')
plt.axvline(sig1high,color='#cc0000',ls='dashed')
plt.axvline(sigmed, color="#cc0000",ls="solid")
plt.axvline(sig2low,color='#cc0000',ls='dotted')
plt.axvline(sig2high,color='#cc0000',ls='dotted')
#plt.title(name)
plt.xlabel("log diameter (m)", fontsize=15)
plt.ylabel("Number of Monte Carlo results", fontsize=15)
#plt.legend(loc="upper right")
plt.tight_layout()
plt.savefig(path+"/plots/diam_hist.pdf")
plt.close()

#### Make albedo histogram ####

weighted_alb,sigmed_alb,sig1low_alb,sig1high_alb,sig2low_alb,sig2high_alb = med1sigma(pV, weights)
logAlb = [numpy.log10(alb) for alb in weighted_alb]
#print("len(D), len(weighted_diameter)", len(D), len(weighted_diameter))
#print(sigmed,sig1low,sig1high)
sigmed = numpy.log10(sigmed_alb)
sig1low = numpy.log10(sig1low_alb)
sig1high = numpy.log10(sig1high_alb)
sig2low = numpy.log10(sig2low_alb)
sig2high = numpy.log10(sig2high_alb)

weighted_alb_full,sigmed_alb_full,sig1low_alb_full,sig1high_alb_full,sig2low_alb_full,sig2high_alb_full = med1sigma(pV_full, weights_full)
logAlb_full = [numpy.log10(alb_full) for alb_full in weighted_alb_full]
#print("len(D_full), len(weighted_diameter_full)", len(D_full), len(weighted_diameter_full))
#print(sigmed_full,sig1low_full,sig1high_full)
sigmed_full = numpy.log10(sigmed_alb_full)
sig1low_full = numpy.log10(sig1low_alb_full)
sig1high_full = numpy.log10(sig1high_alb_full)
sig2low_full = numpy.log10(sig2low_alb_full)
sig2high_full = numpy.log10(sig2high_alb_full)

if sig2high-sig2low > 0.2:
    histstep=0.01
elif sig2high-sig2low > 0.02:
    histstep=0.001
else:
    histstep=0.0001

if sig2high_full-sig2low_full > 0.2:
    histstep_full=0.01
elif sig2high_full-sig2low_full > 0.02:
    histstep_full=0.001
else:
    histstep_full=0.0001

plt.figure(figsize=[9,9])
histlimlow=int(min(logAlb)*1000)/1000.
histlimhigh=(int(max(logAlb)*1000)+1)/1000.

histlimlow_full=int(min(logAlb_full)*1000)/1000.
histlimhigh_full=(int(max(logAlb_full)*1000)+1)/1000.

#plt.hist(logD,bins=numpy.arange(histlimlow,histlimhigh,histstep),
#         histtype='step',color='black',label="Best-fit D = %5.3f km"
#            % 10**sigmed)
plt.hist(logAlb,bins=numpy.arange(histlimlow,histlimhigh,histstep),
         histtype='step',color='black', edgecolor='darkblue')
plt.hist(logAlb_full,bins=numpy.arange(histlimlow_full,histlimhigh_full,histstep_full),
         histtype='step',color='black', lw = 0.25)
print("length of logD:",len(logD), "& length of logD_full:",len(logD_full))

plt.axvline(sig1low,color='darkblue',ls='dashed')
plt.axvline(sig1high,color='darkblue',ls='dashed')
plt.axvline(sigmed, color="darkblue",ls="solid")
plt.axvline(sig2low,color='darkblue',ls='dotted')
plt.axvline(sig2high,color='darkblue',ls='dotted')
#plt.title(name)
plt.xlabel("log Visual Albedo", fontsize=15)
plt.ylabel("Number of Monte Carlo results", fontsize=15)
#plt.legend(loc="upper right")
plt.tight_layout()
plt.xlim(-0.85, -0.0001)
plt.savefig(path+"/plots/albedo_hist.pdf")
plt.close()


#### Make the gamma histogram plot ####
weighted_gamma,sigmed,sig1low,sig1high,sig2low,sig2high = med1sigma(gamma,weights)

weighted_gamma_full,sigmed_full,sig1low_full,sig1high_full,sig2low_full,sig2high_full = med1sigma(gamma_full,weights_full)

logJ=[numpy.log10(g) for g in weighted_gamma]
sig1low = numpy.log10(sig1low)
sig1high = numpy.log10(sig1high)
sigmed = numpy.log10(sigmed)
print("gamma:", len(gamma), sigmed, sig1low, sig1high)
print("len(gamma), len(weighted_gamma)", len(gamma), len(weighted_gamma))

logJ_full=[numpy.log10(g_full) for g_full in weighted_gamma_full]
sig1low_full = numpy.log10(sig1low_full)
sig1high_full = numpy.log10(sig1high_full)
sigmed_full = numpy.log10(sigmed_full)
print("gamma:", len(gamma_full), sigmed_full, sig1low_full, sig1high_full)
print("len(gamma), len(weighted_gamma)", len(gamma_full), len(weighted_gamma_full))

if sig2high-sig2low > 0.2:
    histstep=0.01
elif sig2high-sig2low > 0.02:
    histstep=0.001
else:
    histstep=0.0001

if sig2high_full-sig2low_full > 0.2:
    histstep_full=0.01
elif sig2high_full-sig2low_full > 0.02:
    histstep_full=0.001
else:
    histstep_full=0.0001

plt.figure(figsize=[9,9])
histlimlow=int(min(logJ)*1000)/1000.
histlimhigh=(int(max(logJ)*1000)+1)/1000.

histlimlow_full=int(min(logJ_full)*1000)/1000.
histlimhigh_full=(int(max(logJ_full)*1000)+1)/1000.

#plt.hist(logJ,bins=numpy.arange(histlimlow,histlimhigh,histstep),
#         histtype='step',color='black',
#         label="$\sqrt{\kappa \\rho C}$=%5.1f $Jm^{-2} s^{-0.5} K^{-1})$"
#               % 10**(sigmed))
#plt.hist(logJ,bins=numpy.arange(histlimlow,histlimhigh,histstep),
#         histtype='step',color='black')
plt.hist(logJ, bins='fd', histtype='step', color='black', edgecolor = 'deeppink')
plt.hist(logJ_full, bins='fd', histtype='step', color='black', lw = 0.25)

plt.axvline(sig1low,color='deeppink',ls='dashed')
plt.axvline(sig1high,color='deeppink',ls='dashed')
plt.axvline(sig2low,color='deeppink',ls='dotted')
plt.axvline(sig2high,color='deeppink',ls='dotted')
# plt.axvline(numpy.log10(gamma_med),color="r",ls="solid")
plt.axvline(sigmed,color="deeppink",ls="solid")

#plt.title(name)
plt.xlabel("log $\sqrt{\kappa \\rho C} (\,J\,m^{-2}\, s^{-0.5}\, K^{-1})$", fontsize=15)
plt.ylabel("Number of Monte Carlo results", fontsize=15)
plt.xlim(1,5)
#plt.legend(loc="upper right")
plt.tight_layout()
plt.savefig(path+"/plots/gamma_hist.pdf")
plt.close()

#print("gamma=%5.2f; %5.2f to %5.2f at the 16 and 84th percentile" %
      # (gamma_med, 10**(sig1low), 10**(sig1high)))
#    (objects[name].gamma_median, sig1low, sig1high))

# EDIT below if you want to make fc histogram
#  #### Make Cratering Fraction histogram ####
#  weighted_fc,sigmed,sig1low,sig1high,sig2low,sig2high = med1sigma(fc, weights)
#  logD = [numpy.log10(d) for d in weighted_diameter]
#  print("len(D), len(weighted_diameter)", len(D), len(weighted_diameter))
#  print(sigmed,sig1low,sig1high)
#  sigmed = numpy.log10(sigmed)
#  sig1low = numpy.log10(sig1low)
#  sig1high = numpy.log10(sig1high)
#  sig2low = numpy.log10(sig2low)
#  sig2high = numpy.log10(sig2high)
#
#  weighted_fc_full,sigmed_full,sig1low_full,sig1high_full,sig2low_full,sig2high_full = med1sigma(fc_full,weights_full)
#  logD_full = [numpy.log10(d_full) for d_full in weighted_diameter_full]
#  print("len(D_full), len(weighted_diameter_full)", len(D_full), len(weighted_diameter_full))
#  print(sigmed_full,sig1low_full,sig1high_full)
#  sigmed_full = numpy.log10(sigmed_full)
#  sig1low_full = numpy.log10(sig1low_full)
#  sig1high_full = numpy.log10(sig1high_full)
#  sig2low_full = numpy.log10(sig2low_full)
#  sig2high_full = numpy.log10(sig2high_full)
#
#  if sig2high-sig2low > 0.2:
#      histstep=0.01
#  elif sig2high-sig2low > 0.02:
#      histstep=0.001
#  else:
#      histstep=0.0001
#
#  if sig2high_full-sig2low_full > 0.2:
#      histstep_full=0.01
#  elif sig2high_full-sig2low_full > 0.02:
#      histstep_full=0.001
#  else:
#      histstep_full=0.0001
#
#  plt.figure(figsize=[9,9])
#  histlimlow=int(min(logD)*1000)/1000.
#  histlimhigh=(int(max(logD)*1000)+1)/1000.
#
#  histlimlow_full=int(min(logD_full)*1000)/1000.
#  histlimhigh_full=(int(max(logD_full)*1000)+1)/1000.
#
#  #plt.hist(logD,bins=numpy.arange(histlimlow,histlimhigh,histstep),
#  #         histtype='step',color='black',label="Best-fit D = %5.3f km"
#  #            % 10**sigmed)
#  plt.hist(logD,bins=numpy.arange(histlimlow,histlimhigh,histstep),
#           histtype='step',color='black', edgecolor='#cc0000')
#  plt.hist(logD_full,bins=numpy.arange(histlimlow_full,histlimhigh_full,histstep_full),
#           histtype='step',color='black', lw = 0.25)
#  print("length of logD:",len(logD), "& length of logD_full:",len(logD_full))
#
#  plt.axvline(sig1low,color='#cc0000',ls='dashed')
#  plt.axvline(sig1high,color='#cc0000',ls='dashed')
#  plt.axvline(sigmed, color="#cc0000",ls="solid")
#  plt.axvline(sig2low,color='#cc0000',ls='dotted')
#  plt.axvline(sig2high,color='#cc0000',ls='dotted')
#  #plt.title(name)
#  plt.xlabel("log diameter (m)", fontsize=15)
#  plt.ylabel("Number of Monte Carlo results", fontsize=15)
#  #plt.legend(loc="upper right")
#  plt.tight_layout()
#  plt.savefig(path+"/plots/fc_hist.pdf")
#  plt.close()

#### Make diameter vs thermal intertia scatter plot ####
plt.figure(figsize=[9, 9])
ax=plt.gca()
plt.scatter(D,gamma,s=1.5,c='black',marker='o',edgecolor='none')
ax.set_xscale('log')
ax.set_yscale('log')

(Hbin,xbin,ybin)=numpy.histogram2d(D,gamma,bins=50)
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

#print(sig2,foo[sig2],sig1,foo[sig1])

plt.contour(xbino,ybino,Hout,colors='#cc0000',levels=[foo[sig2],foo[sig1]])
plt.text(max(D)*.7,0.7,"Contours contain 68%\nand 95.5% of all points")
#plt.text(max(D)*0.7,0.3,"Best-fit $\Gamma$ =%5.2f" % objects[name].gamma_median)
#plt.text(max(D)*0.7,0.15,"Best-fit $D$ =%5.2f km" % objects[name].diameter_median)
#plt.ylim(-2,0)
# plt.ylim(0.01,1)
plt.title(name)
plt.xlabel("Diameter (m)", fontsize=15)
plt.ylabel("Gamma", fontsize=15)
plt.tight_layout()
#plt.savefig(path+"/plots/Dvsgamma.pdf")
plt.close()

### 2D HISTOGRAM D vs ALBEDO ###

# weighted_diameter,sigmed,sig1low,sig1high,sig2low,sig2high = med1sigma(D, weights)
# weighted_diameter_full,sigmed_full,sig1low_full,sig1high_full,sig2low_full,sig2high_full = med1sigma(D_full,weights_full)
# weighted_alb,sigmed_alb,sig1low_alb,sig1high_alb,sig2low_alb,sig2high_alb = med1sigma(pV, weights)
# weighted_alb_full,sigmed_alb_full,sig1low_alb_full,sig1high_alb_full,sig2low_alb_full,sig2high_alb_full = med1sigma(pV_full, weights_full)
#
# def TwoDim(x, y, x2, y2, ylim_min, ylim_max, xlab, ylab, fignameend):
#     '''
#     This function creates a 2D hexagonal plot
#     '''
#
#     fig, ax = plt.subplots(figsize = [9, 9], tight_layout = True)
#     #l = ax.fill_between(x, y)
#     plt.hexbin(numpy.log10(x), y, cmap="gist_gray", gridsize = (60 ,30), bins = 5000)
#     plt.hexbin(numpy.log10(x2), y2, cmap="magma", gridsize=(40, 30), bins=1000, mincnt=1)
#     plt.scatter()
#     # another option for colormap: plt.nipy_spectral()
#     plt.xlabel(str(xlab), fontsize=15)
#     plt.ylabel(str(ylab), fontsize=15)
#     plt.ylim(ylim_min, ylim_max)
#     plt.colorbar(orientation="horizontal", pad = 0.1, shrink = 0.75)
#
#     #return plt.savefig("./plots/" + str(xlab) + "vs" + str(ylab) + "_2D.pdf"
#     #return plt.show()
#     return plt.savefig("./plots/DvsAlb_hex_"+ str(fignameend) + ".png", dpi = 1200)

#TwoDim(weighted_diameter_full, weighted_alb_full, weighted_diameter, weighted_alb, 0.05, 1, "log Diameter (m)", "Albedo", "patch_overlayed2")
#TwoDim(weighted_diameter, weighted_alb, 0.05, 1, "log Diameter (m)", "Albedo", "patch_weighted")
#TwoDim(D, pV, 0.05, 1, "log Diameter (m)", "Albedo", "patch_unweighted")
#TwoDim(D_full, pV_full, 0.05, 1, "log Diameter (m)", "Albedo", "full_unweighted")
#TwoDim(weighted_diameter_full, weighted_alb_full, 0.05, 1, "log Weighted Diameter (m)", "Weighted Albedo", "full_weighted")

params = Table()
params["name"] = [name for name, obj in objects.items()]
params["diameter"] = [obj.diameter_median for name, obj in objects.items()]
params["diameter_16th_percentile"] = [obj.diameter_16th_percentile for name, obj in objects.items()]
params["diameter_84th_percentile"] = [obj.diameter_84th_percentile for name, obj in objects.items()]
params["gamma"] = [obj.gamma_median for name, obj in objects.items()]
params["gamma_16th_percentile"] = [obj.gamma_16th_percentile for name, obj in objects.items()]
params["gamma_84th_percentile"] = [obj.gamma_84th_percentile for name, obj in objects.items()]
params["albedo"] = [obj.albedo_median for name, obj in objects.items()]
params["albedo_16th_percentile"] = [obj.albedo_16th_percentile for name, obj in objects.items()]
params["albedo_84th_percentile"] = [obj.albedo_84th_percentile for name, obj in objects.items()]
params["crater_fraction"] = [obj.crater_fraction_median for name, obj in objects.items()]
params["crater_fraction_16th_percentile"] = [obj.crater_fraction_16th_percentile for name, obj in objects.items()]
params["crater_fraction_84th_percentile"] = [obj.crater_fraction_84th_percentile for name, obj in objects.items()]
params["ir_fraction"] = [obj.ir_fraction_median for name, obj in objects.items()]
params["ir_fraction_16th_percentile"] = [obj.ir_fraction_16th_percentile for name, obj in objects.items()]
params["ir_fraction_84th_percentile"] = [obj.ir_fraction_84th_percentile for name, obj in objects.items()]
params["period"] = [obj.period_median for name, obj in objects.items()]
params["period_16th_percentile"] = [obj.period_16th_percentile for name, obj in objects.items()]
params["period_84th_percentile"] = [obj.period_84th_percentile for name, obj in objects.items()]
params["pole_peak_ra"] = [obj.pole_peak_ra for name, obj in objects.items()]
params["pole_peak_dec"] = [obj.pole_peak_dec for name, obj in objects.items()]
params["pole_mean_ra"] = [obj.pole_mean_ra for name, obj in objects.items()]
params["pole_mean_dec"] = [obj.pole_mean_dec for name, obj in objects.items()]
params["pole_bar"] = [obj.pole_bar for name, obj in objects.items()]
params["theta"] = [obj.theta_median for name, obj in objects.items()]
params["theta_16th_percentile"] = [obj.theta_16th_percentile for name, obj in objects.items()]
params["theta_84th_percentile"] = [obj.theta_84th_percentile for name, obj in objects.items()]

os.chdir(path)
# params.write("neo_properties.tbl", overwrite=True, format="ipac")
ofile = name+"_results.tbl"
print(ofile)
params.write(ofile, overwrite=True, format="ipac")