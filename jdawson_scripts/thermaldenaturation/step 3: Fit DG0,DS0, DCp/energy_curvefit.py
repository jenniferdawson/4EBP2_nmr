#!/usr/bin/env python
#Bootstrap fit 

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy

from scipy.optimize import curve_fit
from global_fitting import bootstrap

Tref = 293.15 #reference temperature in Kelvin

try:
    Listfile = sys.argv[1]       #Input file name for pKa values
    residue = sys.argv[2]        #Specify residue name
except:
    print "Usage: %s 'Input filename list Residue name'" %sys.argv[0]
    sys.exit()

Tref = 293.15 #Reference temperature in K

def twostate(t,DG0,DS0,DCp):  #energy for two state thermal denaturation
    return ( DG0 - DS0*(t-Tref) + DCp*(t-Tref-t*np.log(t/Tref)) )


Directory = residue + "/pHtitr_bytemp_curvefit/"
OutputFile = Directory + residue + "_twostate_thermalfit_refined.txt"
OutputFig = Directory + residue + "_twostate_thermalfit_refined.png"
#OutputFile = "test.txt"
f1 = open(OutputFile,"w")

#loading data
Tref = 293.15 #reference temp in kelvin
data = np.loadtxt(Listfile,delimiter="\t")
temp = data[:,0] + 273.15
DG = data[:,1]
DG_err = data[:,2]

def thermal_energy(T,G,S,C):
    return (G-S*(T-Tref)+C*(T-Tref-T*np.log(T/Tref)))

#curvefit
iter = 1000  #number of variations 
para = np.zeros((iter,3))
yfit = np.zeros((iter,len(DG)))
for i in range(0,iter):
    DGi = DG + np.multiply(DG_err,np.random.randn(len(DG)))
    p = [0.5,0.,0.]
    popt, pcov = curve_fit(thermal_energy,temp,DGi,p)
    para[i,:] = popt
    yfit[i,:] = thermal_energy(temp,*popt)
med_para = np.median(para,axis=0)
med_y = np.median(yfit,axis=0)
std_DG0 = np.median(abs(para[:,0]-med_para[0]))/0.67449
std_DS0 = np.median(abs(para[:,1]-med_para[1]))/0.67449
std_DCp = np.median(abs(para[:,2]-med_para[2]))/0.67449
print med_para[0],std_DG0
print med_para[1],std_DS0
print med_para[2],std_DCp
print med_y
resid = med_y - DG
chisq = np.dot(resid,resid)

f1.write("Thermal denaturation parameters:\n")
f1.write("DG0 = %8.6f error = %8.6f kcal/mol\n" %(med_para[0],std_DG0))
f1.write("DS0 = %8.6f error = %8.6f kcal/mol\n" %(med_para[1],std_DS0))
f1.write("DCp = %8.6f error = %8.6f kcal/mol\n" %(med_para[2],std_DCp))
f1.write("chi sq = %8.6f\n\n" %(chisq))

line = np.arange(0,temp[len(temp)-1]-273.15+10,5)
yline = np.zeros(len(line))
plt.errorbar(temp-273.15,DG,yerr=DG_err,fmt='ko')
plt.plot(temp-273.15,med_y,'k-')
plt.plot(line,yline,'k-')
plt.xlabel("Temperature")
plt.ylabel("DG unfolding kcal/mol")
plt.title("Thermal denaturation of G48")
plt.savefig(OutputFig)
#plt.show()
f1.close()
