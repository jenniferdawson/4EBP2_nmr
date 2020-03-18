#!/usr/bin/env python
#Bootstrap fit for two-state temperature melt data

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy

from scipy.optimize import curve_fit
from global_fitting import bootstrap

Tref = 293.15 #reference temperature in Kelvin

try:
    Listfile = sys.argv[1]       #Input file name 
    residue = sys.argv[2]        #Specify residue name
    pHval =sys.argv[3]		 #Specify pH
except:
    print "Usage: %s 'Input filename list Residue name'" %sys.argv[0]
    sys.exit()

Directory = ""
OutputFile = "test.txt"
#OutputFile = Directory + residue + pHval + "fit.txt"
#OutputFig = Directory + residue + pHval + "fit.png"
f1 = open(OutputFile,"w")

#loading data
data = np.loadtxt(Listfile,delimiter="\t")
temp = data[0]
w1H = data[2]
w15N = data[1]
w = w1H+0.154*w15N
pKa = 5.09
pKa_err = 0.08
pH = float(pHval) - 0.0028*(temp - 298.15) #temperature dependence value from BioBuffer info, referenced to 25C

iter = 10 #number of fits
iter_boot = 10 #number of bootstrap rounds

y = np.zeros((iter,len(w)))
para_mean = np.zeros((iter,6))
pKai = pKa + pKa_err*np.random.randn(iter)
print pKai
for i in range(0,iter):
    pKa_set = pKai[i]
    def threestate(x,w1,w2,w3,G,S,C):  #Define function: three state pHtitration, then temperature melt 
	K1 = 10**-(pH-pKa_set)
    	DGcalc = G - S*(x-Tref) + C*(x-Tref-x*np.log(x/Tref))
    	K2 = np.exp(-DGcalc/(1.987*10**-3*x))
    	return (w1+w2*K1+w3*K1*K2)/(1+K1+K1*K2)

    #Bootstrap
    p = [w[0],w[0],w[len(w)-1],2.,0.,0.]  #initial guess of fitting parameters
    global_fit = False
    lmfit = False
    para = bootstrap(p,w,temp,threestate,global_fit,lmfit,iterations=iter_boot)
    para_mean[i,:] = np.mean(para,axis=0)
    print "Blah"
    print para_mean[i,:]
    y[i,:] = threestate(temp,*para_mean) 

median_y = np.median(y,axis=0)  
median_para = np.median(para_mean,axis=0)
std_w1 = np.median(abs(para[:,0]-median_para[0]))/0.67449 #median absolute deviation converted to stdev
std_w2 = np.median(abs(para[:,1]-median_para[1]))/0.67449
std_w3 = np.median(abs(para[:,2]-median_para[2]))/0.67449
std_DG0 = np.median(abs(para[:,3]-median_para[3]))/0.67449
std_DS0 = np.median(abs(para[:,4]-median_para[4]))/0.67449
std_DCp = np.median(abs(para[:,5]-median_para[5]))/0.67449

print
print median_para,std_w1,std_w2,std_DG0,std_DS0
print median_DCp, std_DCp

resid = median_y -w
chisq = np.dot(resid,resid)
print chisq

#f1.write("%s %s Two state thermal denaturation fit:\n" %(residue,pHval))
#f1.write("w1 = %8.6f error = %8.6f ppm\nw2 = %8.6f error = %8.6f ppm\n" %(median_w[0],std_w1,median_w[1],std_w2,))
#f1.write("DG0 = %8.6f error = %8.6f kcal/mol\nDS0 = %8.6f error = %8.6f kcal/mol*K\n" %(median_para[0],std_DG0,median_para[1],std_DS0,))
#f1.write("DCp = %8.6f error = %8.6f kcal/mol*K\n" %(median_para[2],std_DCp,))
#f1.write("chi sq = %8.6f\n\n" %(chisq))
#f1.write("\n\n#Data and fit\nTemp (deg C)\tw (ppm)\tw fit (ppm)\n")
#for i in range(0,len(w)):
#    f1.write("%5.2f\t%8.6f\t%8.6f\n" %(temp[i]-273.15,w[i],median_y[i]))

#Plotting data
p1, = plt.plot(temp,w,"ko")
p2, = plt.plot(temp,median_y,"g")
plt.xlabel('Temperature (K)')
plt.ylabel('w weighted (ppm)')
plt.title(residue+" "+pHval+' melt')
plt.legend([p1, p2], ["Data", "Two state fit"],loc='best')
#plt.savefig(OutputFig)
plt.show()


f1.close()
