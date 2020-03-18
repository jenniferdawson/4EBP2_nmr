#!/usr/bin/env python
#Bootstrap fit for two-state temperature melt data

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from global_fitting import bootstrap

try:
    Listfile = sys.argv[1]       #Input file name
    residue = sys.argv[2]        #Specify residue name
    tempval = sys.argv[3]          #Specify pH value
except:
    print "Usage: %s 'Input filename list Residue name'" %sys.argv[0]
    sys.exit()
Directory = residue + "/pHtitr_bytemp/"
#OutputFile = "test.txt"
OutputFile = Directory + residue + "_"+ tempval + "C_refinefit.txt"
OutputFig = Directory + residue + "_"+ tempval + "C_refinefit.png"
f1 = open(OutputFile,"w")

#constants
iter = 1000 #Number of iterations
Tref = 293.15 #reference temperature in Kelvin
pKa = 4.6
pKa_err = 0.3

#loading data
temp=float(tempval)
data_pH = np.loadtxt(Listfile,delimiter="\t")
pH = data_pH[0]-0.0028*(temp-25.) #adjusting for temperature dependence of buffer
w1H = data_pH[2]
w15N = data_pH[1]
#w = w1H+0.154*w15N
xdata = np.concatenate((pH+20.,pH),axis=0)
ydata = np.concatenate((w1H,w15N),axis=0)

def twostate(x,w1h,w1n,w2h,w2n,pKa1):  #Define function: two state temperature melt 
    wcalc = np.zeros(len(x))
    for i in range(0,len(x)):
	if x[i]>20:  #pH titration
    	    K1 = 10**-(x[i]-20.-pKa1)
	    wcalc[i] = (w1h+w2h*K1)/(1+K1)
        else:
    	    K1 = 10**-(x[i]-pKa1)
            wcalc[i] = (w1n+w2n*K1)/(1+K1)
    return wcalc

def twostate_weighted(x,w1,w2,pKa1):
    K1=10**-(x-pKa1)
    return (w1+w2*K1)/(1+K1)

#Bootstrap: 2state model
#w1H,w15N fit
p = [w1H[len(w1H)-1],w1H[0],w15N[len(w15N)-1],w15N[0],5.]  #initial guess of fitting parameters
global_fit = False
lmfit = False
bs = bootstrap(p,ydata,xdata,twostate,global_fit,lmfit,iterations=iter)
para = np.median(bs,axis=0)
y = np.zeros((iter,len(xdata)))
for i in range(0,iter):
    pi = [bs[i,0],bs[i,1],bs[i,2],bs[i,3],bs[i,4]]
    y[i,:] = twostate(xdata,*pi)
median_y = np.median(y,axis=0)
median_w1H = median_y[0:len(pH)]
median_w15N = median_y[len(pH):len(ydata)]
resid = median_y -ydata
#resid = median_y -w
chisq = np.dot(resid,resid)
pKai = pKa + pKa_err*np.random.randn(iter)
K2 = 10**(bs[:,4]-pKai)-1
E2 = np.array([0.])
for i in range(0,iter): 
    if K2[i] > 0: 
    	energy = -1.987*10**-3*293.15*np.log(K2[i]) 
    	E2 = np.vstack((E2,energy))
E2 = np.delete(E2,0,0)
median_energy = np.median(E2)

#Median of fit parameters: 2state
std_w1h = np.median(abs(bs[:,0]-para[0]))/0.67449
std_w1n = np.median(abs(bs[:,1]-para[1]))/0.67449
std_w2h = np.median(abs(bs[:,2]-para[2]))/0.67449
std_w2n = np.median(abs(bs[:,3]-para[3]))/0.67449
std_pKa = np.median(abs(bs[:,4]-para[4]))/0.67449
std_energy = np.median(abs(E2-median_energy))/0.67449
#print "w"
#print para[0], std_w1h, para[1], std_w2h
#print para[2], std_w1n, para[3], std_w2n
#print "pKa"
#print para[4], std_pKa
#print "chisq"
#print chisq
f1.write("Two state pH titration refined fit:\n")
f1.write("pKa = %8.6f error = %8.6f energy = %8.6f error = %8.6f kcal/mol\n" %(para[4],std_pKa,median_energy,std_energy))
f1.write("w1h = %8.6f error = %8.6f ppm   w1n = %8.6f error %8.6f ppm\n" %(para[0],std_w1h,para[1],std_w1n))
f1.write("w2h = %8.6f error = %8.6f ppm   w2n = %8.6f error %8.6f ppm\n" %(para[2],std_w2h,para[3],std_w2n))
f1.write("chi sq = %8.6f\n\n" %(chisq))

#Plotting data
plt.figure(1)
plt.subplot(211)
plt.plot(pH,w1H,"ko")
plt.plot(pH,median_w1H,"k-")
plt.xlabel('pH')
plt.ylabel('w1H (ppm)')
plt.gca().invert_yaxis()
#plt.gca().invert_xaxis()
plt.title(residue+" "+tempval+'C pH titration')
plt.subplot(212)
p1,=plt.plot(pH,w15N,"ko")
p2,=plt.plot(pH,median_w15N,"k-")
plt.xlabel('pH')
plt.ylabel('w15N (ppm)')
plt.gca().invert_yaxis()
#plt.gca().invert_xaxis()
plt.legend([p1, p2], ["Data", "Two state fit"],loc='best')
plt.savefig(OutputFig)
#plt.show()

f1.close()
