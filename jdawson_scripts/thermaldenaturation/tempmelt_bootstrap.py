#!/usr/bin/env python
#Bootstrap fit for two-state temperature melt data

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from global_fitting import bootstrap

iter = 1000 #Number of iterations
Tref = 293.15 #reference temperature in Kelvin
w1_init = 28.6
w2_init = 25.5
DG_init = 0.
DS_init = 0.
DCp_init = 0.

try:
    Listfile = sys.argv[1]       #Input file name
    residue = sys.argv[2]        #Specify residue name
    pHval = sys.argv[3]          #Specify pH value
except:
    print "Usage: %s 'Input filename list Residue name'" %sys.argv[0]
    sys.exit()


def twostate(x,w1,w2,G,S,C):  #Define function: two state temperature melt 
    #DGcalc = G - S*(x-Tref) + C*(x-Tref-x*np.log(x/Tref))
    DGcalc = G - S*(x-Tref) - C*(x-Tref)**2/(2*Tref)
    A = np.exp(-DGcalc/(1.987*10**-3*x))
    return (w1+w2*A)/(1+A)

Directory = residue + "/temp_2state/"
OutputFile = Directory + residue + "_pH"+ pHval + "fit.txt"
OutputFig = Directory + residue + "_pH"+ pHval + "fit.png"
f1 = open(OutputFile,"w")

#loading data
data = np.loadtxt(Listfile,delimiter="\t")
temp = data[0]
w1H = data[2]
w15N = data[1]
w = w1H+0.154*w15N

#Bootstrap: 2state model
p = [w[0],22.6,DG_init,DS_init,DCp_init]  #initial guess of fitting parameters
global_fit = False
lmfit = False
bs = bootstrap(p,w,temp,twostate,global_fit,lmfit,iterations=iter)
#bs = bootstrap(p,y,x,Exp,global_fit,lmfit,yerr=y_err,iterations=iter)


#Median of fit parameters: 2state
#med_w1 = np.median(bs[:,0]) #using median and median absolute deviation because they are less sensitive to outliers (here, bad fits)
#std_w1 = np.median(abs(bs[:,0]-med_w1))/0.67449 #median absolute deviation converted to stdev
#med_w2 = np.median(bs[:,1])
#std_w2 = np.median(abs(bs[:,1]-med_w2))/0.67449
#med_DG0 = np.median(bs[:,2])  
#std_DG0 = np.median(abs(bs[:,2]-med_DG0))/0.67449 
#med_DS0 = np.median(bs[:,3])  
#std_DS0 = np.median(abs(bs[:,3]-med_DS0))/0.67449
#med_DCp = np.median(bs[:,4])  
#std_DCp = np.median(abs(bs[:,4]-med_DCp))/0.67449
#print
#print med_w1, std_w1
#print med_w2, std_w2
#print med_DG0, std_DG0
#print med_DS0, std_DS0
#print med_DCp, std_DCp
#y = np.zeros((len(bs),len(w)))
#for i in range (0,len(bs)):
#   p_i = [bs[i,0],bs[i,1],bs[i,2],bs[i,3], bs[i,4]]
#   y[i,:] = twostate(temp,*p_i)
#median_y = np.median(y,axis=0)
#resid = median_y - w
#chisq = np.dot(resid,resid)


#Filtering endpts: Keeping w1 and w2 within reality
#para = np.array([0.,0.,0.,0.,0.])
#for i in range(0,iter):
#    if (w1_init-4.) <= bs[i,0] <= (w1_init+4.):
#	if (w2_init-4.) <= bs[i,1] <= (w2_init+4.):
#	    a = np.array(bs[i,:])
#	    para = np.vstack([para,a])
#para = np.delete(para,0,0)
para = bs

#Median of fit parameters: 2state
med_w1 = np.median(para[:,0]) #using median and median absolute deviation because they are less sensitive to outliers (here, bad fits)
std_w1 = np.median(abs(para[:,0]-med_w1))/0.67449 #median absolute deviation converted to stdev
med_w2 = np.median(para[:,1])
std_w2 = np.median(abs(para[:,1]-med_w2))/0.67449
med_DG0 = np.median(para[:,2])  
std_DG0 = np.median(abs(para[:,2]-med_DG0))/0.67449 
med_DS0 = np.median(para[:,3])  
std_DS0 = np.median(abs(para[:,3]-med_DS0))/0.67449
med_DCp = np.median(para[:,4])  
std_DCp = np.median(abs(para[:,4]-med_DCp))/0.67449
print
print med_w1, std_w1
print med_w2, std_w2
print med_DG0, std_DG0
print med_DS0, std_DS0
print med_DCp, std_DCp
y = np.zeros((len(para),len(w)))
for i in range (0,len(para)):
   p_i = [bs[i,0],bs[i,1],bs[i,2],bs[i,3],bs[i,4]]
   y[i,:] = twostate(temp,*p_i)
median_y = np.median(y,axis=0)
resid = median_y - w
chisq = np.dot(resid,resid)

f1.write("Two state pH titration fit:\n")
f1.write("w1 = %8.6f error = %8.6f ppm\nw2 = %8.6f error = %8.6f ppm\n" %(med_w1,std_w1,med_w2,std_w2,))
f1.write("DG0 = %8.6f error = %8.6f kcal/mol\nDS0 = %8.6f error = %8.6f kcal/mol*K\n" %(med_DG0,std_DG0,med_DS0,std_DS0,))
f1.write("DCp = %8.6f error = %8.6f kcal/mol*K\n" %(med_DCp,std_DCp,))
f1.write("chi sq = %8.6f\n\n" %(chisq))

#Plotting data
p1, = plt.plot(temp,w,"ko")
p2, = plt.plot(temp,median_y,"g")
plt.xlabel('Temperature (K)')
plt.ylabel('w weighted (ppm)')
plt.title(residue+" "+pHval+' melt')
plt.legend([p1, p2], ["Data", "Two state fit"],loc='best')
plt.savefig(OutputFig)
#plt.show()

f1.close()
