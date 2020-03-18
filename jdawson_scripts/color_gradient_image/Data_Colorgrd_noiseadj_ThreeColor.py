#!/usr/bin/env python
#
#Data_Colorgrd.py: 1. Makes string (for use in pymol selection) specifying the analyzed residues in input file
#              2. Outputs command setting color of each analyzed residue in a color gradient based on data values
# 
#JED June 13,2011

import os
import sys

outputFile = "./AnalyzedResidues_pHdepE_test.pml"

#Color settings=[R,G,B].  Gradient from color1 to color2
color1 = [0.02 , 0.50 , 0.72]   #blue
color2 =  [1.00 , 1.00 , 1.00]   #white
color3 = [1.00 , 0.00 , 0.00]   #red
try:
    DataFile = sys.argv[1]       #Specify input data file
except:
    print "Usage: %s 'Input file'" %sys.argv[0]
    sys.exit()

outputFile = open(outputFile, "w")
data = []
Residue = []
select = "select AnalyzedResidues, resi "
if (os.path.isfile(DataFile)):
    f = open(DataFile, "r")
    for line in (f.readlines()):
	lineStriped = line.strip()
	if (lineStriped[0] != '#' and lineStriped[0].isdigit() != 0):
            line = lineStriped.split()
            name = int(line[0])
            select = select + str(name) + "+"
            data = data[:] + [float(line[1])]
            Residue = Residue[:] + [name]      
    f.close()
select = select.rstrip("+")
outputFile.write("%s\n" %(select))
outputFile.write("select N, name n\n")
outputFile.write("select AnalyzedResidues, AnalyzedResidues&N\n")
outputFile.write("show spheres, AnalyzedResidues\n")
#MaxData = max(data)
MaxData=4.0
#MinData = min(data)
MinData=-4.0

print data[0]
for x in range(len(data)):
    if data[x]<= 0.:
    	Diff = -1*MinData
    	if data[x]<MinData:
            scale=0.00
        else:
            scale = (data[x]-MinData)/Diff
        c1 = color3[0] + (color2[0]-color3[0])*scale
        c2 = color3[1] + (color2[1]-color3[1])*scale
        c3 = color3[2] + (color2[2]-color3[2])*scale
    else:
    	Diff = MaxData
    	if data[x]>MaxData:
            scale=0.00
        else:
            scale = (MaxData-data[x])/Diff
        c1 = color1[0] + (color2[0]-color1[0])*scale
        c2 = color1[1] + (color2[1]-color1[1])*scale
        c3 = color1[2] + (color2[2]-color1[2])*scale
    #outputFile.write("select Res%s, resi %s\n" %(Residue[x],Residue[x]))
    #outputFile.write("select Res%s=Res%s&N\n" %(Residue[x],Residue[x]))
    outputFile.write("set_color Res%s= [%1.3f , %1.3f , %1.3f]\n" %(Residue[x],c1,c2,c3))
    outputFile.write("color Res%s, resi %s\n" %(Residue[x],Residue[x]))

outputFile.close()

#[R,G,B] values for different colors in pymol (From pymol manual)
#optimized rgb values for cmyk output:
#dblue= [0.05 , 0.19 , 0.57]
#blue=  [0.02 , 0.50 , 0.72]
#mblue= [0.5  , 0.7  , 0.9 ]
#lblue= [0.86 , 1.00 , 1.00]
 
#green= [0.00 , 0.53 , 0.22]
#lgreen=[0.50 , 0.78 , 0.50]
#yellow=[0.95 , 0.78 , 0.00]
#orange=[1.00 , 0.40 , 0.0 ]
 
# these are trivial
#red=   [1.00 , 0.00 , 0.00]
#mred=  [1.00 , 0.40 , 0.40]
#lred=  [1.00 , 0.80 , 0.80]
#vlred= [1.00 , 0.90 , 0.90]
#white= [1.00 , 1.00 , 1.00]
#vlgray=[0.95 , 0.95 , 0.95]
#lgray= [0.90 , 0.90 , 0.90]
#gray=  [0.70 , 0.70 , 0.70]
#dgray= [0.50 , 0.50 , 0.50]
#dgray=[0.30 , 0.30 , 0.30]
#black= [0.00 , 0.00 , 0.00]
##
