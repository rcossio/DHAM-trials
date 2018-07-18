# Example run command
# python dham-one-vector.py -i umbrella-sampling.log -m -0.7 -M 0.7 -nb 70 -T 310.0 -it 10000 -g dham2.free.dat -p dham2.prob.dat 

import numpy as np
import argparse
import sys

# --------------------------------------
#     Parse                
# --------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument("-i" , dest="metadataFile"    , required=True)
parser.add_argument("-m",  dest="Xmin"            , required=True)
parser.add_argument("-M",  dest="Xmax"            , required=True)
parser.add_argument("-T" , dest="temperature"     , required=True)
parser.add_argument("-nb", dest="numberOfBins"    , required=True)
parser.add_argument("-g" , dest="freeEnergyFile"  , required=True)
parser.add_argument("-it" ,dest="iterations"      , required=True)
parser.add_argument("-p",  dest="probabilityFile" , required=True)


args = parser.parse_args()
metadataFile    = args.metadataFile 
temperature     = float(args.temperature)
numberOfBins    = int(args.numberOfBins)
freeEnergyFile  = open(args.freeEnergyFile,'w')
probabilityFile = open(args.probabilityFile,'w')
Xmin            = float(args.Xmin)
Xmax            = float(args.Xmax) 
iterations      = int(args.iterations)

#---------------------------------------
#	Parameters
#---------------------------------------
Kb           = 0.0019872041    # in kcal/mol/K
KbT          = Kb*temperature  # in kcal/mol


#----------------------------------------
#	Reading metadataFile
#----------------------------------------
print "#    Using metadataFile:", metadataFile
print "#    temperature = ",temperature," K "
print "#    Using    KbT=%14.12f"%KbT

Filename=[]
HarmonicCenter=[]
HarmonicConstant=[]

for line in open(metadataFile):
    if line.strip()[0] =="#": continue
    Filename.append(line.split()[0])
    HarmonicCenter.append(float(line.split()[1]))
    HarmonicConstant.append(float(line.split()[2]))

NumberOfSimulations = len(HarmonicCenter)
print "#    Number of simulations: ",NumberOfSimulations


#---------------------------------------------------
#	Reads all the simulations 
#---------------------------------------------------
Data=[]
SamplesPerSimulation=np.zeros(NumberOfSimulations,dtype=int)
print "#    Simulation:     Filename:     Samples: "

for k in range(NumberOfSimulations):
    xp_dumm=[]
    for line in open(Filename[k]):
	val = float(line.split()[1])
	if (val > Xmax) or (val < Xmin):
		continue
        xp_dumm.append(val)
    SamplesPerSimulation[k]= len(xp_dumm)
    Data.append(xp_dumm)       
    print "#      %6i"%(k+1),"   ",Filename[k],"   %14i"%SamplesPerSimulation[k]

# Create bins and bin centers
BinLength= (Xmax-Xmin)/float(numberOfBins)
Xmin -= BinLength
Xmax += BinLength
BinCenters = Xmin+BinLength*(0.5+np.arange(numberOfBins))
print "#    Binning Data:   Xmin=",Xmin,"   Xmax=",Xmax,"   Bin length=",BinLength 
      
# calculo de los u_i^(k) (potencial de bias-k en la caja i)
up=np.zeros((numberOfBins,NumberOfSimulations),dtype=np.float64)
for k in range(NumberOfSimulations):
    for i in range(numberOfBins):
          up[i,k]= 0.5*HarmonicConstant[k]*(BinCenters[i]-HarmonicCenter[k])**2

# Compute Population matrices and Transition matrices
# Populations and Transition are computed by asigning to each sample its nearest bin center

print "# Calculating population matrix"
PopulationMatrix          = np.zeros((numberOfBins,NumberOfSimulations),dtype=int)
NumberOfTransitionsMatrix = np.zeros((numberOfBins,numberOfBins),dtype=int)
for k in range(NumberOfSimulations):
    NearestBinCenter=[]
    for i in range(SamplesPerSimulation[k]):
        NearestBinCenter.append(np.argmin(np.abs(BinCenters-Data[k][i])))
        PopulationMatrix[NearestBinCenter[i],k] += 1

    for j in range(1,len(NearestBinCenter)):
        NumberOfTransitionsMatrix[NearestBinCenter[j],NearestBinCenter[j-1]] += 1


# calcula matriz de markov sin normalizar
print "#Calculating markov matrix"
MarkovMatrix = np.zeros((numberOfBins,numberOfBins),dtype=np.float64)
for i in range(numberOfBins):
    for j in range(numberOfBins):
        deno= 0.0
        for k in range(NumberOfSimulations):
            expo= 0.5*(up[j,k]-up[i,k])/KbT
            if PopulationMatrix[i,k] >= 1:
                deno= deno + float(PopulationMatrix[i,k])*np.exp(-expo,dtype=np.float64)
        if NumberOfTransitionsMatrix[j,i] == 0:
            MarkovMatrix[j,i]= 0.0
        else:
            MarkovMatrix[j,i]= float(NumberOfTransitionsMatrix[j,i])/deno

# Normalize the Markov matrix

print "#Normalizing matrix"
for i in range(numberOfBins):
    suma = np.sum(MarkovMatrix[:,i])
    if suma > 0.0: 
        MarkovMatrix[:,i] = MarkovMatrix[:,i]/suma
    else:
	sys.exit("Error. There is a null row of the Markov matrix. This means there are empty bins (which can be caused by the origin of the simulations).\n")

# Diagonalize the Markov matriz 
print "# Diagonalizing"
vr = np.ones(numberOfBins)
for i in range(iterations):
    vnew = MarkovMatrix.dot(vr)
    if i%50==0:
	vnew /= np.linalg.norm(vnew)
	vr   /= np.linalg.norm(vr)
        print "#",i, np.max(vnew-vr)
    vr = vnew

# Find the eigenvector correspondent with the largest eigenvalue and correct the sign if necesary (positive elements required)
if vr[0] < 0.0:
    vr= -vr
vr /= np.sum(vr)

for i in range(vr.size):
	probabilityFile.write("%14.6g %14.6g \n" %(BinCenters[i],vr[i]))
probabilityFile.close()

# Refer PFM to minimum value and print 
PMF   = -KbT*np.log(vr)
PMF  -= np.min(PMF)

freeEnergyFile.write("#  BinValue        PMF\n")
for i in range(numberOfBins):
	freeEnergyFile.write("%14.6g %14.6g \n"%(BinCenters[i],PMF[i]))
freeEnergyFile.close()

