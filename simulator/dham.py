#python DHAM.py $min $max $bins $temperature xxx.wham.inp > wham.$m.out

# Libraries
import numpy as np
import sys
import math

# Parameters
Xmin         = float(sys.argv[1])
Xmax         = float(sys.argv[2]) 
NumberOfBins = int(sys.argv[3])
Temperature  = float(sys.argv[4]) 
metadatafile = sys.argv[5]
Iterations   = 2000
Kb           = 0.0019872041 # in kcal/mol/K
KbT          = Kb*Temperature

#Reading metadatafile

print "#    Using metadatafile:", metadatafile
print "#    Temperature = ",Temperature," K "
print "#    Using    KbT=%14.12f"%KbT

Filename=[]
HarmonicCenter=[]
HarmonicConstant=[]

for line in open(metadatafile):
    if line.strip()[0] =="#": continue
    Filename.append(line.split()[0])
    HarmonicCenter.append(float(line.split()[1]))
    HarmonicConstant.append(float(line.split()[2]))

NumberOfSimulations = len(HarmonicCenter)
print "#    Number of simulations: ",NumberOfSimulations


# Reads all the simulations 
Data=[]
SamplesPerSimulation=np.zeros(NumberOfSimulations,dtype=int)
print "#    Simulation:     Filename:     Samples: "

for k in range(NumberOfSimulations):
    xp_dumm=[]
    for line in open(Filename[k]):
        xp_dumm.append(float(line.split()[1]))
    SamplesPerSimulation[k]= len(xp_dumm)
    Data.append(xp_dumm)       
    print "#      %6i"%(k+1),"   ",Filename[k],"   %14i"%SamplesPerSimulation[k]

# Create bins and bin centers
BinLength= (Xmax-Xmin)/float(NumberOfBins)
Xmin -= BinLength
Xmax += BinLength
NumberOfBins +=2
BinCenters = Xmin+BinLength*(0.5+np.arange(NumberOfBins))
print "#    Binning Data:   Xmin=",Xmin,"   Xmax=",Xmax,"   Bin length=",BinLength 
      
# calculo de los u_i^(k) (potencial de bias-k en la caja i)
up=np.zeros((NumberOfBins,NumberOfSimulations),dtype=np.float64)
for k in range(NumberOfSimulations):
    for i in range(NumberOfBins):
          up[i,k]= 0.5*HarmonicConstant[k]*(BinCenters[i]-HarmonicCenter[k])**2

# Compute Population matrices and Transition matrices
# Populations and Transition are computed by asigning to each sample its nearest bin center

print "# Calculating population matrix"
PopulationMatrix          = np.zeros((NumberOfBins,NumberOfSimulations),dtype=int)
NumberOfTransitionsMatrix = np.zeros((NumberOfBins,NumberOfBins),dtype=int)
for k in range(NumberOfSimulations):
    NearestBinCenter=[]
    for i in range(SamplesPerSimulation[k]):
        NearestBinCenter.append(np.argmin(np.abs(BinCenters-Data[k][i])))
        PopulationMatrix[NearestBinCenter[i],k] += 1

    for j in range(1,len(NearestBinCenter)):
        NumberOfTransitionsMatrix[NearestBinCenter[j],NearestBinCenter[j-1]] += 1


# calcula matriz de markov sin normalizar
print "#Calculating markov matrix"
MarkovMatrix = np.zeros((NumberOfBins,NumberOfBins),dtype=np.float64)
for i in range(NumberOfBins):
    print "#... row "+str(i+1)+"/"+str(NumberOfBins)
    for j in range(NumberOfBins):
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
for i in range(NumberOfBins):
    suma = np.sum(MarkovMatrix[:,i])
    if suma > 0.0: 
        MarkovMatrix[:,i] = MarkovMatrix[:,i]/suma
    else:
        MarkovMatrix[:,i] = np.zeros(nbin,dtype=np.float64)

# Diagonalize the Markov matriz 
print "# Diagonalizing"
vr = np.random.rand(NumberOfBins)
for i in range(Iterations):
    vnew = MarkovMatrix.dot(vr)
    if i%50==0: 
        print "#",i, np.max(vnew-v)
    vr = vnew/np.linalg.norm(vnew)

# Find the eigenvector correspondent with the largest eigenvalue and correct the sign if necesary (positive elements required)
if vr[0] < 0.0:
    vr= -vr

# Refer PFM to minimum value and print 
PMF1 = -KbT*np.log(vr)

#Results
print "#BinNumber       BinValue        DhamPMF"
for j in range(1,NumberOfBins-1):
  print "%14.8f %14.8f"%(BinCenters[j],PMF1[j]-min(PMF1))

