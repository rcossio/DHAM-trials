# Example run command
# python biased-simulator.py -s samples.dat -sp sampledprob.dat -p exactprob.dat -g freenergy.dat -k 100.0 -x0 0.0 

import numpy as np
import argparse


# --------------------------------------
#     Parse                
# --------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument("-s" , dest="samples"           , required=True)
parser.add_argument("-sp", dest="sampledProbability", required=True)
parser.add_argument("-p" , dest="exactProbability"  , required=True)
parser.add_argument("-g" , dest="freeEnergy"        , required=True)
parser.add_argument("-k" , dest="forceConstant"     , required=True)
parser.add_argument("-x0", dest="referenceValue"    , required=True)

args = parser.parse_args()
samplesFile            = open(args.samples,'w')
sampledProbabilityFile = open(args.sampledProbability,'w')
exactProbabilityFile   = open(args.exactProbability,'w')
freeEnergyFile         = open(args.freeEnergy,'w')
k                      = float(args.forceConstant)
x0                     = float(args.referenceValue)

#--------------------------------------------------
T        = 310.0        #in K
Kb       = 0.0019872041 #in kcal/mol/K
beta     = 1/(Kb*T)       #in kcal/mol
Nsamples = 1000000

dX = 1e-4
X  = np.arange(-1.0,1.0+1e-10,dX)
bias = (k/2.)*(X-x0)**2
G = 0.8*((X*3)-4*(X*3)**2+(X*3)**4)
G += bias
G -= np.min(G)
P  = np.exp(-beta*G)
P /= np.sum(P)

s = np.random.choice(X,size=Nsamples,p=P)

h,e = np.histogram(s,bins=100,density=True)
e  = (e[0:-1]+e[1:])/2

for i in range(G.size):
        freeEnergyFile.write("%14.6g %14.6g\n" %(X[i],G[i]))
freeEnergyFile.close()

for i in range(P.size):
        exactProbabilityFile.write("%14.6g %14.6g\n" %(X[i],P[i]/dX))
exactProbabilityFile.close()

for i in range(s.size):
	samplesFile.write("%14.6g %14.6g\n" %(i,s[i]))
samplesFile.close()

for i in range(h.size):
        sampledProbabilityFile.write("%14.6g %14.6g\n" %(e[i],h[i]))
sampledProbabilityFile.close()
