# Example run command
# python deltaU-convergence.py -T 310.0 -k 1500.0 -m -15.0 -M 15.0 -nb 40 -x0 13.0 -f0 samples.13.0.dat -x1 14.0 -f1 samples.14.0.dat 

import sys
import numpy as np
import os
import argparse


# --------------------------------------
#     Parse                
# --------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument("-T" , dest="temperature"   , required=True)
parser.add_argument("-k" , dest="forceConstant" , required=True)
parser.add_argument("-x0", dest="x0"            , required=True)
parser.add_argument("-x1", dest="x1"            , required=True)
parser.add_argument("-f0", dest="f0"            , required=True)
parser.add_argument("-f1", dest="f1"            , required=True)
parser.add_argument("-nb", dest="Nbins"         , required=True)
parser.add_argument("-m" , dest="Umin"         , required=True)
parser.add_argument("-M" , dest="Umax"         , required=True)


args = parser.parse_args()
T         = float(args.temperature)
Kspring   = float(args.forceConstant)
x0        = float(args.x0)
x1        = float(args.x1)
FileData0 = args.f0
FileData1 = args.f1
Nbins     = int(args.Nbins)
Umin      = float(args.Umin)
Umax      = float(args.Umax)

# Parameters
Kb=0.0019872041              # in kcal/(mol.K)
beta = 1/(Kb*T)
stepU = (Umax-Umin) / float(Nbins)
axis=np.arange(Umin,Umax+1e-8,stepU)  # dU range


# Functions
def U0(x):
    U = (Kspring/2)*(x-x0)**2
    return U
def U1(x):
    U = (Kspring/2)*(x-x1)**2
    return U

def RC(dU):
    rc = dU/(Kspring *(x0-x1) ) + ( x0 +x1 ) / 2
    return rc 

def P(FileData): 
    dU = []
    for line in open(FileData):
         x = float(line.split()[1])
         dU.append( U1(x)-U0(x) )
    dU = np.array(dU,dtype=float)
   
    p, edges = np.histogram(dU,bins=axis,density=True)
    p /= np.sum(p)
    return p


# Calculate probabilities
p0 = P(FileData0)
p1 = P(FileData1)
dU = 0.5*axis[0:-1]+0.5*axis[1:]

#Transform to PMF and report
print "#       deltaU    reac.coord.        prob0        prob1       PMF-0         PMF-1        function"
for i in range(dU.size):
    if p0[i] == 0.0 or p1[i] == 0.0: 
        continue

    x = RC(dU[i])
    PMF0 = -(1/beta)*np.log(p0[i]) - U0(x) 
    PMF1 = -(1/beta)*np.log(p1[i]) - U1(x) 
    q    = np.log(p1[i]/p0[i])+beta*dU[i]
    
    print "%14.6g %14.6g %14.6g %14.6g %14.6g %14.6g %14.6g" %( dU[i], x, p0[i], p1[i], PMF0, PMF1, q) 
