# Notes about sampling

The simulator.py script obtaines samples of a system given a 1D free energy profile. Parameters should be modified within the script. Run with the following command, 

```
python simulator.py -s samples.dat -sp sampledprob.dat -p exactprob.dat -g freeenergy.dat
```

Interesting things to try are:
 | | 
--- | --- | ---
Harmonic well            | G = 4*X**2
Double symmetric well    | G = -4*(X*3)**2+(X*3)**4
Double asymmetric well   | G = 0.8*((X*3)-4*(X*3)**2+(X*3)**4)
