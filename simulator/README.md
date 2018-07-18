# Notes about sampling

The simulator.py script obtains samples of a system given a 1D free energy profile. Parameters should be modified within the script. Run with the following command, 

```
python simulator.py -s samples.dat -sp sampledprob.dat -p exactprob.dat -g freeenergy.dat
```

Interesting things to try are:
* Harmonic well           ``` G = 4*X**2 ```
* Double symmetric well   ``` G = -4*(X*3)**2+(X*3)**4 ```
* Double asymmetric well  ``` G = 0.8*((X*3)-4*(X*3)**2+(X*3)**4)```


The biased-script.py script samples the free energy profile plus a biasing harmonic force to be set from the command line,

```
python biased-simulator.py -s samples.dat -sp sampledprob.dat -p exactprob.dat -g freeenergyb.dat -k 100.0 -x0 0.0 
```

To simulate an umbrella sampling method you can run the script umbrella-sampling.sh, it will also write a .log file to be used for DHAM software. Run with,

```
bash umbrella-sampling.sh
```

