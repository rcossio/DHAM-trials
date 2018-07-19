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
python biased-simulator.py -s samples.dat -sp sampledprob.dat -p exactprob.dat -g freeenergy.dat -k 100.0 -x0 0.0 
```

To simulate an umbrella sampling method you can run the script umbrella-sampling.sh, it will also write a .log file to be used for DHAM software. Run with,

```
bash umbrella-sampling.sh
```

Now we apply wham method to re obtain the free energy profile.
The first method uses a diagonalization routine provided by Numpy,
```
python dham.py -i umbrella-sampling.log -m -0.7 -M 0.7 -nb 70 -T 310.0 -g dham.free.dat -w dham.evals.dat -p dham.prob.dat
```

or proposing a vector and multiplying it iteratively to obtain the eigenvector with eigenvalue 1,
```
python dham-one-vector.py -i umbrella-sampling.log -m -0.7 -M 0.7 -nb 70 -T 310.0 -it 100 -g dham2.free.dat -p dham2.prob.dat
```

Aditionally, (and orthogonally to WHAM), we can check the convergence between adjacent simulations, checking the convergente between deltaU levels with the following command.

```
python deltaU-convergence.py -T 310.0 -k 1500.0 -m -15.0 -M 15.0 -nb 40 -x0 0.0 -f0 samples.15.dat -x1 0.05 -f1 samples.16.dat > tmp.dat
```
where the plots of column 1 vs 7 or 2 vs 7 should be constant in some range. 
To compare the whole possible pairs of simulations, we can run the script 

```
bash convergence.sh
```

