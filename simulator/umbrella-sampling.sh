
k=100.0
xrange=$(seq -0.5 0.05 0.5)
logfile=umbrella-sampling.log

echo "#SamplesFile  ForceConstant  ReferenceValue" > $logfile
i=1
for x0 in $xrange
do
	python biased-simulator.py -s samples.$i.dat -sp sampledprob.$i.dat -p exactprob.$i.dat -g freeenergy.$i.dat -k $k -x0 $x0
	echo "samples.$i.dat      $k      $x0" >> $logfile
	((i++))
done
