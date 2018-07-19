
k=1500.0
xrange=$(seq -0.7 0.05 0.7)
logfile=umbrella-sampling.log

echo "#SamplesFile  ReferenceValue  ForceConstant " > $logfile
i=1
for x0 in $xrange
do
	python biased-simulator.py -s samples.$i.dat -sp sampledprob.$i.dat -p exactprob.$i.dat -g freenergy.$i.dat -k $k -x0 $x0
	echo "samples.$i.dat      $x0       $k" >> $logfile
	((i++))
done
