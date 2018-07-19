
k=1500.0
xrange=$(seq -0.70 0.05 0.65)

i=1
for x0 in $xrange
do
	x1=$(python -c "print ${x0}+0.05")
	j=$(($i+1))
	python deltaU-convergence.py -T 310.0 -k $k -m -15.0 -M 15.0 -nb 40 -x0 $x0 -f0 samples.$i.dat -x1 $x1 -f1 samples.$j.dat > dUconv.$i.dat
	((i++))
done
