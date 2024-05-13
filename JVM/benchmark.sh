particle_counts=(10000 100000)
thread_counts=(1 3 5 7 11 23 47)

rm times.txt
for parts in "${particle_counts[@]}"
do
	for threads in "${thread_counts[@]}"
	do
		echo $parts $threads >> times.txt
		for cnt in {1..7}
		do
			{ time java -Djava.util.concurrent.ForkJoinPool.common.parallelism=$threads nbody 10 $parts ; } 2>> times.txt
		done
	done
done
