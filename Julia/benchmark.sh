particle_counts=(10000 100000)
thread_counts=(2 4 6 8 12 24 48)

rm times_mut.txt
for parts in "${particle_counts[@]}"
do
	for threads in "${thread_counts[@]}"
	do
		echo $parts $threads >> times_mut.txt
		for cnt in {1..7}
		do
			{ time julia --threads $threads nbodymut.jl 10 $parts ; } 2>> times_mut.txt
		done
	done
done
rm times_imut.txt
for parts in "${particle_counts[@]}"
do
	for threads in "${thread_counts[@]}"
	do
		echo $parts $threads >> times_imut.txt
		for cnt in {1..7}
		do
			{ time julia --threads $threads nbodyimut.jl 10 $parts ; } 2>> times_imut.txt
		done
	done
done
