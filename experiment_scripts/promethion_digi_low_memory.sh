digi_path=../khmer_python_env/bin/normalize-by-median.py
k=15
m=35e6
input=#INSERT PATH/TO/ERR3152367.FASTQ
for cutoff in 19
do
    output=ERR3152367_classified_k${k}_m${m}_c${cutoff}.fastq.keep
    timing_file=ERR3152367_k${k}_m${m}_c${cutoff}_timing.txt
	/usr/bin/time -o ${timing_file} python ${digi_path} -k ${k} --cutoff ${cutoff} -M ${m} -o ${output} ${input}
done
