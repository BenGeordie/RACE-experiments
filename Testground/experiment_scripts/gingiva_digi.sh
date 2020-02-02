digi_path=/home/brc7/tools/khmer_python_env/bin/normalize-by-median.py
k=20
m=3e8
input=SRR1804823_classified.fastq
for cutoff in 1 2 3 9 15
do
    output=SRR1804823_classified_k${k}_m${m}_c${cutoff}.fastq.keep
    timing_file=SRR1804823_k${k}_m${m}_c${cutoff}_timing.txt
    
    /usr/bin/time -o ${timing_file} python ${digi_path} -k ${k} --cutoff ${cutoff} -M ${m} -o ${output} ${input}
done
