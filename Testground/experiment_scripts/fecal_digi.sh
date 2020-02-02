digi_path=/home/brc7/tools/khmer_python_env/bin/normalize-by-median.py
k=12
m=3e8
input=/home/brc7/DiversitySampling/datasets/fecal/SRR2175724_classified.fastq

for cutoff in 1 2 3 5 9 13 19 25
do
    output=SRR2175724_classified_k${k}_m${m}_c${cutoff}.fastq.keep
    timing_file=SRR2175724_k${k}_m${m}_c${cutoff}_timing.txt
    
    /usr/bin/time -o ${timing_file} python ${digi_path} -k ${k} --cutoff ${cutoff} -M ${m} -o ${output} ${input}
done
