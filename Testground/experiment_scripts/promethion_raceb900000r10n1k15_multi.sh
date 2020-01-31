config="b900000r10n1k15"
run="ERR3152367"
k=15
fastpath=datasets/promethion/classified_fastq
range=900000
r=10
n=1
input=${fastpath}/${run}_classified.fastq
timing=race-${run}-b${range}-r${r}-k${k}-timing.txt
taus=0.1,0.3,0.6,1,3,6,10,30,60,100,300,600
race_path=/home/brc7/DiversitySampling/RACE-with-LSH-for-Squiggles/Testground/bin/RACEMulti
output=datasets/results/results_ERR315/${config}/race-${run}-b${range}-r${r}-n1-k${k}.txt
/usr/bin/time -o ${timing} ${race_path} ${taus} ${range} ${r} ${n} ${k} ${input} ${output}
