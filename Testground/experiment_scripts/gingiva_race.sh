run="SRR1804823"
fastpath=datasets/gingiva
k=20
input=/home/brc7/DiversitySampling/datasets/gingiva/SRR1804823_classified.fastq
timing=race-${run}-k${k}-timing.txt
race_path=/home/brc7/DiversitySampling/RACE-with-LSH-for-Squiggles/Testground/bin/RACEMulti
taus=0.001,0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,4,5
time ${race_path} ${taus} 5000000 200 1 ${k} ${input} race-${run}-k${k}-n1/race-${run}-k${k}-n1 > ${timing}

