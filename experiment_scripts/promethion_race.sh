run="ERR3152367"
k=15
fastpath=datasets/promethion/classified_fastq
range=900000
r=10
n=1
input=#INSERT PATH/TO/ERR3152367.FASTQ
timing=race-${run}-b${range}-r${r}-k${k}-timing.txt
taus=0.1,0.3,0.6,1,3,6,10,30,60,100,300,600
race_path=../bin/RACEMulti
output=race-${run}-b${range}-r${r}-n1-k${k}.txt
/usr/bin/time -o ${timing} ${race_path} ${taus} ${range} ${r} ${n} ${k} ${input} ${output}
