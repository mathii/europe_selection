LDIR=~/selection/data/
DATA=~/selection/counts/all.reads.freq

while read a b
do
	mkdir -p ~/selection/analysis/poly/${a}
	OUT=~/selection/analysis/poly/${a}
	for what in vsCEU place time
	do
		python ~/spindrift/Qx.py -q ${DATA} -n 10000 \
		-p ${LDIR}/polypairs_${what}.txt \
		-o ${OUT}/${a}_${what}_reads -v \
		-g ${b} 2> ${OUT}/${a}_${what}_reads.log
	done
done < ~/selection/data/traits.txt
	
