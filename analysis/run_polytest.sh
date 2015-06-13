LDIR=~/selection/data/
DATA=~/data/v6/use/v61kg_polynames

while read a b
do
mkdir -p ~/selection/analysis/poly/${a}
OUT=~/selection/analysis/poly/${a}
for what in vsCEU place time
do
python ~/spindrift/Qx.py -d ${DATA} -n 10000 \
-p ${LDIR}/polypairs_${what}.txt \
-o ${OUT}/${a}_${what} -v \
-g ${b} -i ${LDIR}/polypops_inbred.txt \ 
2> ${OUT}/${a}_${what}.log
done
done < ~/selection/data/traits.txt
	
