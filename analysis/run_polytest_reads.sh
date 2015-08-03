dataset=$1
V=v8
LDIR=~/selection/code/files/${V}/
DATA=~/selection/counts/${V}/all.reads${dataset}.freq

while read a b
do
mkdir -p ~/selection/analysis/${V}/poly/${a}
OUT=~/selection/analysis/${V}/poly/${a}/
for what in pairs allpops
do
bsub -q reich -oo ${OUT}/${a}_${what}_reads${dataset}.out "python \
~/spindrift/Qx.py -q ${DATA} -n 10000 \
-p ${LDIR}/polypairs_${what}${dataset}.txt \
-o ${OUT}/${a}_${what}_reads${dataset} -v \
-g ${b} 2> ${OUT}/${a}_${what}_reads${dataset}.log"
done

#    what=vsSteppeMN
#    for extra in Steppe_MN Steppe_MN_no_Srubs
#    do
#    DATA=~/selection/counts/${V}/all.reads.${extra}.freq
#    python ~/spindrift/Qx.py -q ${DATA} -n 10000 \
#    -p ${LDIR}/polypairs_${what}.txt \
#    -o ${OUT}/${a}_${what}_${extra}_reads -v \
#    -g ${b} 2> ${OUT}/${a}_${what}_${extra}_reads.log
#    done

#    python ~/spindrift/Qx.py -q ${DATA} -n 10000 \
#    -p ${LDIR}/polyall.txt \
#    -o ${OUT}/${a}_reads -v \
#    -g ${b} 2> ${OUT}/${a}_all_reads.log

done < ~/selection/data/traits.txt
	
