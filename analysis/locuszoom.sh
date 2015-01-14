#Run this probably in the analysis directory where the results are. 

LZ=~/Packages/locuszoom/bin/locuszoom
RES=~/selection/analysis/gscan/scan_results.txt

for GENE in MCM6 SLC45A2 NADSYN1 HERC2
do
	echo ${GENE}
	${LZ} --metal ${RES} --pop EUR --build hg19 --source 1000G_March2012 \
	--refgene ${GENE} --pvalcol corrected.p --markercol ID --flank 500kb &
done
