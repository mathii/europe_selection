#Run this probably in the analysis directory where the results are. 

LZ=~/Packages/locuszoom/bin/locuszoom
RES=~/selection/analysis/v8/gscan/scan_results_read2.corrected.txt

GENE=MCM6
echo ${GENE}
${LZ} --metal ${RES} --pop EUR --build hg19 --source 1000G_March2012 \
 --refgene ${GENE} --pvalcol corrected.p --markercol ID --flank 1000kb


for GENE in SLC45A2 TLR1 FADS1 ATXN2 NADSYN1 HERC2 GRM5 SLC22A4
do
	echo ${GENE}
	${LZ} --metal ${RES} --pop EUR --build hg19 --source 1000G_March2012 \
	--refgene ${GENE} --pvalcol corrected.p --markercol ID --flank 500kb
done

SNP=rs1979866
echo ${SNP}
${LZ} --metal ${RES} --pop EUR --build hg19 --source 1000G_March2012 \
--refsnp ${SNP} --pvalcol corrected.p --markercol ID --flank 500kb

for GENE in EGFL8 HLA-B HLA-DPB1 HLA-DQB2 HLA-A ZKSCAN3
do
	echo ${GENE}
	${LZ} --metal ${RES} --pop EUR --build hg19 --source 1000G_March2012 \
	--refgene ${GENE} --pvalcol corrected.p --markercol ID --flank 500kb
done
