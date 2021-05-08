# Align filtered FASTQ to FASTA loci - repeat for each sample
FQ="ADH1_3_100bp.fq.maxlen.fq.gz"
REF="ADH1_3prime_100bp_250bp.fa"
OUT=${1%%.fq}
#minimap2 -x map-pb -a --eqx -L -O 5,56 -E 4,1 -B 5 --secondary=no -z 400,50 -r 2k -Y -R "@RG\tID:${FQ%%.fq}\tSM:${FQ%%.fq}" $REF $FQ  | samtools sort | samtools view -q 10 -bhS > ${OUT}_q10.bam
#samtools index ${OUT}_q10.bam
# Convert BAM to custom tabular indel format - repeat for each sample
# pysam required
python indelsum_nofilt.py ADH1_3_100bp.fq.maxlen_q10.bam ADH1_3prime_100bp_250bp.fa  >ADH1_3_100bp.fq.maxlen_q10_out_nofilt.txt
# add nick region start and end as last columns based on lookup file
for l in *nofilt.txt; do 
    name=`echo $l| sed 's/.fq.*//'| sed 's/_q10.*//'` 
    cut1=`grep "$name" cut_positions.txt| cut -f2`
    cut2=`grep "$name" cut_positions.txt| cut -f3`
    sed "s/$/\t${cut1}\t${cut2}/" $l | sed "s/^/$name\t/" > ${l%%.txt}.nickregion.txt
done
# Combine all loci tables
cat *_out_nofilt.nickregion.txt >all_loci_maxlen_q10_nickregion_master.tsv
# Output custom table filtered based on indelregion and indel length
# Also output statistics per locus to stdout
python3 analyse_master.py all_loci_maxlen_q10_nickregion_master.tsv > all_loci_maxlen_q10_nickregion_master_statistics.tsv


