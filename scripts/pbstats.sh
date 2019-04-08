FILE=$1
BAM=$2

INDELS=`awk '($5 == "del" || $5 == "ins")' $FILE`

Total=`samtools view $BAM | wc -l`
del_allele=`awk '($5 == "del" || $5 == "ins")' $FILE | awk '$5 == "del"' | awk '!a[$4]++' | cut -f1 | sort | uniq  | wc -l`
ins_allele=`awk '($5 == "del" || $5 == "ins")' $FILE | awk '$5 == "ins"' | awk '!a[$4]++' | cut -f1 | sort | uniq  | wc -l`
WT_reads=`awk '$3 == "WT"' $FILE | cut -f1 | sort| uniq | wc -l`
Max_del_len=`awk '($5 == "del" || $5 == "ins")' $FILE | awk '$5 == "del"' | cut -f8 | sort -nr | head -1`
Max_ins_len=`awk '($5 == "del" || $5 == "ins")' $FILE | awk '$5 == "ins"' | cut -f8 | sort -nr | head -1`
Mean_del_len=`awk '($5 == "del" || $5 == "ins")' $FILE | awk '$5 == "del"' | cut -f8 | awk '{ sum += $1 } END { if (NR > 0) print sum / NR }'`
Mean_ins_len=`awk '($5 == "del" || $5 == "ins")' $FILE | awk '$5 == "ins"' | cut -f8 | awk '{ sum += $1 } END { if (NR > 0) print sum / NR }'`
Del_reads=`awk '($5 == "del" || $5 == "ins")' $FILE | awk '$5 == "del"' | awk '!a[$4]++' | wc -l`
Ins_reads=`awk '($5 == "del" || $5 == "ins")' $FILE | awk '$5 == "ins"' | awk '!a[$4]++' | wc -l`
Uniq_small_indel_reads=`echo "$Total - $Del_reads - $Ins_reads - $WT_reads" | bc -l`

if [ -z "$Total" ]
then
    Total="0"
fi
if [ -z "$del_allele" ]
then
    del_allele="0"
fi
if [ -z "$ins_allele" ]
then
    ins_allele="0"
fi
if [ -z "$WT_reads" ]
then
    WT_reads="0"
fi
if [ -z "$Max_del_len" ]
then
    Max_del_len="0"
fi
if [ -z "$Mean_del_len" ]
then
    Mean_del_len="0"
fi
if [ -z "$Max_ins_len" ]
then
    Max_ins_len="0"
fi
if [ -z "$Mean_ins_len" ]
then
    Mean_ins_len="0"
fi
if [ -z "$Del_reads" ]
then
    Del_reads="0"
fi
if [ -z "$Ins_reads" ]
then
    Ins_reads="0"
fi
if [ -z "$Uniq_small_indel_reads" ]
then
    Uniq_small_indel_reads="0"
fi


printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" "$Total" "$WT_reads" "$Uniq_small_indel_reads" "$Del_reads" "$del_allele" "$Mean_del_len" "$Max_del_len" "$Ins_reads" "$ins_allele" "$Mean_ins_len" "$Max_ins_len"
