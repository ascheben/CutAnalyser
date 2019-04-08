#usage getstats.sh [min_deletion_len_for_totals_count] [min_homology_len] [treatments_file] [mhfinder_out_table]
#TotalInsertion Alleles/Reads include the insertions smaller than user-defined minimum length

mindel=$1 #only affects counts of total deletion reads and total deletion alleles
minmh=$2 
treatments=$3
infile=$4

{
if [[ "$mindel" -gt "$minmh" ]]
then
    printf 'Warning: Minimum deletion length is greater than minimum microhomology length.\nMicrohomologies associated with deletions shorter than the minimum deletion length will be excluded!\n'
fi
}

printf 'Group\tMH_Read_Ratio\tMH_Allele_Ratio\tPTD_Read_Ratio\tPTD_Allele_Ratio\tReads\tAlleles\tSingle_deletion_MH_reads\tMult_deletion_MH_reads\tSingle_Del_MH_Alleles\tMul_Del_MH_Alleles\tReads_with_deletion>min\tAlleles_with_deletion>min\tTandem_Dup\tTandem_Dup_Alleles\tTandem_Dup_Reads\tIns_Alleles>min\tIns_Reads>min\tTotal_Ins_Alleles\tTotal_Ins_Reads\n'
while read group; do 

    
    z=`awk -v i=$group '$1==i' $infile | awk '!a[$2]++' | awk '{sum+=$6}END{print sum}'`
    y=`awk -v i=$group '$1==i' $infile | awk '!a[$2]++'| wc -l`


	# "Single and mult MH read count" Sum of reads without mismatch, with MH and with deletion length and MH length >= minimums
    #a=`awk -v i=$group '$1==i' $infile |sed 's/\./0/g' | awk -v minh=$minmh -v mind=$mindel '($16=="N" && $9>=mind && $13>=minh)' | awk '{sum+=$6}END{print sum}'` 

	# "Single deletion MH read count" 
    b=`awk -v i=$group '$1==i' $infile |sed 's/\./0/g' | awk -v minh=$minmh -v mind=$mindel '($16=="N" && $9>=mind && $13>=minh && $2==$17)' | awk '{sum+=$6}END{print sum}'`

	# "Mult deletion MH read count uniquified"
    c=`awk -v i=$group '$1==i' $infile |sed 's/\./0/g' | awk -v minh=$minmh -v mind=$mindel '($16=="N" && $9>=mind && $13>=minh && $2!=$17)' | awk '!a[$2]++'| awk '{sum+=$6}END{print sum}'`
	
	# "Single deletion MH alleles"
    d=`awk -v i=$group '$1==i' $infile |sed 's/\./0/g' | awk -v minh=$minmh -v mind=$mindel '($16=="N" && $9>=mind && $13>=minh && $2==$17)' | wc -l`
	
	# "Uniquified mult deletion MH Alleles"
    e=`awk -v i=$group '$1==i' $infile |sed 's/\./0/g' | awk -v minh=$minmh -v mind=$mindel '($16=="N" && $9>=mind && $13>=minh && $2!=$17)' | awk '!a[$2]++'| wc -l`
	
	# "Reads with deletions (mult only counted once)"
    f=`awk -v i=$group -v mind=$mindel '($1==i && $7=="del" && $9>=mind)' $infile | awk '!a[$2]++'| awk '{sum+=$6}END{print sum}'`
	
    # "Alleles with deletions (mult only counted once)"
    g=`awk -v i=$group -v mind=$mindel '($1==i && $7=="del" && $9>=mind)' $infile | awk '!a[$2]++'| wc -l`

    # "Perfect tandem duplications (counting multiple tandem duplication in one allele)"
    h=`awk -v i=$group '($1==i && $18!="NA")' $infile | awk '{sum+=$18} END{print sum;}'`
 
    # "Alleles with at least one perfect tandem duplication"
    i=`awk -v i=$group '$1==i' $infile | awk '($18!="NA" && $18!="0")' | wc -l`

    # Reads with at least one tandem duplication reads
    j=`awk -v i=$group '$1==i' $infile | awk '$18!="NA"' | awk '$18!="0"' | awk '{sum+=$6}END{print sum}'`

    #Alleles with an insert >= minimum insertion length set for mhfinder.py
    k=`awk -v i=$group '($1==i && $18!="NA")' $infile  | wc -l`

    # Reads with insert >= minimum insertion length set for mhfinder.py
    l=`awk -v i=$group '$1==i' $infile | awk '$18!="NA"' | awk '{sum+=$6}END{print sum}'`
   
    # Total insertion alleles
    m=`awk -v i=$group '$1==i' $infile | awk '$7=="Ins"' | wc -l`

    # Total insertion reads
    n=`awk -v i=$group '$1==i' $infile | awk '$7=="Ins"' | awk '{sum+=$6}END{print sum}'`


#    if [ -z "$a" ]
#    then
#        a="0"
#    fi
    if [ -z "$b" ]
    then
        b="0"
    fi

    if [ -z "$c" ]
    then
        c="0"
    fi
    if [ -z "$d" ]
    then
        d="0"
    fi
    if [ -z "$e" ]
    then
        e="0"
    fi
    if [ -z "$f" ]
    then
        f="0"
    fi
    if [ -z "$g" ]
    then
        g="0"
    fi
    if [ -z "$h" ]
    then
        h="0"
    fi
    if [ -z "$i" ]
    then
        i="0"
    fi
    if [ -z "$j" ]
    then
        j="0"
    fi
    if [ -z "$k" ]
    then
        k="0"
    fi
    if [ -z "$l" ]
    then
        l="0"
    fi
    if [ -z "$m" ]
    then
        m="0"
    fi
    if [ -z "$n" ]
    then
        n="0"
    fi
    if [ -z "$o" ]
    then
        o="0"
    fi
    if [ -z "$p" ]
    then
        p="0"
    fi
    if [ -z "$q" ]
    then
        q="0"
    fi
    if [ -z "$u" ]
    then
        u="0"
    fi

    #UniqueMH_Del_Reads/Total_Del_Reads 
    if [ "$f" -ne 0 ]
    then
        o=`echo "($b + $c) / $f" | bc -l`
    else
        o=0
    fi
    #UniqMH_Alleles/TotalAlleles
    if [ "$g" -ne 0 ]
    then
        p=`echo "($d + $e) / $g " | bc -l`
    else
        p=0
    fi
    #PTD_Reads/Insertions_reads
    if [ "$l" -ne 0 ]
    then
        q=`echo "$j / $l" | bc -l`
    else
        q=0
    fi
    #PTD_Alleles/Insertions_alleles
    if [ "$k" -ne 0 ]
    then
        u=`echo "$i / $k" | bc -l`
    else
        u=0
    fi

    printf '%s\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$group" "$o" "$p" "$q" "$u" "$z" "$y" "$b" "$c" "$d" "$e" "$f" "$g" "$h" "$i" "$j" "$k" "$l" "$m" "$n"
#When variable assignments fail due to missing content in input table, variable values can carry over between steps in the loop
#A quick hack to avoid this is to set the values to zero at the end of the loop
    z=0
    y=0
    b=0
    c=0
    d=0
    e=0
    f=0
    g=0
    h=0
    i=0
    j=0
    k=0
    l=0
    m=0
    n=0
    o=0
    p=0
    q=0

done<$treatments
