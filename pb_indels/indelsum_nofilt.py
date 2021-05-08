#!/usr/bin/python

import sys
import pysam
from Bio import SeqIO

#Remove reads substantially greater than amplicon size
#cat ADH1_3.fastq | paste  - - - - | awk 'length($2)  <= 900' | sed 's/\t/\n/g' > het_5_ADH1_3_less901.fastq
#map reads to reference and filter for mapping quality
#minimap2 -x map-pb -a --eqx -L -O 5,56 -E 4,1 -B 5 --secondary=no -z 400,50 -r 2k -Y -R "@RG\tID:het5_100\tSM:het5_100" ../ref/Het_5prime_100bp_250bp.fa ../raw/clean/het_5_less801_100bp.fq  | samtools sort | samtools view -q 20 -bhS > het_5_less801_100bp_q20.bam
#Get insertion read number
#python indelsum.py ADH1_3_less901_100bp_q20.bam ../ref/ADH1_3prime_100bp_250bp.fa | cut -f3- | sort | uniq -c | sort -nr | awk '$2 == "ins"' | awk '{sum += $1} END {print sum}'
#Get deletion allele number
#python indelsum.py ADH1_3_less901_100bp_q20.bam ../ref/ADH1_3prime_100bp_250bp.fa | cut -f3- | sort | uniq -c | sort -nr | awk '$2 == "del"' | wc -l


bamFile = sys.argv[1]
ref = sys.argv[2]
bamFP = pysam.Samfile(bamFile, "rb")

for seq_record in SeqIO.parse(ref, "fasta"):
    refseq = str(seq_record.seq) #store reference sequence

readcount = 0
for read in bamFP:
    mut_type = ""
    if not read.is_unmapped:   #if it's mapped
        readcount = readcount + 1
        cigarLine=read.cigar
        align = read.get_aligned_pairs() #a list of aligned read (query) and reference positions
        pos_in_cigar = 0
        seq = read.get_forward_sequence() #store aligned sequence in forward orientation
        name = read.query_name

        refaln = [] #reference alignment
        qualn = [] #query alignment

        for pospair in align:
            rpos = pospair[1]
            qpos = pospair[0]
            gap = "-"

            if rpos is None:
                refaln.append(gap)
            else:
                zrpos = rpos - 1 #zero indexed position in string
                refaln.append(refseq[zrpos])
            if qpos is None:
                qualn.append(gap)
            else:
                zqpos = qpos - 1
                qualn.append(seq[zqpos])
        str_refaln = ''.join(refaln)
        str_qualn = ''.join(qualn)
        if len(str_refaln) != len(str_refaln):
            print("Error: ref and query lengths unequal")

        for (cigarType,cigarLength) in cigarLine:
            #try:
            pos_in_cigar = pos_in_cigar + cigarLength
            if cigarType == 1 or cigarType == 2 : #insertion or deletion
                mut_type = "indel"

                if cigarLength > 0: #minimum length
                    in_start_pos = pos_in_cigar - cigarLength - 1 #query position 1bp before indel
                    refstart = align[in_start_pos][1]

                    try:
                        refend = refstart + cigarLength
                    except:
                        newalign = align[0:in_start_pos]
                        newalign = list(reversed(newalign))
                        for newpospair in newalign:
                            if newpospair[1] is not None:
                                refstart = newpospair[1]
                                break

                    #if align[in_start_pos][1] is not None:
                    #    refstart = align[in_start_pos][1]  #ref position 1bp before indel
                    #else:    
                    ##Edge case: when deletion is preceded by insertion and align[in_start_pos][1] is None
                    ##Solution: walk back from align[in_start_pos][1] in list until align[n][1] is not None
                    #    newalign = align[0:in_start_pos]
                    #    newalign = list(reversed(newalign))
                    #    for newpospair in newalign:
                    #        if newpospair[1] is not None:
                    #            refstart = newpospair[1]
                    #            break

                    refend = refstart + cigarLength

                    if cigarType == 1: #insertion
                        indel = "ins"

                    elif cigarType == 2: #deletion
                        indel = "del"

                    print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (str_refaln, str_qualn, str(readcount),name,indel,str(refstart),str(refend),str(cigarLength)))
        
    else:
        print("%s\t%s\tUN\tNA\tNA\tNA" % name)

    if mut_type is not "indel":
        print("NA\tNA\t%s\t%s\tWT\tNA\tNA\tNA" % (str(readcount),name))


