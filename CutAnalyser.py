#!/usr/bin/env python3
import sys
import csv
import re
import os
import string
import random
from collections import Counter
from argparse import ArgumentParser
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO

###########################################
# Authors: Armin Scheben & Felix Wolter   #
# Date: 2018/08/18                        #
###########################################

def repeat_to_length(string_to_expand, length):
    '''Expand a string n times.abs'''
    return (string_to_expand * (int(length/len(string_to_expand))+1))[:length]

def random_string(string_length):
    '''Generate a random string of fixed length '''
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(string_length))

def gaps_to_bases(gaps):
    '''Convert list of gap positions to
    list of base positions.
    '''
    bases = []
    for i in range(len(gaps) - 1):
        # get zero based positions not including last position
        # example: AGGTCT[2:4] is GT  
        start = gaps[i][1]
        stop = gaps[i+1][0]
        bases.append([start,stop])
    return bases

def realign(seq1, seq2):
    '''Globally realign sequences using custom alignment
    penalties. Returns list of two realigned gapped sequences.
    '''
    # Strip gap characters from sequence
    cseq1 = seq1.replace('-', '')
    cseq2 = seq2.replace('-', '')
    # Strong penalise mismatches and reduce match reward
    alignments = pairwise2.align.globalms(cseq1, cseq2, 0, -6, -4, -0.05, one_alignment_only=True,penalize_extend_when_opening=True)
    #print("Realignments\n", format_alignment(*alignments[0]))
    rseqs = [alignments[0][0], alignments[0][1]]
    return rseqs

def correct_realign(seq1, seq2):
    '''Make alignment corrections of errors caused
    by custom penalties.
    '''
    # Get insertion site(s) from reference
    g1 = getgaps(seq1)
    # Get deletion site(s) from mutant
    g2 = getgaps(seq2)
    # If there is only one gap in reference and mutant
    if len(g1) == 1 and len(g2) == 1:
        # If insertion starts before deletion
        if g1[0][0] < g2[0][0]:
            diff = g2[0][0] - g1[0][1]
        # if insertion starts after deletion
        elif g1[0][0] > g2[0][0]:
            diff = g1[0][0] - g2[0][1]

        # Index of insertion start in reference
        sg1 = g1[0][1]
        # Index of deletion end in mutatn
        sg2 = g2[0][0]
        #diff = sg2 - sg1

        if diff > 1 and diff < 4:
            for i in range(diff - 1):
                # Insert gap at beginning of insertion
                seq1 = seq1[:sg1] + '-' + seq1[sg1:]
                # Insert gap at end of deletion
                seq2 = seq2[:sg2+1] + '-' + seq2[sg2+1:]
                # Makes no difference if insertion and deletion are in reverse order
    elif len(g1) == 2:
        # If there are two insertions, which is likely an error
        first_start = g1[0][0]
        first_end = g1[0][1]
        second_start = g1[1][0]
        # If there are one or two bases stranded between the gaps 
        if (second_start - first_end) < 4:
            stranded_bases = seq1[first_end:second_start - 1]
            gap_insert = repeat_to_length("-",len(stranded_bases))
            seq1 = seq1[:first_start-1] + stranded_bases + gap_insert + seq1[first_start-1:first_end] + seq1[second_start-1:]
            seq2 = seq2[:first_start-1] + gap_insert + seq2[first_start-1:]

    cseqs = [seq1,seq2]
    return cseqs

def find_seq(seq,fasta):
    '''Locally align two sequences and return start and stop
    coordinates of alignment.
    '''
    coords = ()
    seq = seq.replace('-', '')
    # Read in fasta sequence of target sequence and surrounding region
    record = SeqIO.read(fasta, "fasta")
    refseq = record.seq
    # Penalty for gap extend is set to zero
    # This helps address the issue of chimeric inserts from different regions
    alignments = pairwise2.align.localms(seq, refseq, 2, -5, -4, 0,one_alignment_only=True)
    #print("Searching for: ", seq)
    #print(format_alignment(*alignments[0]))
    # Try block here because alignments may be empty
    try:
        balign = alignments[0][0]
        bscore = alignments[0][2]
        bstart = alignments[0][3]
        bstop = alignments[0][4]
        # Max alignment score with match reward of 2
        max_score = len(seq) * 2
        scorelen = bscore / (len(seq) * 2)
        #print("Region Length*2 / Score = ", scorelen)
        # Ignore alignments with a penalty >= 5
        # Short sequences <6bp are dangerous
        #if bscore > (max_score - 5):
        #if bscore > (max_score / 2):
        coords = (bstart,bstop,scorelen,balign)
    except:
        pass

    return coords


def tandup(alseq,refseq,minin, ref_fasta):
    '''Detect perfect tandem duplications associated with one or more insertions
    in the modified sequence.By default a minimum insert length of 4 is used.
    '''
    gapchar = "-"
    #gl = getgaps(refseq)
    tdcount = 0
    alresults = {}
    len_inseqs = []
    #realigned = realign(refseq,alseq)
    #print("## ARMIN ALIGNMENT ##")
    #print(realigned[0])
    #print(realigned[1])
    #corrected = correct_realign(realigned[0],realigned[1])
    #print("## ARMIN ALIGNMENT CORRECTED ##")
    #print(corrected[0])
    #print(corrected[1], "\n")
    #cor_alseq = corrected[1]
    #cor_refseq = corrected[0]
    gl = getgaps(refseq)
    #print("\n### START NEW ALLELE ###\n")
    #print("CAS-ANALYSER ALIGNMENT")
    #print("REF:",refseq)
    #print("ALN",alseq)
    #print("\nARMIN ALIGNMENT")
    #print("REF:",refseq)
    #print("MUT:",alseq)
    for i in range(len(gl)):
        inseq = alseq[gl[i][0]-1:gl[i][1]]
        #print("\nInsert sequence: ", inseq)
        len_inseqs.append(len(inseq))
        # Check for tandem duplication if insertion is big and no deletion >2bp
        if len(inseq)>=int(minin) and "---" not in alseq:

            # If final or only gap
            #if i + 1 == len(gl):
                #flank = alseq[gl[i][1]:]
            # If not final or only gap only get flank sequence up to next gap
            #else:
            #    flank = alseq[gl[i][1]:gl[i+1][0]-1]

            # Get all sequence before gap as left flank
            rflank = alseq[gl[i][1]:]
            lflank = refseq[:gl[i][0]-1]
            # strip leading gap characters 
            rflank = rflank.strip('-')
            lflank = lflank.strip('-')
            # Get last len(inseq) characters of left flank
            lflank = lflank[-len(inseq)-15:]
            rflank = rflank[:len(inseq)+15]
            #print("Flank  sequence: ", rflank)
            #print("Left flank sequ: ", lflank)
            # Find perfect tandem duplications
            dup = os.path.commonprefix([inseq, rflank])
            dupl = os.path.commonprefix([inseq, lflank])
            # Perfect tandem duplications
            if len(dup) == len(inseq) or len(dupl) == len(inseq):
                tdcount = tdcount + 1
                ##print("\nINSERTION TYPE: " + ins_type)
            # Check for imperfect tandem duplications using alignment
            else:
                ralignments = pairwise2.align.localms(inseq, rflank, 2, -5, -4, -0.1, one_alignment_only=True)
                lalignments = pairwise2.align.localms(inseq, lflank, 2, -5, -4, -0.1, one_alignment_only=True)

                #print("\nAlignment with flanks: \n")
                #print(format_alignment(*ralignments[0]))
                #print(format_alignment(*lalignments[0]))
                # Number of gaps or mismatches in alignment
                ###print("length of insert and score ", len(inseq), alignments[0][2])
                #gm_num = len(inseq) - alignments[0][2]
                # Fraction of bases aligned
                ralscore = ralignments[0][2]
                rscorelen = ralscore / (len(inseq) * 2)
                lalscore = lalignments[0][2]
                lscorelen = lalscore / (len(inseq) * 2)
                #print("Right Flank Align Length*2 / Score = ", rscorelen)
                #print("Left Flank Align Length*2 / Score = ", lscorelen)
                # Identify source of insert sequence through regional alignment
                coords = find_seq(inseq, ref_fasta)
                region_scorelen = coords[2]
                region_alignment = coords[3]
                gaps_region = getgaps(region_alignment)
                #print(gaps_region)
                #print(gaps_to_bases(gaps_region))
                #print("Ins_Alignment_vs_Region", region_alignment)
                #print("Insertion vs ref fasta alignment position", coords)
                # Retrieve position of full right insertion flank
                coords_flank = find_seq(alseq[gl[i][1]:], ref_fasta)
                #print("Flank vs ref fasta alignment position", coords_flank)
                rflank_region =  getgaps(coords_flank[3])

                if rscorelen > 0.9 or lscorelen > 0.9:
                    tdcount = tdcount + 1
                else:
                    # Provide dictionary with region alignmnet results
                    # random two digit number so key names for multiple inserts are not same
                    rands = random_string(2)
                    key1 = "Score_rflank_" + rands
                    key2 = "Score_lflank_" + rands
                    key3 = "Score_region_" + rands
                    key4 = "Pos_rflank_" + rands
                    key5 = "Pos_ins_" + rands
                    alresults[key1] = rscorelen
                    alresults[key2] = lscorelen
                    alresults[key3] = region_scorelen
                    alresults[key4] = gaps_to_bases(rflank_region)
                    alresults[key5] = gaps_to_bases(gaps_region)

                #start_flank = coords_flank[0]
                #start_ins = coords[0]
                #stop_ins = coords[1]
                #dist_r = abs(start_flank - start_ins)
                #dist_l = abs(stop_ins - start_flank)

                #if flank_scorelen < region_scorelen:
                #    ins_type = "Patch insertion"
                #if flank_scorelen >= region_scorelen or dist_r < 15 or dist_l < 15:
                #    ins_type = "Imperfect tandem duplication"
                #print("\nINSERTION TYPE: " + ins_type)
                # Count as TD if less than two mismatches or over 90% matches
                #if gm_num < 2 or al_fr > 0.9:
                #if gm_num < 2:
                #    tdcount = tdcount + 1
                #print("Tandem duplication insert and flank sequences\n",inseq,flank)
               #     coords = find_seq(inseq, ref_fasta)
               #     coords_flank = find_seq(flank, ref_fasta)
                   # print("Tandum duplication, insert and flank positions\n",coords, coords_flank)

                #else:
                    #coords = find_seq(inseq, ref_fasta)
                    #coords_flank = find_seq(flank, ref_fasta)
                    #print("Non tandem dup: Insert and flank sequences\n",inseq, flank)

    if max(len_inseqs) < int(minin): #Count short insert lens as NA
        tdcount = "NA"
    alresults["tdcount"] = tdcount

    return alresults

def check_if_mult_deletions(seq):
    #contiguous series of "-" is a deletion
    multigap = re.compile('-+[ACTGactg]+-')
    multigapsearch = multigap.search(seq)
    return multigapsearch

def getgaps(Seq):
    '''
    Returns a nested list of gap start and end coordinates in a
    nucleotide sequence. Also works for single gaps.
    '''
    gaplist = []
    j=1
    Len=len(Seq)
    start=0
    seq=Seq
    if Seq[0]!='N' and "-" in Seq:
        while len(seq)>0:
            gap_s=seq.find('-')
            start=start+gap_s+1
            seq=Seq[start-1:]
            pool = [seq.find('A'),seq.find('T'),seq.find('G'),seq.find('C')]
            if max(pool)!=-1:
                nul=[ n for n in pool if n>=0]
                end= start+ min(nul)-1
                tmpl = [start,end]
                gaplist.append(tmpl)
                start=end
                seq=Seq[start:]
                j+=1
                if seq.find('-') == -1:
                    break
            else:
                end=Len
                tmpl = [start,end]
                gaplist.append(tmpl)
                j+=1
                break
    return gaplist

def cut_mismatch(aftergaprefseq, aftergapalignseq):
    refbases = []
    alignbases = []
    matchbases = []
    mismatch = "N"
    mm_results = []
    for base in aftergaprefseq: #add each base to list
        refbases.append(base)
    for base in aftergapalignseq:
        alignbases.append(base)
    if len(refbases) == len(alignbases): #ensure sequences are same length
        for i in range(len(refbases)):
            if refbases[i] == alignbases[i]:
                matchbases.append(refbases[i])#append only matching bases, all mismatch positions discarded
            else:
                mismatch = "Y"
    else:
        print("Error: aftergaprefseq and aftergapalignseq have different lengths!")
    matchbaseseq = ''.join(matchbases)
    mm_results.extend((matchbaseseq,mismatch))

    return mm_results


#def find_single_mh(alignseq,refseq, min_mh_len):
def find_single_mh(alignseq,refseq):
    result = []
    gap = re.compile('-+')
    gapseq = gap.search(alignseq).group()
    gapstart = gap.search(alignseq).start()
    gapend = gap.search(alignseq).end()

    gaprefseq = refseq[gapstart:gapend]
    aftergaprefseq = refseq[gapend:]
    aftergapalignseq = alignseq[gapend:]
    #with mismatch (wmm) sites used
    mhseq_wmm = os.path.commonprefix([gaprefseq, aftergaprefseq])

    mm_out = cut_mismatch(aftergaprefseq, aftergapalignseq)
    matchseq = mm_out[0]
    mm = mm_out[1]
    mhseq_womm = os.path.commonprefix([gaprefseq, matchseq])
    if mhseq_wmm is not '' or mhseq_womm is not '':
        #if len(mhseq_wmm) >= min_mh_len or  len(mhseq_womm) >= min_mh_len:
        result.append(mhseq_wmm)
        result.append(mhseq_womm)
    #whether a mismatch occurs should be reported for all alleles with a deletion
    result.append(mm)#add whether a mismatch occurs in aftergaprefseq relative to aftergapalignseq
    return result

def find_mult_mh(alignseq,refseq):
    gaplist = getgaps(alignseq)
    mh_list =[]
    for i in range(len(gaplist)):
        if i + 1 == len(gaplist):
            gaprefseq = refseq[gaplist[i][0]-1:gaplist[i][1]]
            aftergaprefseq = refseq[gaplist[i][1]:]
            aftergapalignseq = alignseq[gaplist[i][1]:]
        else:
            gaprefseq = refseq[gaplist[i][0]-1:gaplist[i][1]]
            aftergaprefseq = refseq[gaplist[i][1]:gaplist[i+1][0]-1]
            aftergapalignseq = alignseq[gaplist[i][1]:gaplist[i+1][0]-1]

        mhseq_wmm = os.path.commonprefix([gaprefseq, aftergaprefseq])
        mm_out = cut_mismatch(aftergaprefseq, aftergapalignseq)
        matchseq = mm_out[0]
        mm = mm_out[1]
        mhseq_womm = os.path.commonprefix([gaprefseq, matchseq])

        if mhseq_wmm is not '' or mhseq_womm is not '':
            mhpair = [mhseq_wmm,mhseq_womm,mm]
            mh_list.append(mhpair)
        else:
            mh_list.append(mm)

    return mh_list


def main(align,path,name,tool,min_mh_len,minin,ref_fasta,debug,stats):
    count = 0
    single_count = 0
    tot_single_count = 0
    mult_count = 0
    tot_mult_count = 0
    all_mh_wgaps = {}
    all_mh_wogaps = {}
    with open ('%s/%s_mh.txt'%(path, name), 'w') as mh_table:
    #wgaps = mismatches included
    #wogaps = mismatches removed 
        with open(align, 'r') as f:
            if debug: print("Opening " + align)
            reader = csv.reader(f, dialect='excel', delimiter='\t')
            #next(reader, None)  # skip the headers
            headers = next(reader, None) #stores header for later use
            headers = '\t'.join(headers)
            mh_table.write('%s\tDel_len\tMH_check\tMH\tMH_wmm\tMH_mm_len\tMH_nomm\tMH_nomm_len\tMismatch\tReadID\tTanDup\tInsDict\n'%(headers))
            for row in reader:
                if tool == "CE":
                    count = count + 1
                    alseq = row[0]
                    refseq = row[1]
                    unmodified = row[3]
                    #n_inserted = float(row[6]) 
                    n_deleted = float(row[5])
                    read_count = row[8]
                    if debug: print("ID is " + str(count) + ", TargetSeq is " + "\n" + alseq + "\n" + ", RefSeq is " + "\n" + refseq + "\n" + ", Edit type is: " + unmodified + ", Read count is: " + str(read_count) + "\n")
                elif tool == "CA":
                    count = count + 1
                    alseq = row[2]
                    refseq = row[1]
                    realigned = realign(refseq,alseq)
                    corrected = correct_realign(realigned[0],realigned[1])
                    alseq = corrected[1]
                    refseq = corrected[0]
                    edit_type = row[5] # "Ins|del|WT or Sub"
                    read_count = row[4]
                    if debug: print("ID is " + str(count) + ", TargetSeq is " + alseq + ", RefSeq is " + refseq + ", Edit type is: " + edit_type + ", Read count is: " + str(read_count) + "\n")
                else:
                    print("Error! Tool must be CA or CE")

                if len(alseq) != len(refseq):
                    print("Warning: Aligned sequence and reference have different lengths! Microhomology calculations may be incorrect")

                printrow = "\t".join(row)

                deletion = False
                insert = False

                if tool == "CE":
                    if unmodified == "False" and n_deleted > 0:
                        deletion = True
                elif tool == "CA":
                    if edit_type == "del":
                        deletion = True
                    elif edit_type == "Ins":
                        insert = True

                if not deletion and not insert:
                        if debug: print("No deletion identified in allele " + str(count) + "\n")
                        mh_table.write('%s\tNA\tN\tN\t.\t.\t.\t.\tNA\t%s\tNA\tNA\n'%(printrow, count))
                else:
                    if insert:
                            td = tandup(alseq,refseq,minin,ref_fasta) #get number of tandem duplications
                            td_count = td["tdcount"]
                            mh_table.write('%s\tNA\tN\tN\t.\t.\t.\t.\tNA\t%s\t%s\t%s\n'%(printrow,count,td_count,td))
                    else:
                        td = "NA" #rare cases when there are inserts and dels do not have CasAn category
                        td_count = "NA"
                        del_ss = getgaps(alseq) #gets nested list of start and stop positions of gaps in the target sequence
                        if check_if_mult_deletions(alseq): #if not empty
                            #Starting Multi Gap MH analysis
                            if debug: print("Identified multiple deletion in allele " + str(count) + "\n")
                            tot_mult_count = tot_mult_count + 1
                            res = find_mult_mh(alseq,refseq)
                            mh_wgaps = []
                            mh_wgaps_l = []
                            mh_wogaps = []
                            mh_wogaps_l = []
                            mult_count = mult_count + 1
                            id_count = 0 #using id_count we can get the len of each del from del_len even for multiple deletions
                            for mh in res:
                                my_del_ss = del_ss[id_count] #gets start and end pos for the relevant deletion in the nested list
                                del_len = my_del_ss[1] - my_del_ss[0] + 1  #gets deletion length by subtracting end pos by start pos of gap plus 1 
                                id_count = id_count + 1
                                res_s = ','.join(str(v) for v in res) #stringify for print on L211
                                printcount = str(count) + "_" + str(id_count)
                                myfirst = str(mh[0])
                                if mh[0] is "Y" or mh[0] is "N": #if no microhomology only the mismatch pseudobool will be in list item
                                    mh_table.write('%s\t%s\tY\tN\t.\t.\t.\t.\t%s\t%s\t%s\t%s\n'%(printrow,del_len,mh[0],printcount,td_count,td))
                                else:
                                    mh_wgaps.append(mh[0])

                                    if debug: print("Microhomology found in allele " + str(count) + ". Microhomology is: " + res_s + "\n")
                                    mykey = str(count) + "_" + str(read_count)
                                    if mh[0] != '':
                                        all_mh_wgaps[mykey] = mh[0]
                                        if debug: print("Appending the following non-empty microhomology (with mismatches): " + str(mh[0]) + "\n")
                                    mh_wgaps_l.append(str(len(mh[0])).replace('0','.'))
                                    mh_wogaps.append(mh[1])
                                    if mh[1] != '':
                                        all_mh_wogaps[mykey] = mh[1]
                                        if debug: print("Appending the following non-empty microhomology (without mismatches): " + str(mh[1]) + "\n")
                                    mh_wogaps_l.append(str(len(mh[1])).replace('0','.'))

                                    printlen_wmm = str(len(mh[0]))
                                    printlen_womm = str(len(mh[1]))
                                    if mh[0]:
                                       printmh_wmm = mh[0]
                                    else:
                                        printmh_wmm = "."
                                    if mh[1]:
                                       printmh_womm = mh[1]
                                    else:
                                        printmh_womm = "."

                                    mh_table.write('%s\t%s\tY\tY\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(printrow,del_len,printmh_wmm,printlen_wmm,printmh_womm,printlen_womm,mh[2],printcount,td_count,td)) #mh[2] is the mismatch bool 

                        else:
                            #Starting single gap MH analysis
                            if debug: print("Identified single deletion in allele " + str(count) + "\n")
                            del_len = del_ss[0][1] - del_ss[0][0] + 1 #get deletion length by subtracting end pos by start pos of gap plus 1 
                            tot_single_count = tot_single_count + 1
                            res = find_single_mh(alseq,refseq)
                            if res[0] is "Y" or res[0] is "N": #if no microhomology only the mismatch pseudobool will be in list item
                                mh_table.write('%s\t%s\tY\tN\t.\t.\t.\t.\t%s\t%s\t%s\t%s\n'%(printrow,del_len, res[0], count, td_count,td))
                                if debug: print("No microhomology in allele " + str(count) + ". Printing empty columns to row.\n")
                            else:
                                res_s = ','.join(str(v) for v in res)
                                if debug: print("Microhomology found in allele " + str(count) + ". Microhomology is: " + res_s + "\n")
                                mykey = str(count) + "_" + str(read_count)
                                single_count = single_count + 1
                                if res[0] != '':
                                    all_mh_wgaps[mykey] = res[0]
                                    if debug: print("Appending the following non-empty microhomology (with mismatches): " + str(res[0]) + "\n")
                                if res[1] != '':
                                    all_mh_wogaps[mykey] = res[1]
                                    if debug: print("Appending the following non-empty microhomology (without mismatches): " + str(res[1]) + "\n")
                                printlen_wmm = str(len(res[0]))
                                printlen_womm = str(len(res[1]))
                                if res[0]:
                                   printmh_wmm = res[0]
                                else:
                                    printmh_wmm = "."
                                if res[1]:
                                   printmh_womm = res[1]
                                else:
                                    printmh_womm = "."
                                mh_table.write('%s\t%s\tY\tY\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(printrow,del_len,printmh_wmm,printlen_wmm,printmh_womm,printlen_womm,res[2],count,td_count,td)) #res[2] is the mismatch bool 

    print("Done! " + str(count) + " alleles written to " + path + "/" + name + "_mh.txt")
    #Summary statistics
    if stats == 1:
        with open ('%s/%s_mh_stats.txt'%(path, name), 'w') as mh_stats:
            totalmh = mult_count + single_count
            mh_stats.write("Total number of single deletion alleles: " + str(tot_single_count) + "\n")
            mh_stats.write("Microhomologies identified in single deletion alleles: " + str(single_count) + "\n")
            mh_stats.write("Total number of multiple deletion alleles: " + str(tot_mult_count) + "\n")
            mh_stats.write("Microhomologies identified in multiple deletion alleles: " + str(mult_count) + "\n")
            mh_stats.write("Total alleles with microhomologies found (with and without mismatches): " + str(totalmh) + "\n")

            mylist1 = []
            mylist2 = []
            lenset1 = set()
            lenset2 = set()
            countdic1 = {}
            countdic2 = {}
            print("Starting analysis of homologies with mismatches")
            for key, value in all_mh_wgaps.items():
                nbases = len(value)
                nreads = int(key.split("_")[1])
                if debug: print("Key is " + key + ", Value is: " + value + ", NReads is : " + str(nreads) + "\n")
                mylist1.append(nbases)
                if nbases in lenset1:
                    current_nreads = int(countdic1[nbases])
                    incremented_nreads = current_nreads + nreads
                    if debug:
                        print("Current reads for " + str(nbases) + " are : " + str(current_nreads))
                        print("Incremented reads for " + str(nbases) + " are : " + str(incremented_nreads))
                    countdic1[nbases] = incremented_nreads
                else:
                    lenset1.add(nbases)
                    countdic1[nbases] = nreads

            print("Completed analysis of homologies with mismatches. Starting analysis of homologies without mismatches.")
            for key, value in all_mh_wogaps.items():
                nbases = len(value)
                nreads = int(key.split("_")[1])
                if debug: print("Key is " + key + ", Value is: " + value + ", NReads is : " + str(nreads) + "\n")
                mylist2.append(nbases)
                if nbases in lenset2:
                    current_nreads = int(countdic2[nbases])
                    incremented_nreads = current_nreads + nreads
                    if debug:
                        print("Current reads for " + str(nbases) + " are : " + str(current_nreads))
                        print("Incremented reads for " + str(nbases) + " are : " + str(incremented_nreads))
                    countdic2[nbases] = incremented_nreads
                else:
                    lenset2.add(nbases)
                    countdic2[nbases] = nreads

            counts1 = Counter(mylist1)
            counts2 = Counter(mylist2)
            mh_stats.write("Allele counts of microhomology length with mismatches included" + "\n")
            for key, count in sorted(counts1.most_common()):
                mh_stats.write('%s: %s\n'%(key, count))
            mh_stats.write("Read counts of microhomology length with mismatches included" + "\n")
            for key, count in countdic1.items():
                mh_stats.write('%s: %s\n'%(key, count))
            mh_stats.write("Allele counts of microhomology length with mismatches excluded" + "\n")
            for key, count in sorted(counts2.most_common()):
                mh_stats.write('%s: %s\n'%(key, count))
            mh_stats.write("Read counts of microhomology length with mismatches excluded" + "\n")
            for key, count in countdic2.items():
                mh_stats.write('%s: %s\n'%(key, count))

            print("All done!")

parser = ArgumentParser(description='Search Cas-Analyser or CRISPResso allele outputs for microhomologies associated with deletions')
parser.add_argument('-v', '--version', action='version', version='1.1')
parser.add_argument('-i', dest='table', help='the tab-separated allele output file')
parser.add_argument('-p', dest='prefix', help='prefix of outputs. Default: prefix of the input file')
parser.add_argument('-r', dest='reference', help='FASTA format sequence to search for source of inserted sequences')
parser.add_argument('-o', dest='output', help='Output directory')
parser.add_argument('-t', dest='tool', help='Tool used to generate input alleles; can be Cas-Analyser or CRISPResso: CA|CE')
parser.add_argument('-m', dest='min_mh_len', help='Minimum length of microhomology in basepairs [int]')
parser.add_argument('-n', dest='minin', help='Minimum length of insert in basepairs [int]. Default is 4')
parser.add_argument('-s', dest='stats', help='Generate microhomology statistics file')
parser.add_argument('--verbose', action='store_true', help='Verbosely prints all results at each step')
args = parser.parse_args()

if None not in [args.table, args.output]:
    args.output = os.path.abspath(args.output)
    if not os.path.isdir(args.output):
        print('\nOops! It seems the path to output directory does not existent. Please check!\n')
        sys.exit(1)
    if args.tool == None:
        print('\nPlease set the tool used to generate the allele table!\n')
        sys.exit(1)
    if args.reference == None:
        print('\nPlease provide FASTA reference file!\n')
        sys.exit(1)
    if args.verbose:
        debug = 1
    else:
        debug = None
    if args.stats:
        stats = 1
    else:
        stats = None
    if args.prefix == None:
        args.prefix = args.table.rsplit('.',1)[0].split('/')[-1]
        #main(args.table, args.output, args.prefix, args.tool, args.min_mh_len, debug,stats)
    if args.minin == None:
        args.minin = 4
    main(args.table, args.output, args.prefix, args.tool, args.min_mh_len, args.minin,args.reference,debug,stats)

else:
    print
    parser.print_help()
    print
    sys.exit(1)


