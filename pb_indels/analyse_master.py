#!/usr/bin/env python3
import csv
import sys
from collections import defaultdict

in_table = sys.argv[1]
out_name = in_table[0:-3] + "filt.txt"

# Set maximum distance from expected indel region based on nick sites
# Note MAXDIST=1 means that there is no distance and regions overlap
MAXDIST=10000
# Set minimum length of indels to exclude PacBio errors
MINLEN=3

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
            pool = [seq.find('A'),seq.find('T'),seq.find('G'),seq.find('C'),seq.find('a'),seq.find('t'),seq.find('g'),seq.find('c')]
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

def solve(r1, r2):
     # sort the two ranges such that the range with smaller first element
     # is assigned to x and the bigger one is assigned to y
     x, y = sorted((r1, r2))
     #now if x[1] lies between x[0] and y[0](x[1] != y[0] but can be equal to x[0])
     #then the ranges are not overlapping and return the differnce of y[0] and x[1]
     #otherwise return 0 
     if x[0] <= x[1] < y[0] and all( y[0] <= y[1] for y in (r1,r2)):
        return y[0] - x[1]
     return 0

def solve_prefix(refrange, indelrange):
     # sort the two ranges such that the range with smaller first element
     # is assigned to x and the bigger one is assigned to y
     x, y = sorted((refrange, indelrange))
     #now if x[1] lies between x[0] and y[0](x[1] != y[0] but can be equal to x[0])
     #then the ranges are not overlapping and return the differnce of y[0] and x[1]
     #otherwise return 0 
     if x[0] <= x[1] < y[0] and all( y[0] <= y[1] for y in (refrange, indelrange)):
        mydist = y[0] - x[1]
        # if indel downstream of nickrange
        if refrange[1] < indelrange[0]:
            mydist = "+" + str(mydist)
        elif mydist > 0:
            mydist = "-" + str(mydist)
        return mydist
     else:
        return 0

def getmedian(mylist):
    mylist.sort()
    n = len(mylist)
    if n % 2 == 0:
        median1 = mylist[n//2]
        median2 = mylist[n//2 - 1]
        median = (median1 + median2)/2
    else:
        median = mylist[n//2]
    return median

# dicts of sets for each locus
total_reads = defaultdict(set)
wt_reads = defaultdict(set)
del_reads = defaultdict(set)
ins_reads = defaultdict(set)
del_alleles = defaultdict(set)
ins_alleles = defaultdict(set)

with open (out_name, 'w') as out_table:
    with open(in_table, 'r') as f:
        reader = csv.reader(f, dialect='excel', delimiter='\t')
        next(reader)
        dedupdict = {}
        for row in reader:
            locus = row[0]
            name = row[4]
            indeltype = row[5]
            total_reads[locus].add(name)
            if indeltype == "WT":
                wt_reads[locus].add(name)
            else:
                locus = row[0]
                name = row[4]
                numid = row[3]
                indelrange = [int(row[6]),int(row[7])]
                nickrange = [int(row[9]),int(row[10])]
                indel_len = int(row[8])
                mindist = solve(nickrange,indelrange)
                # SWITCH CONTROLLING ALLOWED INDEL REGION
                if mindist < MAXDIST and indel_len >= MINLEN:
                    if locus in dedupdict and name in dedupdict[locus]:
                        old_len = dedupdict[locus][name][0]
                        old_dist = dedupdict[locus][name][1]
                        #if indel_len > old_len:
                        if indel_len > old_len and mindist <= old_dist:
                            dedupdict[locus][name] = [indel_len,mindist,int(row[6])]
                    elif locus in dedupdict:
                        dedupdict[locus][name] = [indel_len,mindist,int(row[6])]
                    else:
                        dedupdict[locus] = {name:[indel_len,mindist,int(row[6])]}
#locus   ref alt number  read_name   type    indel_start indel_end   indel_length    nick_region_start   nick_region_end

    with open(in_table, 'r') as f:
        reader = csv.reader(f, dialect='excel', delimiter='\t')
        indeldict = {}
        for row in reader:
            indeltype = row[5]
            if row[8] == "NA":
                indel_len = 0
            else:
                indel_len = int(row[8])

            if indeltype != "WT" and indel_len >= MINLEN:
                locus = row[0]
                refseq = row[1]
                mutseq = row[2]
                name = row[4]
                indeltype = row[5]
                refstart = int(row[6])
                refstop = int(row[7])
                indel_len = int(row[8])
                indelrange = [int(row[6]),int(row[7])]
                nickrange = [int(row[9]),int(row[10])]
                mindist = solve_prefix(nickrange,indelrange)

                if locus in dedupdict:
                    if name in dedupdict[locus] and dedupdict[locus][name][2] == refstart:
                        if indeltype == 'ins':
                            ins_reads[locus].add(name)
                            ins_seqs = []
                            finalnums = []
                            mydict = {}
                            refgaps = getgaps(refseq)
                            for gapnums in refgaps:
                                start = gapnums[0]
                                stop = gapnums[1]
                                seq = mutseq[start-1:stop]
                                if len(seq) == int(indel_len):
                                    ins_seqs.append(seq)
                                    finalnums.append(gapnums)
                                    mydict[finalnums[0][0]] = seq
                            best_key = min(mydict, key=lambda x:abs(x-int(refstart)))
                            best_seq = mydict[best_key]

                            unique_ins = str(refstart) + "__" + str(refstop) + "__" +  best_seq
                            ins_alleles[locus].add(unique_ins)

                            out_table.write(locus + "\t" + "ins\t" + name + "\t" + str(refstart) + "\t" + str(refstop)+  "\t" + str(best_seq) + "\t" + str(mindist) + "\t" +str(row[9]) + "\t" + str(row[10]) + "\n")
                        elif indeltype == 'del':
                            del_reads[locus].add(name)
                            finalnums = []
                            mydict = {}
                            mutgaps = getgaps(mutseq)
                            for gapnums in mutgaps:
                                start = gapnums[0]
                                stop = gapnums[1]
                                seq = refseq[start-1:stop]
                                if len(seq) == int(indel_len):
                                    finalnums.append(gapnums)
                                    mydict[finalnums[0][0]] = seq
                            best_key = min(mydict, key=lambda x:abs(x-int(refstart)))
                            best_seq = mydict[best_key]
                            unique_del = str(refstart) + "__" + str(refstop) + "__" +  best_seq
                            del_alleles[locus].add(unique_del)
                            out_table.write(locus + "\t" + "del\t" + name + "\t" + str(refstart) + "\t" + str(refstop)+  "\t" + str(best_seq) + "\t" + str(mindist) + "\t" +str(row[9]) + "\t" + str(row[10]) + "\n")


# Header for stats table
outline = ["locus","reads","wild_type_reads","reads_with_deletions","reads_with_insertions","alleles_with_deletions","alleles_with_insertions","del_read_fraction","ins_read_fraction","indel_read_fraction","mean_ins_len","median_ins_len","max_ins_len","min_ins_len","mean_del_len","median_del_len","max_del_len","min_del_len"]
print(*outline,sep="\t")


# Collect stats from dictionaries of indel information
results_dict = {}
for key,value in total_reads.items():
    reads = len(value)
    wild_type_reads = len(wt_reads[key])
    reads_with_deletions = len(del_reads[key])
    reads_with_insertions = len(ins_reads[key])
    alleles_with_deletions = len(del_alleles[key])
    alleles_with_insertions = len(ins_alleles[key])
    del_read_fraction = reads_with_deletions / reads
    ins_read_fraction = reads_with_insertions / reads
    indel_read_fraction = del_read_fraction  +  ins_read_fraction
    # sets of indels lengths per locus

    mean_ins_len = "NA"
    median_ins_len = "NA"
    max_ins_len = "NA"
    min_ins_len = "NA"
    mean_del_len = "NA"
    median_del_len = "NA"
    max_del_len = "NA"
    min_del_len = "NA"

    if len(ins_alleles[key]) > 0:
        insertion_alleles = list(ins_alleles[key])
        inslens = []
        for ins in insertion_alleles:
            sequencelen = len(ins.split("__")[2])
            inslens.append(sequencelen)
        mean_ins_len = sum(inslens) / len(inslens)
        median_ins_len = getmedian(inslens)
        max_ins_len = max(inslens)
        min_ins_len = min(inslens)

    if len(del_alleles[key]) > 0:
        deletion_alleles = list(del_alleles[key])
        dellens = []
        for deletion in deletion_alleles:
            sequencelen = len(deletion.split("__")[2])
            dellens.append(sequencelen)
        mean_del_len = sum(dellens) / len(dellens)
        median_del_len = getmedian(dellens)
        max_del_len = max(dellens)
        min_del_len = min(dellens)


# Dict may be useful for further processing
#    results_dict[key] = {
#    "reads" : reads,
#    "wild_type_reads" : wild_type_reads,
#    "reads_with_deletions" : reads_with_deletions,
#    "reads_with_insertions" : reads_with_insertions,
#    "alleles_with_deletions" : alleles_with_deletions,
#    "alleles_with_insertions" : alleles_with_insertions,
#    "del_read_fraction" : del_read_fraction,
#    "ins_read_fraction" : ins_read_fraction,
#    "indel_read_fraction" : indel_read_fraction,
#    "mean_ins_len" : mean_ins_len,
#    "median_ins_len" : median_ins_len,
#    "max_ins_len" : max_ins_len,
#    "min_ins_len" : min_ins_len,
#    "mean_del_len" : mean_del_len,
#    "median_del_len" : median_del_len,
#    "max_del_len" : max_del_len,
#    "min_del_len" : min_del_len
#    }

    # Write output line
    outline = [key, reads, wild_type_reads, reads_with_deletions, reads_with_insertions, alleles_with_deletions, alleles_with_insertions, del_read_fraction, ins_read_fraction, indel_read_fraction, mean_ins_len, median_ins_len, max_ins_len, min_ins_len, mean_del_len, median_del_len, max_del_len, min_del_len]
    print(*outline,sep="\t")
