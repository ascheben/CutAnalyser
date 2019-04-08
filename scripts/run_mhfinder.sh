CUT="CutAnalyser.py"
GETSTATS="scripts/getstats.sh"
MINMH=3
MININ=4
MINDEL=3

#Find MH in each sample
for result in *.txt; do python3 $CUT -i $result -o . -t CA -m $MININ;done
#Add ID column based on file name
for mh in *_mh.txt; do sed "s/^/${mh%%_mh.txt}\t/" $mh > ${mh%%_mh.txt}_mh.id.txt; rm $mh;done
#Concatenate results
for mhid in *_mh.id.txt; do tail -n +2 $mhid >> mh_all.txt; rm $mhid;done
#Extract treatments
cut -f1 mh_all.txt | sort | uniq > treatments.txt
#Output statistics table
bash $GETSTATS $MINDEL $MINMH treatments.txt mh_all.txt
rm treatments.txt mh_all.txt
