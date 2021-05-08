## Detecting InDels in PacBio long amplicon reads 

The pipeline to assess InDel rates can be executed on a set of example reads as shown below.

``` bash pb_pipeline.sh ```

This generates a set of output files with all InDels and their distance from the nick regions as defined in ```data/cut_positions.txt```. The main final output (```all_loci_maxlen_q10_nickregion_master.tsv```) is summarised in a table of statistics (```all_loci_maxlen_q10_nickregion_master_statistics.tsv```).

The file ```results/all_indels.tsv``` contains the results of InDel analysis for all mutational backgrounds.
