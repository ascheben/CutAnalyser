## **_CutAnalyser_**

A python tool to analyse genome editing cut sites using amplicon sequencing data

## Introduction

CutAnalyser is a single-purpose utility script written in python 3 and tested in an Ubuntu environment. Although extensive testing has not been carried out, the program reliably produces the expected output under a range of conditions and with different inputs. Nevertheless, as the program was not written for broad application, users should be cautious.

## Usage

CutAnalyser.py [options] [-o OUTPUT] [-t TOOL] [-m MIN\_MH\_LEN] [-i TABLE]

|Parameter| Description|
| -h | Show help message and exits |
| -v | Show program version number and exits |
| -i | The tab-separated allele file that provides the main input data for the program. The table must contain a column with the wild type allele sequence, the modified allele sequence and a column indicating whether the mutation type is a deletion. |
| -p | Prefix of output file (default: prefix of the input file) |
| -o | Output directory to print result files to |
| -t | Tool used to generate input file provided via -i; can be Cas-Analyser orCRISPResso (CA|CE). |
| -m | Minimum length of microhomology in basepairs [int]. All microhomologies below this length will be excluded. |
| verbose | Verbosely prints all results at each step. This can be useful for ensuring the program is handling the input data in the expected way. |

Example:

`CutAnalyser.py -i CRISPResso\_Alleles\_frequency\_table.txt  -o /home/ws/crispr/mhFinder/ -t CE -m 4`

## Required input

The program requires the allele frequency output of CRISPResso ([https://github.com/lucapinello/CRISPResso](https://github.com/lucapinello/CRISPResso)) or CasAnalyser ([https://github.com/snugel/cas-analyzer](https://github.com/snugel/cas-analyzer)) to carry out the microhomology search in deletion alleles. The information that mhFinder extracts from the allele frequency table is 1) Wild type sequence, 2) Treated sequence, 3) mutation type, 4) Read count for allele. The numbered columns of the CRISPResso Allele\_Frequency.txt table are shown below, with columns used by mhFinder shown in bold. The Unmodified column can be "True" or "False", and the n\_deleted column is an integer showing the number of deleted bases.

| 1 | 2 | 3 | 4| 5 | 6 | 7 | 8 | 9 | 10 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **Aligned\_ Sequence** | **Reference\_Sequence** | NHEJ | **Unmodified** | HDR | **n\_deleted** | n\_inserted | n\_mutated | **#Reads** | %Reads |



The Cas-Analyzer table is extracted from the web app ([http://www.rgenome.net/cas-analyzer](http://www.rgenome.net/cas-analyzer)), which allows user to export the allele frequency tables. The columns are shown below, with columns used by mhFinder again shown in bold. The Type column can be "WT or Sub", "del" or "Ins".

| 1 | 2 | 3 | 4 | 5 | 6 | 7 |
| --- | --- | --- | --- | --- | --- | --- |
| #ID | **WT Sequence** | **RGEN Treated Sequence** | Length | **Count** | **Type** | HDR |

Allele frequencies generated using other software can be used as input for mhFinder, if they are reformatted to resemble Cas-Analyzer or CRISPResso input.

## Outputs

A successful run of mhFinder will produce two output files 1) [prefix]\_mh.txt and 2) [prefix]\_mh\_stats.txt. The first file is a tab-separated file containing all of the columns of the input file with an additional 8 columns that contain microhomology information appended to it (Table 1).

**Table 1**. Columns appended to input table by mhFinder.

| **Column name** | **Description** | **Data type** |
| --- | --- | --- |
| MH\_check | Identifies whether a microhomology search has been carried out on the allele. | Y|N |
| MH | Identifies whether a microhomology was detected. | Y|N |
| MH\_wmm | Shows the microhomology found when mismatch sites are included as part of the homology | nucleotide string |
| MH\_wmm\_len | Length of microhomology found (mismatch sites included) | [int] or . |
| MH\_nonn | Shows the microhomology found when mismatch sites are excluded before the homology search | nucleotide string |
| MH\_nomm\_len | Length of microhomology found (mismatch sites excluded) | [int] or . |
| Mismatch | Identified whether a mismatch was identified. Defaults to not applicable (NA) when there is no deletion. | Y|N|NA |
| ReadID | Running integer IDs given to each allele. Allows alleles with multiple deletions to be split into multiple alleles (e.g., 10\_1, 10\_2 and 10\_3) for easier parsing. | [int] or [int] + "_" + [int] |

The results include "microhomologies with mismatches" and "microhomologies without mismatches". The difference can be explained using the example below showing a reference allele aligned with the corresponding treated allele.

Reference allele:        AGG**CGC**CTT

Treated allele:          AGG**---**GTT

The treated allele in the example above contains a deletion at positions 4-6. The reference sequence at the same positions is &quot;CGC&quot;, and the reference sequence following this sequence is &quot;CTT&quot;. We can therefore see that position 7 in the alignment contains a mismatch between the reference and the treated allele (C/G). We can now carry out a microhomology search 1) with mismatches and 2) without mismatches.

1. Compare **C**GC with **C**TT for homology site by site. Result: a single nucleotide homology of C.
2. Remove mismatch sites between CTT and GTT, leaving TT. Now compare CGC with TT for homology. Result: no microhomology.

The second output produced by mhFinder is an experimental summary statistic for the microhomologies found in the sample (Table 2). Particularly the treatment of alleles with multiple deletions in this summary statistics is experimental. It is advisable to also calculate your own summary statistics with the output table to ensure the expected results are produced. The script getstats.sh (described below) is currently the preferred method for calculating summary statistics, handling alleles with multiple microhomologies better.

**Table 2**. Statistics provided in the mhFinder output. Microhomology statistics depend on the minimum microhomology length set by the user.

| **Statistic name** | **Description** |
| --- | --- |
| Total number of single deletion alleles | Sum of treated alleles with a single deletion (as opposed to no deletion or multiple deletions) |
| Microhomologies identified in single deletion alleles | Number of alleles with single deletion that showed a microhomology. |
| Total number of multiple deletion alleles | Sum of treated alleles with more than one deletion (as opposed to no deletion or one deletion). Each allele is counted once, regardless of the number of deletions it has. |
| Microhomologies identified in multiple deletion alleles | Sum of alleles with multiple deletions and at least one microhomology. An allele with, e.g., two deletions and two microhomologies is only counted once. |
| Total alleles with microhomologies found (with and without mismatches) | Sum of alleles with at least one microhomology (alleles with multiple microhomologies only counted once) |
| Allele counts of microhomology length | Sum of allele counts for all alleles with microhomologies of a given length (only the microhomology closes to the end of the sequence is counted in alleles with multiple microhomologies, others are ignored) |
| Read counts of microhomology length | Sum of allele read counts for all alleles with microhomologies of a given length (only the microhomology closes to the end of the sequence is counted in alleles with multiple microhomologies, others are ignored) |

## Example output analysis

A typical analysis may involve multiple treatments, each of which has its own allele frequency table that is analysed with mhFinder. By adding an additional column with the treatment name, we can merge the tables together to get a single output table.

`sed 's/^/mytreatment\t/' mysample_result_mh.txt > mysample_result_mh.id.txt`

`for l in *mh.id.txt; do tail -n +2 $l >> mh_all.txt;done`

The table will look like this (column numbers, column names and first row of example values shown).

| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8|
| --- | --- | --- | --- | --- | --- | --- | --- |
| group | #ID | WT Sequence | RGEN Treated Sequence | Length | Count | Type | HDR |
| mytreatment1 | 1 | ATGTG | A--TG | 3 | 100 | del | false |

| 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| MH_check | MH | MH_wmm | MH_wmm_len | MH_womm | MH_womm_len | MM | ReadID |
| Y | Y | TG | 2 | TG | 2 | N | 1 |

To get separate statistics for each of the treatment groups, we then extract the unique treatment groups and print them to a file, so that each line has one group.

`tail -n +2 mh_all.txt | cut -f1 | sort | uniq > treatments.txt`

`cat treatments.txt`

`treatment1`
`treatment2`
`treatment3`

We can then obtain summary statistic using the getstats.sh script.

usage: getstats.sh [min_homology_len] [treatments_file] [CutAnalyser_out_table]

Example: `bash getstats.sh 2 treatments.txt mh_all.txt`

This would print statistics for all alleles with microhomologies of at least two nucleotides length. By default the script excludes alleles with mismatches. This behaviour can only be turned off by manually editing the script.

## Microhomology analysis on example file
The results table for the example data provided on ([http://www.rgenome.net/cas-analyzer](http://www.rgenome.net/cas-analyzer)) was downloaded and added to this repository for testing. The file can be analysed with mhFinder:

`CutAnalyser.py -i result.txt -o /home/arminps/ws/test/ -t CA`

Because there is only one treatment, we need to coerce the results_mh.txt file to allow us to run `getstats.sh`. 

`sed -i 's/^/treatment\t/' results_mh.txt`
`echo "treatment" >> treatments.txt`

Now that we have added a treatment column and a mock treatments file. We can run `getstats.sh` on our single treatment.

`bash getstats.sh 1 alltreat.txt result_mh.txt`
