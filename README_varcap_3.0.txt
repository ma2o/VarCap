VARCAP README 3.0
[VARiant Calling And Postfiltering]

1. Preperation of project

Go to the varcap directory and call the 00_setup_project.sh with the absolute path to the raw read file directory
and the desired absolute path to the output directory (If not existent, it will be created.) Alternatively you can also
add the reference as a third option, however the reference can also be added later.

Command: 
bash 00_setup_project.sh /path/to/fastqfolder /path/to/my/outputdir_01 [/path/to/reference.fasta]

The script will create an output directory named e.g. outputdir_01 from the path above, and will create sample subfolders for 
each bam or fastq(.gz) (also paired) file within that project directory. The input files can be within the given directory
or within a subdirectory of the given source directory.The sample subfolders house the main scripts and the variant.config file,
which holds the necessary parameters.
Note that the fastq files must have the suffix *_1.fastq (_1.fastq.gz) for the forward read and *_2.fastq (_2.fastq.gz) for the reverse read.
 
Most values will be estimated/filled automatically, however you have to supply a reference (either at 00_setup_project.sh  script, at A001_A...sh or you can edit the variable assignment REF_FA_ORIGINAL manually.
Go to the newly generated project folder, there you will see the one subfolder per sample and some accompanying scripts, that
will lead you through the variant calling workflow.


2. Quality control and reference processing

Script: 001_A_quality_filter.sh
- Generates quality output for the raw reads (in sample/raw directory)
- Quality filters the reads ( in sample/filter directory)
- Sets and processes reference ( in reference directory)
  If path to reference was not supplied yet, you can add it here

Command: 
bash AAA_quality_check_raw.sh
or with reference:
bash AAA_quality_check_raw.sh [/path/to/reference.fasta]

This script will submit jobs to the cluster, where they will be processed. With qstat -u <username>, you can monitor your jobs.


3. Setup index for mapping and preocede with mapping

Generates a (reusable) index for mapping in the index_bwa subdirectory of the project directory. Each mapping script will check if the index
exists, else it will be created. 
The index will then be used for mapping. 
By default all reads will be mapped, if you just want to map a subset, you can supply the number here.

Command:
bash 002_A_mapping.sh
or with mapping count, to map 100000 reads:
bash 002_A_mapping.sh 100000

Here you can use 2 modes. Either you just map to the defined reference (REF_FA_ORIGINAL) or you also supply a folder with fasta references
(REF_MAPPING), which are used to map reads of contamination/non target organism. If the folder does not exists, the default method will be used.
Afterwards the remaining reads will be mapped to the target reference.


4. Variant calling

The mapped reads are used for variant calling, except for cortex, which uses the decomposed (unmapped) reads.

Command:
bash 003_A_variant_calling.sh


5. Collecting and filtering variants

The output variant files of different callers will be converted to and merged as a vcf file. Different genome and calling properties will be used 
to tag and filter the vcf files. Its output will be in vcfs_raw (intermediate) and vcfs (final).
The vcfs folder will contain 3 vcf files, the raw one, one prefiltered one with tags applied and a fully filtered one, which has been filtered according to the tags. 
There will be also 2 summary plots (pdf) about the SNP abundances and frequencies within the vcfs folder. This way you can get a quick overview of
your variant frequency, type, location and abundance. The filtered vcf file should be examined first. If any discrepancy is observed, then the prefiltered/raw file
can be consulted in order to further inspect that area.

Command:
bash 004_A_generate_vcfs.sh

Output files in vcfs folder:
Vcfs:
Fully filtered file: *filter_2tags.vcf
Filtered only by frequency: *filter_2.vcf
Raw file: *cov.vcf
Pdfs:
Coverage plot: *coverages.pdf
Frequency plot: *frequency.pdf


TAGS:
COL	(coverage low)	coverage 20% below average coverage of genome
COH	(coverage high)	coverage 20% above average coverage of genome
REP	(repetitive region)	variant is within a region that is longer than the insert size and is present at least 2 times with an edit distance of 4, usually filtered
HOP	(homopolymer region)	variant is within a homopolymeric region >= 8, usually filtered
SAR	(SNP accumulating region)	region with more than 5 SNPs within the insert size
BPA	(break position associated)	variant appears next to a break position, and could therefore be a FP
CV2	(caller per variant)	counts how many caller predicted SNPs or small InDels, CV1 usually filtered
CSV2	(caller per structural variant)	counts how many calls predicted large InDels, Inversions, Duplications, Transpositions (includes also different calls from the same tool)
CNSV2	(caller number per structural variant)	counts how many caller predicted large InDels, Inversions, Duplications, Transpositions (counts only different tools)
SVID	(structural variant id)	assignes every structural variant position a unique id


TODO: 406_annotate2vcf.sh


8. Conversion of vcf to ccf (tab delimited file) and mysql import

The following scripts convert the vcf file to column (tab separated) format and add information about sample and experiment (timepoint,condition,name...)
to be used to upload it to an existing mysql database.

bash 407_vcf2ccf.sh
bash 408_ccf2mysql.sh


9. Cleanup and copy results to project dir

TODO: 50_cleanup.sh
