----------------------------------------------
rnaseqmut: detecting mutations in RNA-Seq samples
Author: Wei Li 
Department of Biostatistics and Computational Biology
Dana-Farber Cancer Institute, Harvard School af Public Health
Email: li.david.wei AT gmail.com
----------------------------------------------

0. Introduction

rnaseqmut is a light-weight C++ program to detect variants (or mutations, including SNPs, indels) from RNA-Seq BAM files. It offers the following features:

  (1) Perform de-novo mutation discovery from a given BAM file;
  (2) For a user-defined mutation list, calculate the read coverage, including reads that support reference allele (reference reads) and alternative allele (alternative reads), from a given BAM file; 
  (3) For a series of RNA-Seq samples, filter interesting mutations based on user-defined criteria.

This software package includes a "core" C++ program, rnaseqmut, to call variants from a single BAM file, and a series of optional scripts, written in Python and bash, to identify putative interesting mutations in a group of RNA-Seq samples. 

Besides mutation detection from RNA-Seq, the "core" program (rnaseqmut) can also be used to call mutations from other high-throughput sequencing platforms, including ChIP-seq, DNA-Seq, etc.

I.     Before running rnaseqmut

II.    Usage

III.   Demo: detecting mutations from a series of RNA-Seq BAM files

IV.    Interpreting results

V.     Acknowledgements

VI.    Version History

----------------------------------------------
I.     Before running rnaseqmut
----------------------------------------------

1. System requirements

For the "core" program, you can either use the provided binary directly (for Linux and Mac OS 64 bit system), or compile the program from the source code. To use the binary directly, rename the corresponding executable to "rnaseqmut". For example, if you are Mac OS 64 bit user, go into the bin directory, and type:

  mv rnaseqmut.macos.x64 rnaseqmut

and refer to I.3 for installation. 

To compile the "core" program, only C++ compiler (gcc version >4.1) is needed. 

To use the supplementary scripts including demo, Linux system with Python (>=2.7) support is required.


2. Compiling

To compile rnaseqmut, you need to first compile bamtools. To compile bamtools, go to the src/bamtools dir, and follow the instructions in bamtools help page below:

https://github.com/pezmaster31/bamtools/wiki/Building-and-installing

After bamtools compilation is complete, go back to the src dir and type

  make

to finish the compilation. After compilation, the executable is in the bin/ directory.


3. Installation

After compiling is successful, add the "bin/" (or "scripts/", if necessary) directory to your PATH variable:

  export PATH=$PATH:/the_path_to_the_rnaseqmut_dir/bin:$PATH

4. Data preparation

rnaseqmut requires sorted BAM files with index (.bai file) as input. Use "samtools sort/index" to sort/index BAM files before running rnaseqmut.

rnaseqmut detects mutations from the NM tag in the BAM format (see http://samtools.sourceforge.net/SAMv1.pdf for BAM specification). By default, if the reference genome is provided, it will use it to infer the mutation (only SNPs, no indels). If the reference genome is not provided (or if the -t/--use_mdtag option is set), it will use the MD tag in the BAM format to infer variants; in this case the MD tag is required in the provided BAM format. 

If the NM tag is missing in your BAM file (this is common if reads are aligned using older RNA-Seq read mapping softwares), you may need to realign them using the latest RNA-Seq read aligners (e.g., Tophat).

If the MD tag is missing (like BAM files generated from STAR), you can use samtools to recalculate MD tag. The command line is as follows:

samtools calmd -b in.bam ref.fa > out.bam

where in.bam is the original BAM file, ref.fa is the reference sequence, and out.bam is the new BAM file with MD tag added.


----------------------------------------------
II.    Usage
----------------------------------------------

The usage of the "core" program, rnaseqmut, and a couple of optional scripts, is listed below.

1. rnaseqmut: the core mutation detection program
USAGE: 

   ./rnaseqmut  [-t] [-r <ref_fasta>] [-l <mutation_list>] [-k] [-d] [-s
                <max_mismatch>] [-n] [-i <min_read>] [-m <mut_span>] [--]
                [--version] [-h] <bam_file>


Where: 

   -t,  --use_mdtag
     Use MD Tag to call mutations instead of using reference genome (by
     -r/--ref_fasta option). This option is automatically set if the
     reference genome is not provided, and requires the BAM file contains
     the MD tag.

   -r <ref_fasta>,  --ref_fasta <ref_fasta>
     The (optional) fasta file for the reference genome. When this option
     is set, -d/--with_indel option will be ignored.

   -l <mutation_list>,  --mutation_list <mutation_list>
     The text file of a given, sorted list of mutations. Each line in a
     file records one mutations, with chromosome, location, reference and
     alternative sequence (separated by tab). The output will only include
     mutations within a given mutation list.

   -k,  --with_indel_read
     Do not skip reads with indels. By default all reads with indels are
     skipped as most RNA-Seq are performed by Illumina sequencing, which is
     prone to indel errors.

   -d,  --with_indel
     Do not skip indels. By default all indels are skipped as most RNA-Seq
     are performed by Illumina sequencing, which is prone to indel errors.
     This option will be ignored if -r/--ref_fasta option is provided.

   -s <max_mismatch>,  --max_mismatch <max_mismatch>
     The maximum number of mismatches in a read. Reads with more number of
     mismatches will be discarded. Default 1.

   -n,  --output_n
     Treat the character N as substitutions.

   -i <min_read>,  --min_read <min_read>
     The minimum read count for the mutation to output. Default 1.

   -m <mut_span>,  --mut_span <mut_span>
     The minimum distance of the mutation to the beginning (end) of the
     read. Default 4.

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <bam_file>
     (required)  The bam file from which mutation will be called


2. merge2ndvcf.py: (optional) merging mutations from multiple RNA-Seq scans.

usage: merge2ndvcf.py [-h] [-l LABEL] [-x REGION] [-m] [-v] [-r MIN_READ] ...

merge tab delimited files from multiple RNA-Seq samples (usually be a second
scan results of rnaseqmut program) into a tab delimited file, preserving all
information. All merged files should have exactly the same order and content
of mutations.

positional arguments:
  vcffiles              VCF file names

optional arguments:
  -h, --help            show this help message and exit
  -l LABEL, --label LABEL
                        The labels of each sample, separated by the comma
  -x REGION, --region REGION
                        Only output mutations falling into a specific region,
                        for example chr11:1-10000
  -m, --merged          Merge forward and reverse fields into one
  -v, --output-vcf      Output VCF formatted files
  -r MIN_READ, --min-read MIN_READ
                        Minimum read requirement for ALT reads. Default 4.

3. filtermut.py: (optional) filter mutations based on user-defined criterion. 

usage: filtermut.py [-h] [-c CONTROL] [-i TREATMENT] [-t MIN_RECURRENT]
                    [-d MIN_RECREAD] [-f MIN_RECFRAC] [-a MIN_REF]
                    [-b MAX_ALT] [-p] [-n] [-z] [-l LABELS] [--DP2]
                    [--DP2-out] [-x REGION]

filter mutations in samples.

optional arguments:
  -h, --help            show this help message and exit

Sample definitions:
  -c CONTROL, --control CONTROL
                        The index of control group samples, separated by
                        comma. For example, 0,2,4 defines a control group of 3
                        samples: the 1st, the 3rd and the 5th sample in the
                        order of the input table. Default: empty (do not use
                        any control samples)
  -i TREATMENT, --treatment TREATMENT
                        The index of treatment group samples, separated by
                        comma. Default: complement of control samples (if
                        -c/--control option is not specified, use all samples
                        as treatment samples).
  -t MIN_RECURRENT, --min-recurrent MIN_RECURRENT
                        Print mutations only occuring in at least this number
                        of good treatment samples, defined as those with
                        mutation >=min-recread reads and >=min-recfrac percent
                        frequency. Default 1.
  -d MIN_RECREAD, --min-recread MIN_RECREAD
                        Minimum alt reads defined in treated good samples,
                        default 10.
  -f MIN_RECFRAC, --min-recfrac MIN_RECFRAC
                        Minimum alt reads frequency in treated good samples,
                        default 0.2.
  -a MIN_REF, --min-ref MIN_REF
                        Minimum reference reads in control samples. Default 4.
  -b MAX_ALT, --max-alt MAX_ALT
                        Maximum alternative reads in control samples. Default
                        4.

Input/output options:
  -p, --passall         Do not do any filtering
  -n, --no-header       Do not print header and script used
  -z, --no-vcf          Do not print in vcf format; print as it is.
  -l LABELS, --labels LABELS
                        Labels used for each sample, separated by comma.
                        Default: SAMPLE_x where x is the sample ID. The number
                        of lables MUST be identical to the number of all
                        samples in the original table, not only those defined
                        by the -c/-i parameter.
  --DP2                 DP2 field instad of DP4 field is used in both
                        input/output files
  --DP2-out             DP2 field instad of DP4 field is used in output files.
                        This option is automatically set true if --DP2 is
                        specified.
  -x REGION, --region REGION
                        Only output mutations falling into a specific region,
                        for example chr11:1-10000



----------------------------------------------
III.   Demo: detecting mutations from a series of RNA-Seq BAM files
----------------------------------------------
A demo script is provided in the demo directory, together with 4 sample RNA-Seq BAM files. This demo illustrates the basic usage of rnaseqmut, including calling de-novo mutations, merging mutations from different samples, calling mutations based on a given list of mutations, and filtering mutations.

A Linux system (with Python support) is required to run the full demo. If only core program is needed, you can just run the first 3 or 4 steps. 

In this demo, we provide two normal samples and two tumor samples, and would like to search for mutations that only occurs in one or more tumor samples, but not (or have low frequency) in any of the normal samples. We split this job into five steps:

Step 1, scan the samples individually and get the mutation list for each sample. In this demo we output all possible mutations (those even occur in only 1 RNA-Seq read), but in reality you may just need those with enough read support (controlled by -i/--min_read option).

Step 2, merge the mutation lists from individual samples. This will be the potentially interesting mutations we would like to investigate.

Step 3, scan the samples again, using the provided mutation list in Step 2.

Step 4, merge the mutations in step 3 into a big table. In this table, each row represents a mutation, and columns record the number of supporting reads (and the number of reference reads that do not support this mutation) in all samples. Based on this table, you can use your own criteria to search interesting mutations.

Step 5 illustrates an example of filtering interesting mutations using a python script "filtermut.py". In this demo, we define "control" samples as two normal samples, and would like to look for mutations that satisfy the ALL of the following criteria:

--occur in at least 1 non-control sample (controlled by -t/--min-recurrent option in filtermut.py) with at least 20% frequency (controlled by -f/--min-recfrac option) and 10 supporting reads (-d/--min-recread);
--do not have any supporting reads in control samples (-b/--max-alt); 
--have at least 4 non-supporting reads (reference reads) in control samples (-a/--min-ref). This requirement will exclude mutations with 0 read coverage in control samples thus we have no idea whether these mutations occur in control samples.

The final output is a VCF file recording mutations and read coverages in all 4 samples. You can also output tab-delimited file instead of VCF file for further downstream analysis (use -z/--no-vcf option in Step 5). For interpreting results, see the next section. 

For more details, refer to the demo script and the usage of each programs/scripts.


----------------------------------------------
IV.    Interpreting results
----------------------------------------------
rnaseqmut will send a tab-delimited file, each line representing one possible mutation. For example,

chr1    4776457 A       C       135     56      1       0

This mutation occurs at chr1:4776457, and it is A to C mutation. The next 4 numbers are the number of reference reads (forward), reference reads (backward), alternative reads (forward) and alternative reads (backward).

If you merge mutations in N samples (as in Step 4 in the demo), there will be 4*N numeric columns following the first 4 column with the exact meaning as in single sample. The order is determined by the order of the input files.

For the script "filtermut.py" in demo Step 5, a VCF file format is provided. We use the "DP4" field in standard VCF file format to record the 4 numbers in each sample. For example,

chr1    4782642 .       T       G       1.0     NORMAL1.DP4=266,217,0,0;NORMAL2.DP4=264,215,0,0;TUMOR1.DP4=390,251,17,1;TUMOR2.DP4=249,184,15,0;



----------------------------------------------
V.    Acknowledgements
----------------------------------------------

rnaseqmut depends on several third-party C++ packages, including bamtools (written by Derek Barnett, https://github.com/pezmaster31/bamtools), TCLAP (Templatized C++ Command Line Parser Library, written by Michael E. Smoot, http://tclap.sourceforge.net/), FASTA/FASTQ parser in C (by Heng Li, http://lh3lh3.users.sourceforge.net/parsefastq.shtml). The source codes of these packages are included in rnaseqmut source code.

We thank Chenfei Wang and Robert K. Bradley for their help and feedback.

----------------------------------------------
VI.     Version History
----------------------------------------------
07/22/2020	Update bamtools to newer version; fix a bug to cause segmentation fault.
06/01/2016      rnaseqmut now supports BAM files generated by STAR.
                We now have a solution to handle BAM files without MD tag (like BAM files from STAR). See the manual.

06/09/2015  0.7 Fix a bug of calculating CIGAR strings.
                Fix a bug of reporting indel reads.

03/08/2014	0.6 Add a "-s" option to control the number of mutations allowed in a read.
                For the script, switch from python3 to python 2.7
                Add a binary executable for Mac OS X64 system

12/11/2013	0.5	Fix compilation option so running the program does not require dynamic linking of libbamtools library.
                Fix a bug that some mutations in the provided list are not printed due to the missing chromosome name in BAM file.
                A linux X64 binary is provided.

10/12/2013	0.4.1	Fix a bug that causes the program to quit due to the index of the reference genome.
			          Suppress error messages due to the incompatible CIGAR letter.

09/16/2013	0.4	For BAM files with MD tag missing, rnaseqmut now allows inferring mutations directly from the given reference genome. If reference genome is given, rnaseqmut will use it (instead of MD tag) to call mutations.

09/16/2013	0.3	rnaseqmut allows comparing mutations with reference genome.
                By default, rnaseqmut will skip all indels and all reads with indels; this behavior can be controlled by -k and -d options. 
                Fix several bugs.

09/14/2013	0.2	Fix a bug of failing to read NM tags for some versions of BAM.
                By default, suppress mutations with character "N" (controlled by -n option in rnaseqmut)

09/12/2013	0.1	rnaseqmut software initiated.


