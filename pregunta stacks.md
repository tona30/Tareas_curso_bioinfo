**4.1 Clean the data**

In a typical analysis, data will be received from an Illumina sequencer, or some other type of sequencer as FASTQ files. The first requirement is to demultiplex, or sort, the raw data to recover the individual samples in the Illumina library. While doing this, we will use the Phred scores provided in the FASTQ files to discard sequencing reads of low quality. These tasks are accomplished using the process_radtags program. 

![](/home/tona/Escritorio/TAREA 4/process_radtags1.png) 

 Some things to consider when running this program:

 - process_radtags can handle both single-end or paired-end Illumina sequencing.
 - The raw data can be compressed, or gzipped (files end with a "`.gz`" suffix).
- You can supply a list of barcodes, or indexes, to process_radtags in order for it to demultiplex your samples. These barcodes can be single-end barcodes or combinatorial barcodes (pairs of barcodes, one on each of the paired reads). Barcodes are specified, one per line (or in tab separated pairs per line), in a text file.
      - If, in addition to your barcodes, you also supply a sample name in an extra column within the barcodes file, process_radtags will name your output files according to sample name instead of barcode.
 - If you believe your reads may contain adapter contamination, process_radtags can filter it out.
 - You can supply the restriction enzyme used to construct the library. In the case of double-digest RAD, you can supply both restriction enzymes.
- If instructed, (`-r` command line option), process_radtags will correct barcodes and restriction enzyme sites that are within a certain distance from the true barcode or restriction enzyme cutsite.

**4.1.1 Understanding barcodes/indexes and specifying the barcode type**

Genotype by sequencing libraries sample the genome by selecting DNA adjacent to one or more restriction enzyme cutsites. By reducing the amount of total DNA sampled, most researchers will multiplex many samples into one molecular library. Individual samples are demarcated in the library by ligating an oligo barcode onto the restriction enzyme-associated DNA for each sample. Alternatively, an index barcode is used, where the barcode is located upstream of the sample DNA within the sequencing adaptor. Regardless of the type of barcode used, after sequencing, the data must be demultiplexed so the samples are again separated. The `process_radtags` program will perform this task, but we much specify the type of barcodes used, and where to find them in the sequencing data. 

 There are a number of different configurations possible, each of them is detailed below.

1.If your data are single-end or paired-end, with an inline barcode present only on the single-end (marked in red):

    @HWI-ST0747:188:C09HWACXX:1:1101:2968:2083 1:N:0: TTATGATGCAGGACCAGGATGACGTCAGCACAGTGCGGGTCCTCCATGGATGCTCCTCGGTCGTGGTTGGGGGAGGAGGCA 
    +     
    @@@DDDDDBHHFBF@CCAGEHHHBFGIIFGIIGIEDBBGFHCGIIGAEEEDCC;A?;;5,:@A?=B5559999B@BBBBBA @HWI-ST0747:188:C09HWACXX:1:1101:2863:2096 1:N:0: TTATGATGCAGGCAAATAGAGTTGGATTTTGTGTCAGTAGGCGGTTAATCCCATACAATTTTACACTTTATTCAAGGTGGA 
    +
    CCCFFFFFHHHHHJJGHIGGAHHIIGGIIJDHIGCEGHIFIJIH7DGIIIAHIJGEDHIDEHJJHFEEECEFEFFDECDDD @HWI-ST0747:188:C09HWACXX:1:1101:2837:2098 1:N:0: GTGCCTTGCAGGCAATTAAGTTAGCCGAGATTAAGCGAAGGTTGAAAATGTCGGATGGAGTCCGGCAGCAGCGAATGTAAA
    
   Then you can specify the `--inline_null` flag to process_radtags. This is also the default behavior and the flag can be ommitted in this case. 
   
2.If your data are single-end or paired-end, with a single index barcode (in blue):

    @9432NS1:54:C1K8JACXX:8:1101:6912:1869 1:N:0:ATGACT TCAGGCATGCTTTCGACTATTATTGCATCAATGTTCTTTGCGTAATCAGCTACAATATCAGGTAATATCAGGCGCA 
    + 
    CCCFFFFFHHHHHJJJJJJJJIJJJJJJJJJJJHIIJJJJJJIJJJJJJJJJJJJJJJJJJJGIJJJJJJJHHHFF @9432NS1:54:C1K8JACXX:8:1101:6822:1873 1:N:0:ATGACT CAGCGCATGAGCTAATGTATGTTTTACATTCCAGAAAGAGAGCTACTGCTGCAGGTTGTGATAAAATAAAGTAAGA 
    + 
    B@@FFFFFHFFHHJJJJFHIJHGGGHIJIIJIJCHJIIGGIIIGGIJEHIJJHII?FFHICHFFGGHIIGG@DEHH @9432NS1:54:C1K8JACXX:8:1101:6793:1916 1:N:0:ATGACT TTTCGCATGCCCTATCCTTTTATCACTCTGTCATTCAGTGTGGCAGCGGCCATAGTGTATGGCGTACTAAGCGAAA 
    +
    @C@DFFFFHGHHHGIGHHJJJJJJJGIJIJJIGIJJJJHIGGGHGII@GEHIGGHDHEHIHD6?493;AAA?;=;=
 Then you can specify the `--index_null` flag to process_radtags.
 3.If your data are single-end with both an inline barcode (in red) and an index barcode (in blue): 
 
    @9432NS1:54:C1K8JACXX:8:1101:6912:1869 1:N:0:ATCACG TCACGCATGCTTTCGACTATTATTGCATCAATGTTCTTTGCGTAATCAGCTACAATATCAGGTAATATCAGGCGCA 
    + 
    CCCFFFFFHHHHHJJJJJJJJIJJJJJJJJJJJHIIJJJJJJIJJJJJJJJJJJJJJJJJJJGIJJJJJJJHHHFF @9432NS1:54:C1K8JACXX:8:1101:6822:1873 1:N:0:ATCACG GTCCGCATGAGCTAATGTATGTTTTACATTCCAGAAAGAGAGCTACTGCTGCAGGTTGTGATAAAATAAAGTAAGA 
    + 
    B@@FFFFFHFFHHJJJJFHIJHGGGHIJIIJIJCHJIIGGIIIGGIJEHIJJHII?FFHICHFFGGHIIGG@DEHH @9432NS1:54:C1K8JACXX:8:1101:6793:1916 1:N:0:ATCACG GTCCGCATGCCCTATCCTTTTATCACTCTGTCATTCAGTGTGGCAGCGGCCATAGTGTATGGCGTACTAAGCGAAA 
    +
    @C@DFFFFHGHHHGIGHHJJJJJJJGIJIJJIGIJJJJHIGGGHGII@GEHIGGHDHEHIHD6?493;AAA?;=;= 
 Then you can specify the `--inline_index` flag to process_radtags. 
 
 4.If your data are paired-end with an inline barcode on the single-end (in red) and an index barcode (in blue):
 
    @9432NS1:54:C1K8JACXX:7:1101:5584:1725 1:N:0:CGATGT ACTGGCATGATGATCATAGTATAACGTGGGATACATATGCCTAAGGCTAAAGATGCCTTGAAGCTTGGCTTATGTT 
    + 
    #1=DDDFFHFHFHIFGIEHIEHGIIHFFHICGGGIIIIIIIIAEIGIGHAHIEGHHIHIIGFFFGGIIIGIIIEE7 @9432NS1:54:C1K8JACXX:7:1101:5708:1737 1:N:0:CGATGT TTCGACATGTGTTTACAACGCGAACGGACAAAGCATTGAAAATCCTTGTTTTGGTTTCGTTACTCTCTCCTAGCAT 
    +
    #1=DFFFFHHHHHJJJJJJJJJJJJJJJJJIIJIJJJJJJJJJJIIJJHHHHHFEFEEDDDDDDDDDDDDDDDDD@
    
    @9432NS1:54:C1K8JACXX:7:1101:5584:1725 2:N:0:CGATGT AATTTACTTTGATAGAAGAACAACATAAGCCAAGCTTCAAGGCATCTTTAGCCTTAGGCATATGTATCCCACGTTA 
    +
    @@@DFFFFHGHDHIIJJJGGIIIEJJJCHIIIGIJGGEGGIIGGGIJIJIHIIJJJJIJJJIIIGGIIJJJIICEH @9432NS1:54:C1K8JACXX:7:1101:5708:1737 2:N:0:CGATGT AGTCTTGTGAAAAACGAAATCTTCCAAAATGCTAGGAGAGAGTAACGAAACCAAAACAAGGATTTTCAATGCTTTG 
    + 
    C@CFFFFFHHHHHJJJJJJIJJJJJJJJJJJJJJIJJJHIJJFHIIJJJJIIJJJJJJJJJHGHHHHFFFFFFFE
    
   Then you can specify the `--inline_index` flag to process_radtags.
   
   5.If your data are paired-end with indexed barcodes on the single and paired-ends (in blue): 
   
    @9432NS1:54:C1K8JACXX:7:1101:5584:1725 1:N:0:ATCACG+CGATGT ACTGGCATGATGATCATAGTATAACGTGGGATACATATGCCTAAGGCTAAAGATGCCTTGAAGCTTGGCTTATGTT 
    +
    #1=DDDFFHFHFHIFGIEHIEHGIIHFFHICGGGIIIIIIIIAEIGIGHAHIEGHHIHIIGFFFGGIIIGIIIEE7 @9432NS1:54:C1K8JACXX:7:1101:5708:1737 1:N:0:ATCACG+CGATGT TTCGACATGTGTTTACAACGCGAACGGACAAAGCATTGAAAATCCTTGTTTTGGTTTCGTTACTCTCTCCTAGCAT 
    +
     #1=DFFFFHHHHHJJJJJJJJJJJJJJJJJIIJIJJJJJJJJJJIIJJHHHHHFEFEEDDDDDDDDDDDDDDDDD@
     
    @9432NS1:54:C1K8JACXX:7:1101:5584:1725 2:N:0:ATCACG+CGATGT AATTTACTTTGATAGAAGAACAACATAAGCCAAGCTTCAAGGCATCTTTAGCCTTAGGCATATGTATCCCACGTTA 
    +
     @@@DFFFFHGHDHIIJJJGGIIIEJJJCHIIIGIJGGEGGIIGGGIJIJIHIIJJJJIJJJIIIGGIIJJJIICEH @9432NS1:54:C1K8JACXX:7:1101:5708:1737 2:N:0:ATCACG+CGATGT AGTCTTGTGAAAAACGAAATCTTCCAAAATGCTAGGAGAGAGTAACGAAACCAAAACAAGGATTTTCAATGCTTTG 
     + 
    C@CFFFFFHHHHHJJJJJJIJJJJJJJJJJJJJJIJJJHIJJFHIIJJJJIIJJJJJJJJJHGHHHHFFFFFFFED
    
   Then you can specify the `--index_index` flag to process_radtags. 
   
   6.If your data are **paired-end** with inline barcodes on the single and paired-ends (in red):
   
    @9432NS1:54:C1K8JACXX:7:1101:5584:1725 1:N:0: ACTGGCATGATGATCATAGTATAACGTGGGATACATATGCCTAAGGCTAAAGATGCCTTGAAGCTTGGCTTATGTT 
    +
     #1=DDDFFHFHFHIFGIEHIEHGIIHFFHICGGGIIIIIIIIAEIGIGHAHIEGHHIHIIGFFFGGIIIGIIIEE7 @9432NS1:54:C1K8JACXX:7:1101:5708:1737 1:N:0: TTCGACATGTGTTTACAACGCGAACGGACAAAGCATTGAAAATCCTTGTTTTGGTTTCGTTACTCTCTCCTAGCAT 
     +
    #1=DFFFFHHHHHJJJJJJJJJJJJJJJJJIIJIJJJJJJJJJJIIJJHHHHHFEFEEDDDDDDDDDDDDDDDDD@
    
    @9432NS1:54:C1K8JACXX:7:1101:5584:1725 2:N:0: AATTTACTTTGATAGAAGAACAACATAAGCCAAGCTTCAAGGCATCTTTAGCCTTAGGCATATGTATCCCACGTTA 
    +
    @@@DFFFFHGHDHIIJJJGGIIIEJJJCHIIIGIJGGEGGIIGGGIJIJIHIIJJJJIJJJIIIGGIIJJJIICEH @9432NS1:54:C1K8JACXX:7:1101:5708:1737 2:N:0: AGTCTTGTGAAAAACGAAATCTTCCAAAATGCTAGGAGAGAGTAACGAAACCAAAACAAGGATTTTCAATGCTTTG 
    + 
    C@CFFFFFHHHHHJJJJJJIJJJJJJJJJJJJJJIJJJHIJJFHIIJJJJIIJJJJJJJJJHGHHHHFFFFFFFED
    
   Then you can specify the `--inline_inline` flag to process_radtags. 
   
 **4.1.2Specifying the barcodes**

The process_radtags program will demultiplex data if it is told which barcodes/indexes to expect in the data. It will also properly name the output files if the user specifies how to translate a particular barcode to a specific output file neam. This is done with the barcodes file, which we provide to process_radtags program. The barcode file is a very simple format — one barcode per line; if you want to rename the output files, the sample name *prefix* is provided in the second column. 

    %cat barcodes_lane3 
    CGATA<tab>sample_01 
    CGGCG     sample_02 
    GAAGC     sample_03 
    GAGAT     sample_04 
    TAATG     sample_05 
    TAGCA     sample_06 
    AAGGG     sample_07 
    ACACG     sample_08 
    ACGTA     sample_09
    
   The sample names can be whatever is meaningful for your project:
   
    %more barcodes_run01_lane01 
    CGATA<tab>spruce_site_12-01 
    CGGCG spruce_site_12-02 
    GAAGC spruce_site_12-03 
    GAGAT spruce_site_12-04 
    TAATG spruce_site_06-01 
    TAGCA spruce_site_06-02 
    AAGGG spruce_site_06-03 
    ACACG spruce_site_06-04 
    
  Combinatorial barcodes are specified, one per column, separated by a tab: 
  
    cat barcodes_lane07 
    CGATA<tab>ACGTA<tab>sample_01 
    CGGCG     ACGTA     sample_02 
    GAAGC     ACGTA     sample_03 
    GAGAT     ACGTA     sample_04 
    CGATA     TAGCA     sample_05 
    CGGCG     TAGCA     sample_06 
    GAAGC     TAGCA     sample_07 
    GAGAT     TAGCA     sample_08
    
    
 1.If you don't want process_radtags to rename your samples, simply do not specify the last column in the barcodes file, and the output files will instead be named after the barcode.
 2.Often, sequencing centers will return data from indexed libraries already demultiplexed. In this case, omit the barcodes file and `process_radtags` will not attempt to demulitplex the data, but can still be used to clean the data.
 
**4.1.3Running process_radtags**
Here is how single-end data received from an Illumina sequencer might look: 

    % ls ./raw 
    lane3_NoIndex_L003_R1_001.fastq.gz 
    lane3_NoIndex_L003_R1_002.fastq.gz 
    lane3_NoIndex_L003_R1_003.fastq.gz 
    lane3_NoIndex_L003_R1_004.fastq.gz 
    lane3_NoIndex_L003_R1_005.fastq.gz 
    lane3_NoIndex_L003_R1_006.fastq.gz 
    lane3_NoIndex_L003_R1_007.fastq.gz 
    lane3_NoIndex_L003_R1_008.fastq.gz 
    lane3_NoIndex_L003_R1_009.fastq.gz 
    lane3_NoIndex_L003_R1_010.fastq.gz 
    lane3_NoIndex_L003_R1_011.fastq.gz 
    lane3_NoIndex_L003_R1_012.fastq.gz 
    lane3_NoIndex_L003_R1_013.fastq.gz
    
 Then you can run process_radtags in the following way: 
 
     % process_radtags -p ./raw/ -o ./samples/ -b ./barcodes/barcodes_lane3 
     \ -e sbfI -r -c -q
     
 I specify the directory containing the input files, `./raw`, the directory I want process_radtags to enter the output files, `./samples`, and a file containing my barcodes, `./barcodes/barcodes_lane3`, along with the restrction enzyme I used and instructions to clean the data and correct barcodes and restriction enzyme cutsites (`-r`, `-c`, `-q`).

Here is a more complex example, using paired-end double-digested data (two restriction enzymes) with combinatorial barcodes, and gzipped input files. Here is what the raw Illumina files may look like: 

     % ls ./raw 
     GfddRAD1_005_ATCACG_L007_R1_001.fastq.gz    
     GfddRAD1_005_ATCACG_L007_R2_001.fastq.gz 
     GfddRAD1_005_ATCACG_L007_R1_002.fastq.gz 
     GfddRAD1_005_ATCACG_L007_R2_002.fastq.gz 
     GfddRAD1_005_ATCACG_L007_R1_003.fastq.gz 
     GfddRAD1_005_ATCACG_L007_R2_003.fastq.gz 
     GfddRAD1_005_ATCACG_L007_R1_004.fastq.gz 
     GfddRAD1_005_ATCACG_L007_R2_004.fastq.gz 
     GfddRAD1_005_ATCACG_L007_R1_005.fastq.gz 
     GfddRAD1_005_ATCACG_L007_R2_005.fastq.gz 
     GfddRAD1_005_ATCACG_L007_R1_006.fastq.gz 
     GfddRAD1_005_ATCACG_L007_R2_006.fastq.gz 
     GfddRAD1_005_ATCACG_L007_R1_007.fastq.gz 
     GfddRAD1_005_ATCACG_L007_R2_007.fastq.gz 
     GfddRAD1_005_ATCACG_L007_R1_008.fastq.gz 
     GfddRAD1_005_ATCACG_L007_R2_008.fastq.gz 
     GfddRAD1_005_ATCACG_L007_R1_009.fastq.gz 
     GfddRAD1_005_ATCACG_L007_R2_009.fastq.gz
     
 Now we specify both restriction enzymes using the `--renz_1` and `--renz_2` flags along with the type combinatorial barcoding used. Here is the command:

     % process_radtags -P -p ./raw -b ./barcodes/barcodes -o ./samples/ \ -
     c -q -r --inline_index --renz_1 nlaIII --renz_2 mluCI
     
 **The output of process_radtags**

The output of the process_radtags differs depending if you are processing single-end or paired-end data. In the case of *single-end* reads, the program will output one file per barcode into the output directory you specify. If the data do not have barcodes, then the file will retain its original name.

If you are processing *paired-end reads*, then you will get four files per barcode, two for the single-end read and two for the paired-end read. For example, given barcode ACTCG, you would see the following four files:

     sample_ACTCG.1.fq
     sample_ACTCG.rem.1.fq
     sample_ACTCG.2.fq
     sample_ACTCG.rem.2.fq
     
 The process_radtags program wants to keep the reads in phase, so that the first read in the `sample_XXX.1.fq` file is the mate of the first read in the `sample_XXX.2.fq` file. Likewise for the second pair of reads being the second record in each of the two files and so on. When one read in a pair is discarded due to low quality or a missing restriction enzyme cut site, the remaining read can't simply be output to the `sample_XXX.1.fq` or `sample_XXX.2.fq` files as it would cause the remaining reads to fall out of phase. Instead, this read is considered a remainder read and is output into the `sample_XXX.rem.1.fq` file if the paired-end was discarded, or the `sample_XXX.rem.2.fq` file if the single-end was discarded. 
 
**Modifying how process_radtags executes**

The process_radtags program can be modified in several ways. If your data do not have barcodes, omit the barcodes file and the program will not try to demultiplex the data. You can also disable the checking of the restriction enzyme cut site, or modify what types of quality are checked for. So, the program can be modified to only demultiplex and not clean, clean but not demultiplex, or some combination.

There is additional information available in process_radtags manual page. 

**4.2 Align data against a reference genome**

If a reference genome is available for your organism, once you have demultiplexed and cleaned your samples, you will align these samples using a standard alignment program, such as GSnap, BWA, or Bowtie2.

Performing this alignment is outside the scope of this document, however there are several things you should keep in mind:

 - *Stacks *can read BAM and SAM files, so you can use any aligner that can generate these standard format output files.
 - *Stacks* can understand gapped alignments. It will interpret the alignment according to the CIGAR string in the SAM/BAM file.
 - Make sure you produce unique alignments, or randomly place repetitive alignments.
 - Be careful of **terminal alignments**. Several alignment programs will create terminal alignments: if an aligner cannot find a full-length alignment it will soft-mask the portion of the read it cannot align and report the fragment that does align. *Stacks *will convert these soft-masked regions to Ns, since they are not part of the alignment. If this occurs extensively at a locus, it can become impossible to call certain polymorphisms or to read the haplotypes from the locus (since the reads end up as random fragments at the locus, different parts of different reads cover the locus, making it impossible to properly read the haplotype). You can confirm that you have terminal alignments by looking at the CIGAR string (which describes the alignment of a read) in the SAM or BAM file. You should see soft-masked regions (an "S" in the CIGAR string) at the end of the read. This alignment behavior can be turned off in some aligners (GSnap).

**4.3 Run the pipeline**

The simplest way to run the pipeline is to use one of the two wrapper programs provided: if you do not have a reference genome you will use denovo_map.pl, and if you do have a reference genome use ref_map.pl. In each case you will specify a list of samples that you demultiplexed in the first step to the program, along with several command line options that control the internal algorithms.

There are two default experimental analyses that you can execute with the wrappers:

 - **A genetic map**: In this case you specify your samples in two groups to denovo_map.pl or ref_map.pl. First, one or two parents in a genetic cross are provided (`-p` command line option), and second, a list of progeny from the cross (`-r` command line option) are supplied.
 - **A population genomic analysis**: In this case you specify all your samples generically to denovo_map.pl or ref_map.pl (`-s` command line option). For sets of populations you can optionally supply a population map -- a file that specifies which samples belong to which populations (See below for more detail). If you don’t supply this file, Stacks will assume all of the specified samples represent a single population.
-  In a genetic map analysis, the catalog -- which is a population-wide view of all loci in the data set -- will be populated only from the parents provided to Stacks, in a population analysis, all individuals are loaded into the catalog.
 - If a set of parents and progeny are supplied to denovo_map.pl or ref_map.pl, both pipelines will execute the genotypes program as the final stage. This will identify mappable loci and generically encode the the type of each of those loci (e.g. an `ab x aa` locus would represent a locus with two alleles in one parent and a single allele in the second parent). The user can re-execute genotypes after the pipeline completes to get genotypes for a particular cross (e.g. F2 or a backcross) and for a particular linkage mapper (e.g. JoinMap or OneMap). Much more information, including examples, can be seen below on the genotypes manual page.
 - If instead a set of generic samples is provided to denovo_map.pl or ref_map.pl, both pipelines will execute the populations program as the final stage. This program will calculate population genomic statistics such as heterozygosity, π, and FIS. A population map can be supplied to populations to generate FST for pairs of populations. Again, populations can be re-executed after the pipeline runs with a variety of command line options to filter the data in a number of ways. In addition, if the loci are reference-aligned, a sliding window algorithm can be turned on to generate a kernel-smoothed average of statistics such as FIS and FST along the genome; a number of different output formats can also be generated. Much more information, including examples, can be seen below on the populations manual page.

