Notes for BisuKit v0.1.2
-----------------------------------------------------------

- Updated 2014-03-24 : Details on Options, Testing BisuKit
- Updated 2014-03-25 : BisuKit code upload
- Updated 2014-05-12 : Code updated, README updated
- Updated 2014-05-13 : Tarball added
- Updated 2014-06-11 : Future updates part of README updated
- Updated 2016-02-01 : Major updates to script (simplified the pipeline) and same version runningon DE of Cyverse
- Updated 2016-02-24 : Major bug fix with RPy2 in BisuKit crashing when chromosome.methylKit file empty
- Updated 2016-06-15 : Major update to comply with the version available on CyVerse
- Updated 2017-01-10 : Citation updated

1. Summary 
==========

BisKit has following features:

- Use Bismark's methylation extractor to get statistics of methylated cytosines
- Identify Differentially Methylated Regions (DMRs) using R packages methylKit and eDMR

Bug reports and comments are highly appreciated.

2. Building 
===========
You can run this directly from command line using python bisukit.py or using CyVerse Discovery Environment (Reference: TBD)
Currently the linux x86_64 platform is supported.  
Email me if you need to run BisuKit on other platforms.

BisuKit requires Python 2.7.9

These are some of the packages that BisKit requires (links and version to be updated) :
- RPy2 for Python
- R libraries (big.memory, GenomicRanges, data.table, ggplot2)

The files you need in path are not actually BAM files but Context specific methylation extracted files
http://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide_v0.14.4.pdf 

Most errors that occur at this stage are related to the PATH settings.
Check that all of these softwares are defined in PATH if error persists.

3. Running BisuKit
==================
Use: python bisukit.py 

Options:
- -h, --help            show this help message and exit

- --genome=GENOME     name and directory of genome fasta file
- --context=CONTEXT   CpG, CHG or CHH
- --specie=SPECIE     Specie code, B73, MM9, HG19
- --cores=CORES       Number of cores available
- -q, --quiet         don't print status messages to stdout

- --name1=NAME1 Name for sample 1
- --name2=NAME2 Name for sample 2

- --ot1=OT1     path of first Top strand methylation extractor file generated from bismark methylation extractor
- --ot2=OT2     path of second Top strand methylation extractor file generated from bismark methylation extractor
- --ob1=OB1     path of first Bottom strand methylation extractor file generated from bismark methylation extractor
- --ob2=OB2     path of second Bottom strand methylation extractor file generated from bismark methylation extractor

- --qvalue=0.01 Q-value of DMRs to filter from all DMRs generated
- --cpg=3       Minimum number of CpGs in a DMR (default to 3)
- --dmc=1       Minimum number of Differentially Methylated Cytosines (DMCs) in a DMR (default to 1)
- --diffmeth=20 Minimum percentage of methylation difference between two samples in a DMR (default to 20)

and if available
- --ctot1=CTOT1     path of first complimentary Top strand methylation extractor file generated from bismark methylation extractor
- --ctot2=CTOT2     path of second complimentary Top strand methylation extractor file generated from bismark methylation extractor
- --ctob1=CTOB1     path of first complimentary Bottom strand methylation extractor file generated from bismark methylation extractor
- --ctob2=CTOB2     path of second complimentary Bottom strand methylation extractor file generated from bismark methylation extractor

3a. Running ZedtoBisuKit
==================
Use: python zedtobisukit.py 

The only difference between BisuKit and ZedtoBisuKit is the input files (OT1,OT2,OB1,OB2,CTOT1,CTOT2,CTOB1,CTOB2 for BisuKit)

- --zedmethratio1=ZEDMETHRATIO1 Output of sample 1 from running Zed-align methratio.py (in format of $name_methratio.txt)
- --zedmethratio2=ZEDMETHRATIO2 Output of sample 2 from running Zed-align methratio.py (in format of $name_methratio.txt)

4. Testing BisuKit
==================
To check whether you have all Python dependencies in place,
```
python check_dependencies.py
```

Typical BisuKit analysis,
```
python bisukit.py --name1 B73_1 --name2 B73_2 --ot1 CpG_OT_293424_GGTAGC_L005_R1_001.fastq_bismark_bt2_pe.txt --ot2 CpG_OT_293426_GTGCTA_L005_R1_001.fastq_bismark_bt2_pe.txt --ob1 CpG_OB_293424_GGTAGC_L005_R1_001.fastq_bismark_bt2_pe.txt --ob2 CpG_OB_293426_GTGCTA_L005_R1_001.fastq_bismark_bt2_pe.txt --genome Zea_mays.AGPv3.23.fa --qvalue 0.01 --dmc 1 --cpg 3 --diffmeth 20 --cores 16 --specie B73 --context CpG 
```

Typical ZedtoBisuKit analysis,
```
python zedtobisukit.py --name1 B73_1 --name2 B73_2 --zedmethratio1 B73_test1_methratio.txt --zedmethratio2 B73_test2_methratio.txt --genome Zea_mays.AGPv3.23_without_scaffold.fa --qvalue 0.01 --dmc 1 --cpg 3 --diffmeth 20 --cores 16 --specie B73 --context CpG
```

5. License
===========

- This is a open source package.
- The program itself may not be modified in any way.
- This license does not allow the use of this program for any commercial purpose. 
- No guarantees are given as to the program's correctness, or the accuracy or completeness of its output.  
- The author accepts no liability for damage or otherwise following from using and interpreting the output of this program.
- The software is supplied "as is", without obligation by the author to provide any services or support.

6. Revision history
====================

- 0.2 (10 January 2017)
- 0.1.2 (15 June 2016)
- 0.1.2 (1 February 2016)
- 0.1a (10 February 2014)  
        - first release

7. FAQ
=======
Q: How do I cite BisuKit (Zedtobisukit)? 
A: Song, J., Zynda, G., Beck, S., Springer, N.M., and Vaughn, M.W. 2016. Bisulfite sequence analyses using CyVerse discovery environment: From mapping to DMRs. Curr. Protoc. Plant Biol. 1:510-529. doi: 10.1002/cppb.20034

8. Future updates
=======
- Add more graphical output .bigWig or .circos
- Add GSNAP compatibility
