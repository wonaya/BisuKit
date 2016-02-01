Notes for BisuKit v0.1.2
-----------------------------------------------------------

- Updated 2014-03-24 : Details on Options, Testing BisuKit
- Updated 2014-03-25 : BisuKit code upload
- Updated 2014-05-12 : Code updated, README updated
- Updated 2014-05-13 : Tarball added
- Updated 2014-06-11 : Future updates part of README updated
- Updated 2016-02-01 : Major updates to script (simplified the pipeline) and same version runningon DE of Cyverse

1. Summary 
==========

BisKit has following features:

- Use Bismark's methylation extractor to get statistics of methylated cytosines
- Identify Differentially Methylated Regions (DMRs) using R packages methylKit and eDMR

Bug reports and comments are highly appreciated.

2. Building 
===========
```
tar -zxvf bisukit-0.1.2.tar.gz
```

Currently the linux x86_64 platform is supported.  
Email me if you need to run BisuKit on other platforms.

BisuKit requires Python (version check)

These are some of the packages that BisKit requires (links and version to be updated) :
- RPy2 for Python
- R

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
- --cores=CORES       no. of cores to use in running BisKit
- -q, --quiet         don't print status messages to stdout

- --bam1=BAMFILE1     name and directory of first bam file generated from bismark for DMR finding
- --bam2=BAMFILE2     name and directory of second bam file generated from bismark for DMR finding

4. Testing BisuKit
==================
To check whether you have all Python dependencies in place,
```
python check_dependencies.py
```

Typical methylKit/eDMR analysis,
```
python bisukit.py --genome HG19.fa --context CpG --cores 12 --bam1 SRR641618_unsorted.bam --bam2 SRR641628_unsorted.bam --specie HG19
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

- 0.1.2 (1 February 2016)
- 0.1a (10 February 2014)  
        - first release

7. FAQ
=======


8. Future updates
=======

