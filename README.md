Notes for BisuKit v0.1a
-----------------------------------------------------------

- Updated 2014-03-24 : Details on Options, Testing BisuKit
- Updated 2014-03-25 : BisuKit code upload

1. Summary 
==========

BisKit has following features:

- Use Bismark's methylation extractor to get statistics of methylated cytosines
- Generate tiles for bisulfite sequencing data
- Identify Differentially Methylated Regions (DMRs) using R packages methylKit and eDMR
- Generate configuration file for Circos plot generation

Bug reports and comments are highly appreciated.

2. Building 
===========
```
tar -zxvf BisuKit-0.1a.tar.gz
python setup.py
```

Currently the linux x86_64 platform is supported.  
Email me if you need to run BisuKit on other platforms.

BisuKit requires Python (version check)

These are some of the packages that BisKit requires (links and version to be updated) :
- Samtools
- Bedtools
- Bismark
- RPy2 for Python
- R

Most errors that occur at this stage are related to the PATH settings.
Check that all of these softwares are defined in PATH if error persists.

3. Running BisuKit
==================
Options:
- -h, --help            show this help message and exit

Required for all function:
- --run=RUN           type of run : tile / methylation_extractor / methylkit / circos
- --genome=GENOME     name and directory of genome fasta file
- --context=CONTEXT   CpG, CHH, CHG or all
- --specie=SPECIE     Specie code, B73, MM9, HG19
- --cores=CORES       no. of cores to use in running BisKit
- -q, --quiet         don't print status messages to stdout

Tile specific options:
These are tile function specific options
- --window=WINDOW     tiling window size (in bp)
- --separate_strands=STRANDED    merge forward and reverse strands into tile=no, output for separate strands=yes
- --path_to_samtools=SAMTOOLS    full path to samtools (eg. /usr/local/bin/samtools)
- --path_to_multibamcov=MULTIBAMCOV  full path to multiBamCov (eg. /usr/local/bin/bedtools/bin/multiBamCov)
- --sam=SAMFILE       name and directory of first sam file generated for tiling analysis, use unsorted sam file
- --circos=CIRCOS     Generate circos format file, default=no

Methylation extractor specific options:
These are Bismark Methylation extractor options
- --paired=PAIRED     paired-end or single-end reads, yes or no
- --path_to_bismark=BISMARK full path to bismark (eg. /usr/local/bin/bismark/bismark)

MethylKit specific options:
These are MethylKit options
- --sam1=SAMFILE1     name and directory of first sam file generated from bismark for DMR finding
- --sam2=SAMFILE2     name and directory of second sam file generated from bismark for DMR finding
- --specie=SPECIE     Specie code, B73, MM9, HG19
- --largemem=LARGEMEM If using large memory notes = yes (tested on 1TB node)
- --nonb73=NONB73     For maize genome, if looking at nonB73 genotype (Maize specific)

4. Testing BisuKit
==================
To check whether you have all Python dependencies in place,
```
python check_dependencies.py
```

Typical tiling analysis, 
```
python bisukit.py --run tile --genome HG19.fa --context CpG --window 100 --paired no --cores 12 --sam SRR641618_unsorted.sam --largemem yes --circos yes
```

Typical methylKit/eDMR analysis,
```
python bisukit.py --run methylkit --genome HG19.fa --context CpG --paired no --cores 12 --sam1 SRR641618_unsorted.sam --sam2 SRR641628_unsorted.sam --specie HG19 --largemem yes
```

To run these analyses, you need methylation extracted files. To run Bismark's methylkation extractor
```
python bisukit.py --run methylkation_extractor --path_to_bismark /usr/local/bin/bismark/bismark --sam1 SRR641618_unsorted.sam --sam2 SRR641628_unsorted.sam --paired no --cores 2 
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

- 0.1a (10 February 2014)  
        - first release

7. FAQ
=======

To Be Updated

8. Future updates
=======

To Be Updated
