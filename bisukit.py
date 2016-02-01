#!/usr/bin/env python

### FIX duplicate issue at methylKit level - Jan 25 done
### Test entire process - Jan 28
### tmp into random number generate - Jan 28
### Make into Agave App - Jan 29

import os,sys
import resource
from optparse import OptionParser,OptionGroup
from datetime import datetime
try:
   import subprocess
except ImportError:
   print >> sys.stderr,"Could not import the subprocess module"
   raise
import multiprocessing
import readline, glob
import subprocess
import psutil
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import random
        
def complete(text, state):
    return (glob.glob(text+'*')+[None])[state]

readline.set_completer_delims(' \t\n;')
readline.parse_and_bind("tab: complete")
readline.set_completer(complete)

class setup :
    @staticmethod
    def get_genome_size(genomefile) :
        chr_list = []
        count_list = []
        count = 0
        for genome in open(genomefile, 'r') :
            if genome[0] == ">" :
                if count > 0 :
                    count_list.append(count)
                chr_list.append(genome[1:].split(" ")[0].strip("\n"))
                count = 0
            else :
                count += len(genome)
        count_list.append(count)
        chr_dict = {}
        for chr in chr_list :
            chr_dict[chr] = count_list[chr_list.index(chr)]
        return chr_dict, chr_list
        
    @staticmethod
    def check_cores(cores):
        if cores > int(multiprocessing.cpu_count()) :
            print "Please reduce the no. of cores"
            sys.exit()

class methylkit:
    @staticmethod
    def prep_methylkit(chr, samname, context, genomefile, randid) :
        import bisect
        ### get rid of cycles and deal with one chr at a time?
        
        print chr, "reading chr genome", datetime.now()
        ### open tmp chrom file and save all the coordinates
        chr_whole_list = []
        for a in open("tmp"+str(randid)+"/"+str(chr)+".out", 'r') :
            chr_whole_list.append(a)
        print chr, "saving cytosine positions", datetime.now()
        cg_loc_list = []
        loc = 0
        if context == "CpG" :
            for base in chr_whole_list[0] : 
                loc += 1
                if base == "C" :
                    if chr_whole_list[0][loc-1:loc+1] == "CG" :
                        cg_loc_list.append(loc)
                elif base == "G" :
                    if chr_whole_list[0][loc-2:loc] == "CG" :
                        cg_loc_list.append(loc)
        elif context == "CHG" : 
            chg_list = ["CAG", "CTG", "CCG"]
            chg_rev_list = ["CTG", "CAG", "CGG"]
            for base in chr_whole_list[0] : 
                loc += 1
                if base == "C" : 
                    if chr_whole_list[0][loc-1:loc+2] in chg_list  :
                        cg_loc_list.append(loc)
                elif base == "G" :
                    if chr_whole_list[0][loc-3:loc] in chg_rev_list :
                        cg_loc_list.append(loc)
        elif context == "CHH" : 
            chh_list = ["CAA", "CAC", "CAT", "CCA", "CCC", "CCT", "CTA", "CTC", "CTT"]
            chh_rev_list = ["TTG", "GTG", "ATG", "TGG", "GGG", "AGG", "TAG", "GAG", "AAG"]
            for base in chr_whole_list[0] : 
                loc += 1
                if base == "C" : 
                    if chr_whole_list[0][loc-1:loc+2] in chh_list  :
                        cg_loc_list.append(loc)
                elif base == "G" :
                    if chr_whole_list[0][loc-3:loc] in chh_rev_list :
                        cg_loc_list.append(loc)
        cg_loc_list_np = np.asarray(cg_loc_list)
        del chr_whole_list
        del cg_loc_list
        cg_loc_list = cg_loc_list_np
        max_cg_val = np.amax(cg_loc_list)+1
        
        print chr, "saving cg loc file", datetime.now()
        cg_loc_file = open("tmp"+str(randid)+"/"+str(chr)+"_cg_loc.txt", 'w') 
        for pos in cg_loc_list :
            cg_loc_file.write(str(pos)+"\n")
        cg_loc_file.close()
        ### check chr or chr+chr
        fp = open(str(context)+"_OT_"+str(samname).split(".bam")[0]+".txt")
        for i, line in enumerate(fp):
            if i == 2:
                if line.split("\t")[2][:3] == "chr" :
                    chrid = "chr"
                else :
                    chrid = ""
                break
        fp.close()
        
        ### write into separate split into the list top and bottom
        print chr, "writing chr separate bme file", datetime.now()
        context_file = open("tmp"+str(randid)+"/"+str(samname).split(".bam")[0]+"-"+str(chr)+"_bme_top.out",'w')
        for a_f in open(str(context)+"_OT_"+str(samname).split(".bam")[0]+".txt", 'r') :
            if len(a_f.split("\t")) > 2 and a_f.split("\t")[2] == chrid+str(chr) :
                context_file.write(a_f)
        if os.path.isfile(str(context)+"_CTOT_"+str(samname).split(".bam")[0]+".txt") :
            for a_f in open(str(context)+"_CTOT_"+str(samname).split(".bam")[0]+".txt", 'r') :
                if len(a_f.split("\t")) > 2 and a_f.split("\t")[2] == chrid+str(chr) :
                    context_file.write(a_f)
        context_file.close()
        context_file = open("tmp"+str(randid)+"/"+str(samname).split(".bam")[0]+"-"+str(chr)+"_bme_bottom.out",'w')
        for a_f in open(str(context)+"_OB_"+str(samname).split(".bam")[0]+".txt", 'r') :
            if len(a_f.split("\t")) > 2  and a_f.split("\t")[2] == chrid+str(chr) :
                context_file.write(a_f)
        if os.path.isfile(str(context)+"_CTOB_"+str(samname).split(".bam")[0]+".txt") :
            for a_f in open(str(context)+"_CTOB_"+str(samname).split(".bam")[0]+".txt", 'r') :
                if len(a_f.split("\t")) > 2 and a_f.split("\t")[2] == chrid+str(chr) :
                    context_file.write(a_f)
        context_file.close()
        
        print chr, "making large list", datetime.now(), max_cg_val
        large_list = np.zeros((max_cg_val,4))
        outfile = open("tmp"+str(randid)+"/"+str(context)+"_"+str(samname)+"_chr"+str(chr)+".methylKit" ,'w')
        
        print chr, "reading cytosine positions", datetime.now()
        for a_f in open("tmp"+str(randid)+"/"+str(samname).split(".bam")[0]+"-"+str(chr)+"_bme_top.out", 'r') : 
            if a_f.split("\t")[1] == "+" :
                large_list[int(float(a_f.split("\t")[3]))][0] += 1
            elif a_f.split("\t")[1] == "-" :
                large_list[int(float(a_f.split("\t")[3]))][1] += 1
        for a_f in open("tmp"+str(randid)+"/"+str(samname).split(".bam")[0]+"-"+str(chr)+"_bme_bottom.out", 'r') : 
            if a_f.split("\t")[1] == "+" :
                large_list[int(float(a_f.split("\t")[3]))][2] += 1
            elif a_f.split("\t")[1] == "-" :
                large_list[int(float(a_f.split("\t")[3]))][3] += 1
        
        print chr, "writing .methylKit", datetime.now()
        for cg_loc in cg_loc_list :
            if int(large_list[cg_loc][0])+int(large_list[cg_loc][1]) != 0 :
                outfile.write(str(chr)+"."+str(cg_loc)+"\t")
                outfile.write(str(chr)+"\t")
                outfile.write(str(cg_loc)+"\t")
                outfile.write("F\t") 
                outfile.write(str(int(large_list[cg_loc][0])+int(large_list[cg_loc][1]))+"\t")
                outfile.write(str(float(large_list[cg_loc][0])/float(int(large_list[cg_loc][0])+int(large_list[cg_loc][1]))*100)+"\t")
                outfile.write(str(float(large_list[cg_loc][1])/float(int(large_list[cg_loc][0])+int(large_list[cg_loc][1]))*100)+"\n")
            elif int(large_list[cg_loc][2])+int(large_list[cg_loc][3]) != 0 :
                outfile.write(str(chr)+"."+str(cg_loc)+"\t")
                outfile.write(str(chr)+"\t")
                outfile.write(str(cg_loc)+"\t")
                outfile.write("R\t") 
                outfile.write(str(int(large_list[cg_loc][2])+int(large_list[cg_loc][3]))+"\t")
                outfile.write(str(float(large_list[cg_loc][2])/float(int(large_list[cg_loc][2])+int(large_list[cg_loc][3]))*100)+"\t")
                outfile.write(str(float(large_list[cg_loc][3])/float(int(large_list[cg_loc][2])+int(large_list[cg_loc][3]))*100)+"\n")
        outfile.close()
        del large_list, cg_loc_list
        
    @staticmethod
    def methylkit(file1, file2, type, chr, specie, randid) :
        ### export LD_LIBRARY_PATH="/opt/apps/R/2.15.1/lib64/R/lib:$LD_LIBRARY_PATH"
        ### python setup.py build --r-home /opt/apps/R/2.15.1/ --r-home-modules /work/02114/wonaya/R_libs/
        #os.environ["PATH"] = "/opt/apps/R/2.15.1/bin/:$PATH"
        #os.environ["LD_LIBRARY_PATH"] = "/opt/apps/R/2.15.1/lib64/R/lib:$LD_LIBRARY_PATH"
        import rpy2.robjects as robjects
        from rpy2.robjects.packages import importr
        os.environ["R_LIBS"] = "/home1/02114/wonaya/R/x86_64-unknown-linux-gnu-library/2.15"
        #os.environ["LD_LIBRARY_PATH"] = "/opt/apps/R/2.15.1/lib64/R/lib:$LD_LIBRARY_PATH"
        print chr, "generating methdiff for", file1, file2, "on",type
        os.chdir("tmp"+str(randid)+"/")
        importr('GenomicRanges')
        importr('bigmemory')
        importr('data.table')
        importr('methylKit')
        basedir = str(os.path.abspath(os.path.dirname(__file__)))+"/../R/base.R"
        robjects.r.assign('basedir',basedir)
        robjects.r('''source(basedir)''')
        diffdir = str(os.path.abspath(os.path.dirname(__file__)))+"/../R/diffMeth.R"
        robjects.r.assign('diffdir',diffdir)
        robjects.r('''source(diffdir)''')
        filelist = [file1, file2]
        file1name = str(type)+"_"+file1+"_chr"+str(chr)+".methylKit"
        file1ID = file1.split("_")[0]
        file2name = str(type)+"_"+file2+"_chr"+str(chr)+".methylKit"
        file2ID = file2.split("_")[0]
        filelist = [file1name, file2name]
        robjects.r.assign('rfilelist',filelist)
        robjects.r.assign('type', type)
        robjects.r.assign('file1ID', file1ID)
        robjects.r.assign('file2ID', file2ID)
        robjects.r.assign('specie', specie)
        robjects.r('''suppressMessages(library(base))''')
        robjects.r('''myobj=read(rfilelist, sample.id=list(file1ID,file2ID), assembly=specie,treatment=c(0,1),context=type, resolution="base")''')
        robjects.r('''meth=unite(myobj, destrand=FALSE)''')
        robjects.r('''rm(valid.methylRawObj)''')
        robjects.r('''rm(rfilelist)''')
        robjects.r('''rm(read)''')
        robjects.r('''rm(myobj)''')
        robjects.r('''gc()''')
        robjects.r('''myDiff=calculateDiffMeth(meth)''')
        robjects.r('''rm(meth)''')
        robjects.r('''gc()''')
        difffile = str(file1ID)+"_"+str(file2ID)+"_"+str(type)+"_diff_chr"+str(chr)+".txt"
        robjects.r.assign('difffile', difffile)
        robjects.r('''write.table(myDiff, file=difffile, sep="\t", row.names = FALSE, col.names= FALSE,  quote = FALSE)''')
        os.chdir("..")
        
    @staticmethod
    def eDMR(file1, file2, type, dmc, cpg, qvalue, diffmeth, randid) :
        import rpy2.robjects as robjects
        from rpy2.robjects.packages import importr
        os.environ["R_LIBS"] = "/home1/02114/wonaya/R/x86_64-unknown-linux-gnu-library/2.15"
        file1ID = file1.split("_")[0]
        file2ID = file2.split("_")[0]
        outfile = "tmp"+str(randid)+"/"+str(file1ID)+"_"+str(file2ID)+"_"+str(type)+"_dmr.txt"
        robjects.r.assign('outfile', outfile)
        eDMRdir = str(os.path.abspath(os.path.dirname(__file__)))+"/../R/eDMR_0.5.1.R"
        robjects.r.assign('eDMRdir',eDMRdir)
        robjects.r('''source(eDMRdir)''')
        robjects.r.assign('type', type)
        robjects.r.assign('file1ID', file1ID)
        robjects.r.assign('file2ID', file2ID)
        difffile = "tmp"+str(randid)+"/"+str(file1ID)+"_"+str(file2ID)+"_"+str(type)+"_diff_sorted.txt"
        robjects.r.assign('difffile', difffile)
        robjects.r('''myDiff=read.table(difffile, header=TRUE)''')
        robjects.r('''mydmr=eDMR(myDiff, mode=1, ACF=FALSE)''')
        robjects.r('''write.table(mydmr, file=outfile, sep="\t",quote = FALSE, row.names=FALSE,col.names=FALSE)''')           

parser = OptionParser()
allgroup = OptionGroup(parser, "Required for all function")
parser.add_option_group(allgroup)
allgroup.add_option("--context", dest="context", help="CpG, CHH, CHG or all")
allgroup.add_option("--genome", dest="genome", help="name and directory of genome fasta file")
allgroup.add_option("--specie", dest="specie", help="Specie code, B73, MM9, HG19")
allgroup.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True, help="don't print status messages to stdout")
mkgroup = OptionGroup(parser, "MethylKit specific options","These are MethylKit options\n")
parser.add_option_group(mkgroup)
mkgroup.add_option("--cores", dest="cores", help="no. of cores to use in running BisKit")
mkgroup.add_option("--bam1", dest="bamfile1", help="name and directory of first bam file generated from bismark for DMR finding")
mkgroup.add_option("--bam2", dest="bamfile2", help="name and directory of second bam file generated from bismark for DMR finding")
mkgroup.add_option("--dmc", dest="dmc", help="Theshold DMC count to filter DMRs, default = 1", default=1)
mkgroup.add_option("--cpg", dest="cpg", help="No. of CpG in a DMR, default = 3", default=3)
mkgroup.add_option("--qvalue", dest="qvalue", help="Q-value of DMRs to print", default=0.05)
mkgroup.add_option("--diffmeth", dest="diffmeth", help="Only print out DMRs with methylation difference greater than", default=20)
(options, args) = parser.parse_args()

randid = str(random.randint(10000, 99999))
if not os.path.isdir("tmp"+randid) :
    os.mkdir("tmp"+randid)
else :
    randid = str(random.randint(10000, 99999))
    os.mkdir("tmp"+randid)
    
print "BisuKit starting up", datetime.now()
setup.check_cores(int(options.cores))
cores = int(options.cores)
chr_list = setup.get_genome_size(options.genome)[1]
genome_list = []
whole_list = []
chr_list = []
chr_order_dict = {}
start = 0
for genome in open(options.genome, 'r') :
    if genome[0] == ">" :
        if len(chr_list) == 0 :
            chr_list.append(genome.split(" dna")[0][1:])
            chr_order_dict[genome.split(" dna")[0][1:]] = start
        else :
            whole_list.append("".join(genome_list))
            chr_list.append(genome.split(" dna")[0][1:])
            start += 1
            chr_order_dict[genome.split(" dna")[0][1:]] = start
            genome_list = []
    else :
        genome_list.append(genome.strip("\n"))

whole_list.append("".join(genome_list))
del genome_list
total_rounds = (len(chr_list)/int(options.cores))+1

### write chr sequence into file
print "writing genome to tmp"
for chrom in chr_list[-5:-3] :
    outfile = open("tmp"+str(randid)+"/"+str(chrom)+".out", 'w')
    for a in whole_list[int(chr_order_dict[chrom])] :
        outfile.write(a)
    outfile.close()
del chr_order_dict
del whole_list

print "prep bam1 and bam2", datetime.now()
for chr in chr_list[-5:-3] :
    jobs = []
    s1 = multiprocessing.Process(target=methylkit.prep_methylkit, args=(chr, options.bamfile1, options.context, options.genome, randid, ))
    s2 = multiprocessing.Process(target=methylkit.prep_methylkit, args=(chr, options.bamfile2, options.context, options.genome, randid, ))
    jobs.append(s1)
    jobs.append(s2)
    s1.start()
    s2.start()
    [x.join() for x in jobs]

print "methylkit start", datetime.now()
jobs = []
for round in range(0, total_rounds) :
    jobs = []
    #for chr in chr_list[round*cores:(round+1)*cores]:
    for chr in chr_list[-5:-3] :
        s1 = multiprocessing.Process(target=methylkit.methylkit, args=(options.bamfile1, options.bamfile2, options.context, chr, options.specie, randid, ))
        jobs.append(s1)
        s1.start()
    [x.join() for x in jobs]

### merge meth  
print "merging methylKit output", datetime.now()
whole_meth_file = open("tmp"+str(randid)+"/"+str((options.bamfile1).split("_")[0])+"_"+str((options.bamfile2).split("_")[0])+"_"+str(options.context)+"_diff.txt", 'w') 
#for chr in chr_list :
for chr in  chr_list[-5:-3] :
    a = open("tmp"+str(randid)+"/"+str((options.bamfile1).split("_")[0])+"_"+str((options.bamfile2).split("_")[0])+"_"+str(options.context)+"_diff_chr"+str(chr)+".txt", 'r')
    ### remove duplicates and 0 meth
    alines = a.readlines()
    for line in set(alines) :
        if float(line.split("\t")[-1].strip("\n")) != 0 :
            whole_meth_file.write(line)
whole_meth_file.close()

os.system("sort -nk1 -nk2 tmp"+str(randid)+"/"+str((options.bamfile1).split("_")[0])+"_"+str((options.bamfile2).split("_")[0])+"_"+str(options.context)+"_diff.txt > tmp"+str(randid)+"/"+str((options.bamfile1).split("_")[0])+"_"+str((options.bamfile2).split("_")[0])+"_"+str(options.context)+"_diff.sorted.txt")
whole_meth_file = open("tmp"+str(randid)+"/"+str((options.bamfile1).split("_")[0])+"_"+str((options.bamfile2).split("_")[0])+"_"+str(options.context)+"_diff_sorted.txt", 'w') 
whole_meth_file.write("chr\tstart\tend\tstrand\tpvalue\tqvalue\tmeth.diff\n")
for a in open("tmp"+str(randid)+"/"+str((options.bamfile1).split("_")[0])+"_"+str((options.bamfile2).split("_")[0])+"_"+str(options.context)+"_diff.sorted.txt",'r'):
    whole_meth_file.write(a)
whole_meth_file.close()

print "running eDMR", datetime.now()
methylkit.eDMR(options.bamfile1, options.bamfile2, options.context, options.dmc, options.cpg, options.diffmeth, options.qvalue, randid)
bedgraph = open((options.bamfile1).split("_")[0]+"_"+(options.bamfile2).split("_")[0]+"_"+options.context+"_DMR.bedGraph", 'w')
outfile = open((options.bamfile1).split("_")[0]+"_"+(options.bamfile2).split("_")[0]+"_"+options.context+"_DMR.txt", 'w')
outfile.write("chr\tstart\tend\twidth\tstrand\tmean.meth.diff\tnum.CpGs\tnum.DMCs\tDMR.pvalue\tDMR.qvalue\n")
afile = open("tmp"+str(randid)+"/"+(options.bamfile1).split("_")[0]+"_"+(options.bamfile2).split("_")[0]+"_"+options.context+"_dmr.txt", 'r')
alines = afile.readlines()
for a in alines[1:] :
    if float(a.split("\t")[5]) >= float(options.diffmeth) or float(a.split("\t")[5]) <= float((options.diffmeth))*-1 :
        if float(a.split("\t")[7]) >= float(options.dmc) and float(a.split("\t")[6]) >= float(options.cpg) and float(a.split("\t")[-1]) <= float(options.qvalue) :
            outfile.write(a)
            bedgraph.write("\t".join(a.split("\t")[0:3])+"\t"+a.split("\t")[5]+"\n")
outfile.close()
bedgraph.close()
os.rmdir("tmp"+str(randid)+"/")
print "BisuKit complete", datetime.now()
