#!/usr/bin/python
#-*- coding: utf-8 -*-

import sys
import  os
import re

def find_coord_premir(fd, dico_coord_premirna):
    for line in fd:
        res = re.search("^(hsa-.+)\t([0-9]+)\t([0-9]+)$", line)
        if res:
            premir = res.group(1)
            start = res.group(2)
            end = res.group(3)
#dico_coord_premirna[premir] = (start, end)
            dico_coord_premirna[premir] = start
            
    return dico_coord_premirna

def find_start_premir(fd, dico_coord_premirna):
    for line in fd:
        res = re.search("^(hsa-.+)\t([0-9]+)\t[0-9]+\t(chr[XY0-9]{1,2})$", line)
        if res:
            premir = res.group(1)
            start = res.group(2)
            chrom = res.group(3)

            print chrom+"\t"+str(int(start)-10)+"\t"+str(int(start)+11)+"\t"+premir+"\t0\t+"

            dico_coord_premirna[premir] = chrom+":"+start
            
    return dico_coord_premirna
        
def find_coord_dpi(fd, dico_coord_dpi):
    for i, line in enumerate(fd):
        res = re.search("^chr[XYM0-9]{1,2}:([0-9]+)..([0-9]+),[+-]$", line)
        if res:
            start = res.group(1)
            end = res.group(2)
#dico_coord_dpi[i] = (start, end)
            dico_coord_dpi[i] = start
    return dico_coord_dpi

def find_start_dpi(fd, dico_coord_dpi):
    for i, line in enumerate(fd):
        res = re.search("^(chr[XYM0-9]{1,2}:[0-9]+)..([0-9]+),[+-]$", line)
        if res:
            start = res.group(1)
            dico_coord_dpi[i] = start
    return dico_coord_dpi


def compare_dict(dict1, dict2):
    for k1 in dict1:
        if dict1[k1] in dict2.values():
                return k1

def readfiles(filenames):# read files of a directory
    for f in filenames:
        for line in open(f):
            yield line


def py_grep(regex, lines):
    print regex
    return [line for line in lines if regex in line]

#res = re.compile(regex)
#for line in stream:
#if res.search(line):
#yield line

def printlines(lines):
    for line in lines:
        print line,

#################################################################################
dico_coord_premirna = {}
dico_coord_dpi = {}

dict_tss = {}
with open("intersect_pre-mirANDdpi/intersect_pre-mirANDdpi_TSSsin_10kbupstreampremir.bed","r") as fd: 
    for line in fd:
#res = re.search("^(chr.+\t[0-9]+\t[0-9]+\thsa-.+\t[-+])$", line)
        res = re.search("^(chr.+)\t([0-9]+)\t([0-9]+)\t(hsa-.+)\t0\t([-+])$", line)        
        if res:        
#dict_tss[res.group(1)+":"+res.group(2)+".."+res.group(3)+","+res.group(5)] = res.group(4)
#dict_tss[res.group(1)] = 1 
            if res.group(4) in dict_tss:
                dict_tss[res.group(4)].append(res.group(1)+":"+res.group(2)+".."+res.group(3)+","+res.group(5))
            else:
                dict_tss[res.group(4)] = []
                dict_tss[res.group(4)].append(res.group(1)+":"+res.group(2)+".."+res.group(3)+","+res.group(5))
            
#print len(dict_tss) ## OK
#print dict_tss["hsa-mir-320a"] ## OK
#m = 0
#for i in dict_tss:
#    m += len(dict_tss[i])
#print m

dict_tpm = {}
#with open("hg19_robust_dpi_median_sup0.5.txt", "r") as stream:
with open("hg19.tc.decompose_smoothing_merged.ctssMaxCounts3.clustername_update__description_expression_tpm_rle.matrix", "r") as stream:
    for line in stream:
        res = re.search("^chr.+\t[0-9]+\t[0-9]+\t(chr.+,[-+])\t.+$", line)        
        if res:        
            dict_tpm[res.group(1)] = line
#print len(dict_tpm) ## OK
#print dict_tpm["chr8:22102538..22102543,-"] ## OK
#print dict_tpm["chr8:22103656..22103662,-"] ## key error but "chr8:22103656..22103662,+" OK


for m in dict_tss:
    for e in dict_tss[m]:
        if e in dict_tpm.keys():
            print m+"\t"+dict_tpm[e],

############################# DRAFT ##################
"""
#with open("hg19.tc.decompose_smoothing_merged.ctssMaxCounts3.clustername_update__description_expression_tpm_rle.matrix","r") as stream:
with open("hg19_robust_dpi_median_sup0.5.txt", "r") as stream:
    for m in dict_tss:
        for e in dict_tss[m]:
            lines = py_grep(e, stream)
            printlines(lines)
#if e == "chr8:22102538..22102543,-":
#for w in py_grep(e, stream):
#print w.next()
    
for t in dict_tss:
    if os.popen("grep "+t+" hg19_robust_dpi_sed.txt","r",1).read() != "":
        print dict_tss[t]+"\t"+os.popen("grep "+t+" hg19_robust_dpi_sed.txt","r",1).read()
            
with open("hg19_robust_dpi_mean_sup0.5.txt", "r") as stream:#mean>=0.5 excludes TSSs with low TPM for hsa-mir-320a
    for line in stream:
        for t in dict_tss:
            if t in line:
                print dict_tss[t]+"\t"+line



with open("/Users/admin/Documents/projects_papers/FANTOM5_miRNA_TSS/human.miRNA.TPM.RLE.csv","r") as stream:
    print stream.readline()
    for line in stream:
        for mir in list_mir:
            if mir in line:
                print line
    
#for mir in list_mir:
#   l = os.popen("grep "+mir+" /Users/admin/Documents/projects_papers/FANTOM5_miRNA_TSS/human.miRNA.TPM.RLE.csv","r")
#   print l

#with open("genomic_coord_pre-miRNAs3.txt","r") as fd:
#    dico_coord_premirna = find_start_premir(fd, dico_coord_premirna)

#with open("genomic_coord_dpi.txt", "r") as fd2:
#    dico_coord_dpi = find_start_dpi(fd2, dico_coord_dpi)

#print len(dico_coord_premirna)

#print compare_dict(dico_coord_premirna, dico_coord_dpi)


"""
