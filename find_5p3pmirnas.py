#!/usr/bin/python
#-*- coding: utf-8 -*-
import sys, os, re


def read_bedfile(stream):
    stream.readline()
    for line in stream:
        res = re.search("^(chr.+)\t([0-9]+)\t([0-9]+)\t(hsa-.+)\t0\t([+-])$", line)
        if res:
            start = res.group(2)
            end = res.group(3)
            chrom = res.group(1)
            mir = res.group(4)
            strand = res.group(5)
            if strand =="+":
                print chrom,"\t",int(start)-1501,"\t",int(start)+20,"\t", mir,"\t0\t",strand 
            elif strand =="-":
                print chrom,"\t",int(end)-20,"\t",int(end)+1501,"\t", mir,"\t0\t",strand 


def read_bedfile_2(stream):
    stream.readline()
    for line in stream:
        res = re.search("^MI[0-9]+\t(chr.+)\t([0-9]+)\t([0-9]+)\t(hsa-.+)\t0\t([+-])\t([0-9]+)\t([0-9]+)\thsa-.+$", line)
        if res:
            start_premir = res.group(2)
            end_premir = res.group(3)
            chrom = res.group(1)
            mir = res.group(4)
            strand = res.group(5)
            start_mir = res.group(6)
            end_mir = res.group(7)
            if strand == "+" and int(start_mir) - int(start_premir) >= 40:
                print chrom,"\t",int(start_premir)-1501,"\t",int(start_premir)+20,"\t", mir,"\t0\t",strand 
            elif strand == "-" and int(end_premir) - int(end_mir) >= 40:
                print chrom,"\t",int(end_premir)-20,"\t",int(end_premir)+1501,"\t", mir,"\t0\t",strand 

################################################################################

"""            
#### homo sapiens
            
with open("mir_producing5pAND3p.bed","r") as stream:
    read_bedfile(stream)

#with open("intersect_pre-mirANDdpi.bed","r") as stream:
#    read_bedfile(stream)



dico_start_mir = {}
with open("miRNAs_single_arm.bed","r") as stream:
    dico_start_mir = read_bedfile(stream, dico_start_mir) # no distinction of miR-A-1 and miR-A-2

#print len(dico_start_mir)

dico_start_premir = {}
with open("pre-mirnas_start.bed","r") as stream2:
    dico_start_premir = read_bedfile(stream2, dico_start_premir)

print len(dico_start_premir)
#for mir in dico_start_premir:
#    if mir in dico_start-mir and dico_start_mir[mir]-dico_start_premir[mir]>=40:
#        print mir

"""
list_mir = []
with open("intersect_pre-mirANDdpi/intersect_pre-mirANDdpi_list.txt","r") as stream:
    for line in stream:
#res = re.search("^(hsa-.+)\n$", line)
#if res:
#mir = res.group(1)
#list_mir.append(mir)
        list_mir.append(line.strip("\n"))
    
with open("pre-miRNAs.bed","r") as stream:
    stream.readline()
    for line in stream:
        res = re.search("^(chr.+)\t([0-9]+)\t([0-9]+)\t(hsa-.+)\t0\t([+-])$", line)
        if res:
            start = res.group(2)
            end = res.group(3)
            chrom = res.group(1)
            mir = res.group(4)
            strand = res.group(5)
            if mir in list_mir:
                if strand=="+":
                    print chrom+"\t"+str(int(start)-10001)+"\t"+str(end)+"\t"+mir+"\t0\t"+strand
                    ## keep the pre-mir
                else:
                    print chrom+"\t"+str(start)+"\t"+str(int(end)+10001)+"\t"+mir+"\t0\t"+strand



"""

##### mus musculus

list_mir = []
with open("Xie_et_al_Cell2013/list_Xie_et_al_mmu.txt","r") as stream:
    for line in stream:
        res = re.search("^(mmu-.+)\n$", line)
        if res:
            mir = res.group(1)
            list_mir.append(mir)

    
with open("Xie_et_al_Cell2013/premir_mmu_sed.bed","r") as stream:
    stream.readline()
    for line in stream:
        res = re.search("^(chr.+)\t([0-9]+)\t([0-9]+)\t(mmu-.+)\t0\t([+-])$", line)
        if res:
            start = res.group(2)
            end = res.group(3)
            chrom = res.group(1)
            mir = res.group(4)
            strand = res.group(5)
            if mir in list_mir:
                if strand=="+":
                    print chrom+"\t"+str(int(start)-50)+"\t"+str(int(start)+10)+"\t"+mir+"\t0\t"+strand
                else:
                    print chrom+"\t"+str(int(end)-10)+"\t"+str(int(end)+50)+"\t"+mir+"\t0\t"+strand
dico = {}
with open("all_miRNAs_in_miRBase.txt","r") as stream2:
    for line in stream2:
        for mir in list_mir:
            if mir+"-5p" in line or mir+"-3p" in line:
                dico[mir] = "5p3p arms"
            elif mir in line:
                dico[mir] = "one arm"

for mir in dico:
    print mir+"\t"+dico[mir]

"""
