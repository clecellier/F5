#!/usr/bin/python
#-*- coding: utf-8 -*- 
import sys
import  os
import re
import numpy as np

def parse(stream):
    print stream.readline(), # write hearders
    for line in stream: 
        res = re.search("^(hsa-.+) (hsa.+) hsa.+:chr[0-9X]+:([0-9]+)\.\.([0-9]+),[-+] ([-0-9\.e]+) ([0-9\.e-]+) ([0-9]+) ([0-9]+) (chr[0-9X]+) ([-+])$", line)
        if res:
            pre = res.group(1)
            mir = res.group(2)
            tss_start = res.group(3)
            tss_end = res.group(4)
            rho = res.group(5)
            pval = res.group(6)
            pre_start = res.group(7)
            pre_end = res.group(8)
            chrom = res.group(9)
            strand = res.group(10)
            """
            if strand == "+" and int(tss_start) <= int(pre_start):
                print line,
            elif strand == "+" and int(tss_start) > int(pre_start) and (int(tss_start) - int(pre_start)) <= 15:
                print line,
            elif strand == "-" and int(tss_end) >= int(pre_end):
                print line,
            elif strand == "-" and int(tss_end) < int(pre_end) and (int(pre_end) - int(tss_end)) <= 15:
                print line,

            """
            if strand=="+" and abs(int(tss_end) - int(pre_start)) <= 20:# and float(rho) >= 0.12 :
                print line,
            elif strand=="-" and abs(int(tss_start) - int(pre_end)) <= 20:# and float(rho) >= 0.12:
                print line,
            
            """
            if float(rho) >= 0.12 # 0.12 for hsa-mir-23a or 0.17 for hsa-mir-320a :
                print line,
            """

def compute(stream):
    m_list = parse2(stream)
    x = np.array(m_list, dtype = int)
    median = np.median(x)
    mean = np.mean(x) 
    xmin = np.amin(x)
    xmax = np.amax(x)
    print "min: "+str(xmin)," max: ", str(xmax)," median: ", str(median)," mean: ", str(mean)

def parse2(stream):
    list_size_tss = []
    print stream.readline(), # write hearders
    for line in stream: 
        res = re.search("^(hsa-.+) (hsa.+) hsa.+:chr[0-9X]+:([0-9]+)\.\.([0-9]+),[-+] ([-0-9\.e]+) ([0-9\.e-]+) ([0-9]+) ([0-9]+) (chr[0-9X]+) ([-+])$", line)
        if res:
            pre = res.group(1)
            mir = res.group(2)
            tss_start = res.group(3)
            tss_end = res.group(4)
            rho = res.group(5)
            pval = res.group(6)
            pre_start = res.group(7)
            pre_end = res.group(8)
            chrom = res.group(9)
            strand = res.group(10)
            
            diff = (int(tss_end) - int(tss_start))
            list_size_tss.append(diff)
    
    return list_size_tss        

def parse3(stream): # for .bed file
    list_size_tss = []
    for line in stream:
        res = re.search("^chr[0-9XYM]+\t([0-9]+)\t([0-9]+)\t.+$", line)
        if res:
            start = res.group(1)
            end = res.group(2)
            diff = (int(end) - int(start))
            list_size_tss.append(diff)
    return list_size_tss


#################################################################################
m_file = "intersect_pre-mirANDdpi/correlation_TPM_sup09/correlation_spearman_TPMsup09_coord_noNAN.txt"

with open(m_file,"r") as stream:
    parse(stream)
#compute(stream)
