#!/usr/bin/python
#-*- coding: utf-8 -*-
import sys
import os
import re
import scipy
import numpy as np

from operator import itemgetter

from scipy.stats import pearsonr #values
from scipy.stats import spearmanr #ranks
from scipy.stats import kendalltau #ranks
from matplotlib import pyplot as plt


#######################################################################
"""
This script takes a .tsv file as input containing the TPM of
different miRNAs in different samples as well as DPI_CAGE in columns. The first
column must correspond to the samples:

samples\tmir-1\tmir-2\tTSS_id-1\tTSS-id-2
"""
#######################################################################

def parse(fd, m_dict):
    headers = fd.readline().split(" ") #check the separator
    for h in headers[1:]: # samples is no key
       m_dict[h] = []
    for line in fd:
        TPMs = line.split(" ")
        for h in headers[1:]:#excludes samples name from the list
            m_dict[h].append(TPMs[headers.index(h)])
    return m_dict 

def compare(e1, e2):
    r1 = re.search("(hsa.+):hsa.+", e1)
    r2 = re.search("(hsa.+):chr.+", e2)
    if r1 and r2 and r1.group(1) == r2.group(1):
        return True
    else:
        return False
    
#######################################################################

    
m_file = "intersect_pre-mirANDdpi/intersect_pre-mirANDdpi_TSSsin_10kbupstreampremir_all_dpi_all_data.txt"
m_dict = {} #dict[header]=column

with open(m_file,"r") as fd:
    m_dict = parse(fd, m_dict)

#x = np.array(m_dict["hsa-mir-23a:chr19:13947489..13947497,-"], dtype = float)
#print str(np.amax(x))


print "mirna\tTSS_id\tcorrelation\tp_val"
for e1 in m_dict:
    for e2 in m_dict:
        if compare(e1, e2):
#print m_dict[e1]
#print m_dict[e2]            
            tpm_mir = np.array(m_dict[e1], dtype = float) 
            tpm_dpi = np.array(m_dict[e2], dtype = float)
            if np.amax(tpm_dpi) >= 0.9: # only consider significant TSSs - mir-320a and mir-23a/24-2 as benchmarks
#r_row, p_value = pearsonr(tpm_mir, tpm_dpi)
                r_row, p_value = spearmanr(tpm_mir, tpm_dpi)
#r_row, p_value = kendalltau(tpm_mir, tpm_dpi)
                print e1+"\t"+e2+"\t"+str(r_row)+"\t"+str(p_value)

"""
######################################################################                
## focus on hsa-miR-320a
## plot TPM
plt.figure(1)
t1 = m_dict["hsa.mir.320a"]
x = range(len(t1))
t2 = m_dict["hsa.mir.320a.chr8.22102875..22102886.."]
t3 = m_dict["hsa.mir.320a.chr8.22102848..22102855.."]
t4 = m_dict["hsa.mir.320a.chr8.22102538..22102543.."] #TSS of the 5'-capped pre-miRNA


t2 = np.array(m_dict["hsa.mir.320a.chr8.22102875..22102886.."], dtype = float)
t3 = np.array(m_dict["hsa.mir.320a.chr8.22102848..22102855.."], dtype = float)
t4 = np.array(m_dict["hsa.mir.320a.chr8.22102538..22102543.."], dtype = float) #TSS of the 5'-capped pre-miRNA

print "\thsa.mir.320a.chr8.22102875..22102886..\thsa.mir.320a.chr8.22102848..22102855..\thsa.mir.320a.chr8.22102538..22102543.."
print "median\t",str(np.median(t2)),"\t",str(np.median(t3)),"\t", str(np.median(t4))
print "mean\t", str(np.mean(t2)),"\t",str(np.mean(t3)),"\t", str(np.mean(t4))


#plt.plot(x,t1,"k",x,t2,"b",x,t3,"r",x,t4,"g")
plt.subplot(121)
plt.plot(x,t1,'y')
sub1 = plt.subplot(121)
sub1.set_title('TPM_expression')

plt.subplot(122)
plt.plot(x,t2,'b')
plt.plot(x,t3,'r')
plt.plot(x,t4,'g')
sub2 = plt.subplot(122)
sub2.set_title('TPM_CAGE')
plt.legend(("chr8.22102875..22102886", "chr8.22102848..22102855", "chr8.22102538..22102543"), 'best')
plt.ylabel("TPM")
plt.xlabel("samples")
plt.show()
"""

