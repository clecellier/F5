#!/usr/bin/python
#-*- coding: utf-8 -*-
import os
import re
from operator import itemgetter
from matplotlib import pyplot as plt
import argparse


parser = argparse.ArgumentParser(description="This script generates a background FASTA file from a foreground FASTA file")
parser.add_argument('-f', '--foreground', help='foreground file name', required=True)

parser.add_argument('-b', '--background', help='background file name', required=True)

parser.add_argument('-o', '--output', help='output file name', required=True)

args = vars(parser.parse_args())
########### Create one dictionnary of Sequences and one dictionnary of GC ratio##########################

def readFASTAfile(fd,dico):# or Biopython solution
    regex1= re.compile ('^>(.+)\n$')#Watch out novelmRNA!!!
    regex2= re.compile ('^([ACGT]+)\n$')
    m_file = fd.readlines()
    for line in m_file:
        res1 = regex1.search(line) 
        res2 = regex2.search(line)
        if res1:
            m_miR = res1.group(1)
        elif res2:
            m_seq = res2.group(1)
            if dico.has_key(m_miR):           
                dico[m_miR] += m_seq
            else:
                dico[m_miR] = m_seq

def GCratio(dico,dicoGC):
    for m_seq in dico:
        nbG = dico[m_seq].count('G')
        nbC = dico[m_seq].count('C')
        dicoGC[m_seq] = float(nbG+nbC)/len(dico[m_seq])
      
def FindClosestGC(val,m_list):
    index = len(m_list)-1# lists are 0-based
    while val > m_list[index][1] and index > 0:
        index -=1
    elem = m_list[index][0]
    del m_list[index]
    return elem #key of dicoBGGC
       
######### foreground ##########
dicoFGSeq ={}
dicoFGGC = {}
with open(args['foreground'], "r") as fd:
    readFASTAfile(fd,dicoFGSeq)
    GCratio(dicoFGSeq,dicoFGGC)

######## background ############

dicoBGSeq ={}
dicoBGGC = {}
with open(args['background'], "r") as fd:
    readFASTAfile(fd,dicoBGSeq)
    GCratio(dicoBGSeq,dicoBGGC)
    
m_list = dicoBGGC.items()
m_list.sort(key=itemgetter(1),reverse=True) # list of tuples (sample,GCratio) sorted on decreasing GC ratio

with open(args['output'], "w") as res: 
        for i in dicoFGGC:   
            m_key = FindClosestGC(dicoFGGC[i],m_list)
            if m_key in dicoBGSeq:
                    res.write(">")
                    res.write(m_key)
                    res.write("\n")
                    res.write(dicoBGSeq[m_key])
                    res.write("\n")
#print len(dicoFGGC)
#print len(dicoBGGC)

"""
ListGCfg = []
ListGCbg = []
for i in dicoFGGC:
    m_key = FindClosestGC(dicoFGGC[i],m_list)
    ListGCfg.append(dicoFGGC[i]) 
    ListGCfg.sort(reverse=True)
    ListGCbg.append(dicoBGGC[m_key])
    ListGCbg.sort(reverse=True)


plt.figure(1) 
#x=np.linspace(-5,5,100)
plt.subplot(211)
plt.plot(ListGCfg, 'bs')
#plt.title('Foreground')
plt.text(200, .7, 'Foreground',fontsize=14)
plt.ylabel('percentage GC')
plt.grid(True)
plt.subplot(212)
plt.plot(ListGCbg, 'g^')
#plt.title('Background')
plt.text(200, .7, 'Background',fontsize=14)
plt.ylabel('percentage GC')
plt.xlabel("sequences")
plt.grid(True)
plt.show()
"""
