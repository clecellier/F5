#!/usr/bin/python
#-*- coding: utf-8 -*-
import sys, os, re
import scipy
import numpy as np
 
from operator import itemgetter 
from scipy.stats import pearsonr #values
from scipy.stats import spearmanr #ranks
from scipy.stats import kendalltau #ranks
from matplotlib import pyplot as plt
#######################################################################

def parse_max(stream):
    for line in stream:
#res=re.search("^chr[XYM0-9]{1,2}\t[0-9]+\t[0-9]+\tchr.+,[+-]\t.\t[+-]\t.+$ ", line)
        #if res:
        tpms = map(float, line.split("\t")[6:])
        if max(tpms) >= 1:
            print line

def parse_median(stream):
    for line in stream:
        t = np.array(line.split("\t")[6:], dtype = float)
        if np.median(t) >= 0.5:
            print line

def parse_mean(stream):
    for line in stream:
        t = np.array(line.split("\t")[6:], dtype = float)
        if np.mean(t) >= 0.4:
            print line


#######################################################################
with open("hg19.tc.decompose_smoothing_merged.ctssMaxCounts3.clustername_update__description_expression_tpm_rle.matrix","r") as stream:
#parse_max(stream) # include miR-320a TSSs with rho >1.7 but weak mean and median
    parse_median(stream)# >=1 too stringent
#parse_mean(stream)

