#!/usr/bin/python
#-*- coding: utf-8 -*-

""" Module parsing the master correlation file b/w DPI clusters and miRNAs. """


import sys
import getopt
import re
import numpy as np


def parse(stream):
    print stream.readline(),  # write hearders
    regexp = re.compile(r"(hsa-.+) (hsa.+) hsa.+:chr[0-9X]+:([0-9]+)\.\.([0-9]+),[-+] ([-0-9\.e]+) ([0-9\.e-]+) ([0-9]+) ([0-9]+) (chr[0-9X]+) ([-+])")
    for line in stream:
        res = re.match(regexp, line)
        if res:
            tss_start = res.group(3)
            tss_end = res.group(4)
            pre_start = res.group(7)
            pre_end = res.group(8)
            strand = res.group(10)
            #
            #if strand == "+" and int(tss_start) <= int(pre_start):
            #    print line,
            #elif (strand == "+" and int(tss_start) > int(pre_start) and
            #        (int(tss_start) - int(pre_start)) <= 15):
            #    print line,
            #elif strand == "-" and int(tss_end) >= int(pre_end):
            #    print line,
            #elif (strand == "-" and int(tss_end) < int(pre_end) and
            #        (int(pre_end) - int(tss_end)) <= 15):
            #    print line,
            if strand == "+" and abs(int(tss_end) - int(pre_start)) <= 20:
                # and float(rho) >= 0.12 :
                print line,
            elif strand == "-" and abs(int(tss_start) - int(pre_end)) <= 20:
                # and float(rho) >= 0.12:
                print line,

            #if float(rho) >= 0.12 # 0.12 for hsa-mir-23a or 0.17 for
            #hsa-mir-320a :
            #    print line,


def compute(stream):
    m_list = parse2(stream)
    x = np.array(m_list, dtype=int)
    median = np.median(x)
    mean = np.mean(x)
    xmin = np.amin(x)
    xmax = np.amax(x)
    print "min: {0:f} max: {1:f} median: {2:f} mean: {3:f}".format(xmin, xmax,
                                                                   median,
                                                                   mean)


def parse2(stream):
    list_size_tss = []
    print stream.readline(),  # write hearders
    regexp = re.compile(r"(hsa-.+) (hsa.+) hsa.+:chr[0-9X]+:([0-9]+)\.\.([0-9]+),[-+] ([-0-9\.e]+) ([0-9\.e-]+) ([0-9]+) ([0-9]+) (chr[0-9X]+) ([-+])")
    for line in stream:
        res = re.match(regexp, line)
        if res:
            tss_start = res.group(3)
            tss_end = res.group(4)
            diff = (int(tss_end) - int(tss_start))
            list_size_tss.append(diff)
    return list_size_tss


def parse3(stream):  # for .bed file
    list_size_tss = []
    regexp = re.compile("chr[0-9XYM]+\t([0-9]+)\t([0-9]+)\t.+")
    for line in stream:
        res = re.matci(regexp, line)
        if res:
            start = res.group(1)
            end = res.group(2)
            diff = (int(end) - int(start))
            list_size_tss.append(diff)
    return list_size_tss


def main():
    usage = '''
    %s -f <input file>
    ''' % (sys.argv[0])

    try:
        opts, _ = getopt.getopt(sys.argv[1:], "f:h")
    except getopt.GetoptError:
        sys.exit(str(getopt.GetoptError) + usage)

    m_file = None
    for o, a in opts:
        if o == "-f":
            m_file = a
        else:
            sys.exit(usage)
    if not m_file:
        sys.exit(usage)
    with open(m_file) as stream:
        parse(stream)


############################################################################
#                                   MAIN                                   #
############################################################################
if __name__ == "__main__":
    main()
