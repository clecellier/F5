#!/usr/bin/python
#*-* coding: utf-8 *-*

""" Module constructing miRNA clusters from predicted miRNA TSSs. """

import sys
import getopt


def construct_mir_tss_clusters(asso_file, corr_thr):
    with open(asso_file) as stream:
        clusters = {}
        stream.readline()  # Remove header
        for line in stream:
            spl = line.split(',')
            mir = spl[0]
            tss = "\t".join(spl[1:5])
            corr = eval(spl[5])
            if corr > corr_thr:
                if tss in clusters:
                    clusters[tss].append((mir, corr))
                else:
                    clusters[tss] = [(mir, corr)]
        return clusters


def print_closeby_clusters(clust_distrib, clusters):  # TODO decomplexify
    for chromo in clust_distrib:
        if len(clust_distrib[chromo]) > 1:
            for index, (start, end, strand) in enumerate(
                    clust_distrib[chromo]):
                for indice in xrange(index + 1, len(clust_distrib[chromo])):
                    c_start, c_end, c_strand = clust_distrib[chromo][indice]
                    max_dist = 10000
                    if (abs(c_start - start) < max_dist or
                            abs(c_start - end) < max_dist or
                            abs(c_end - start) < max_dist or
                            abs(c_end - end) < max_dist):
                        tss1 = "{0}\t{1:d}\t{2:d}\t{3}".format(chromo, start,
                                                               end, strand)
                        tss2 = "{0}\t{1:d}\t{2:d}\t{3}".format(chromo, c_start,
                                                               c_end, c_strand)
                        print "\n\nCloseby clusters for TSSs:"
                        print "{0}\t{1:d}\t{2:d}\t{3}".format(chromo, start,
                                                              end, strand)
                        print "{0}\t{1:d}\t{2:d}\t{3}".format(chromo, c_start,
                                                              c_end, c_strand)
                        print "Cluster from TSS {0}:".format(tss1)
                        for mir, corr in clusters[tss1]:
                            print "{0} with correlation {1:f}".format(
                                mir, corr)
                        print "Cluster from TSS {0}:".format(tss2)
                        for mir, corr in clusters[tss2]:
                            print "{0} with correlation {1:f}".format(
                                mir, corr)


def print_clusters(clusters):
    count = 1
    clust_distrib = {}
    for tss in clusters:
        if len(clusters[tss]) > 2:
            print "\n\n########################################"
            print "Cluster n°{0:d} from TSS {1}:".format(count, tss)
            for mir, corr in clusters[tss]:
                print "{0} with correlation {1:f}".format(mir, corr)
            print "#########################################"
            count += 1
            spl = tss.split('\t')
            chromo = spl[0]
            start = eval(spl[1])
            end = eval(spl[2])
            strand = spl[3]
            if chromo in clust_distrib:
                clust_distrib[chromo].append((start, end, strand))
            else:
                clust_distrib[chromo] = [(start, end, strand)]
    print_closeby_clusters(clust_distrib, clusters)


def main():
    usage = '''
    %s -f <mirna-tss association file> -c <correlation threshold>
    ''' % (sys.argv[0])

    try:
        opts, _ = getopt.getopt(sys.argv[1:], "f:c:h")
    except getopt.GetoptError:
        sys.exit(str(getopt.GetoptError) + usage)

    asso_file = None
    corr_thr = None
    for o, a in opts:
        if o == "-f":
            asso_file = a
        elif o == "-c":
            corr_thr = eval(a)
        else:
            sys.exit(usage)
    if not(asso_file and corr_thr):
        sys.exit(usage)

    clusters = construct_mir_tss_clusters(asso_file, corr_thr)
    print_clusters(clusters)


###############################################################################
#                               MAIN
###############################################################################
if __name__ == "__main__":
    main()
