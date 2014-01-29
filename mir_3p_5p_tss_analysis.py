#!/usr/bin/python
#*-* coding: utf-8 *-*

"""
Multiple analysis of 5p/3p expression wrt predicted TSSs.

Analysis of potentially different promoters for 5p and 3p miRs from the miR-TSS
associations.

"""

import re
import argparse


# TODO: reduce variable number
def get_top_mir_tss_asso(mir_tss_file, correlation_thr, pval_thr, pls_thr):
    with open(mir_tss_file) as stream:
        asso = {}
        stream.readline()  # Remove header
        regexp5p3p = re.compile(r"3p\(|5p\(")
        regexpmir = re.compile(r"miRNA\|(hsa-.+-)(3p|5p)\(.+\)")
        for line in stream:
            spl = line.split(',')
            mir = spl[0]
            tss = "\t".join(spl[1:5])
            corr = eval(spl[5])
            pval = eval(spl[6])
            pls = eval(spl[10])
            if (re.search(regexp5p3p, mir) and corr >= correlation_thr and
                    pval <= pval_thr and
                    ((pls_thr and pls >= pls_thr) or (not pls_thr))):
                grp = re.match(regexpmir, mir)
                mir_name = grp.group(1)
                mir_arm = grp.group(2)
                if mir_name in asso:
                    if mir_arm in asso[mir_name]:
                        asso[mir_name][mir_arm].add(tss)
                    else:
                        asso[mir_name][mir_arm] = set([tss])
                else:
                    asso[mir_name] = {}
                    asso[mir_name][mir_arm] = set([tss])
        return asso


def get_all_mir_tss_asso(mir_tss_file):
    with open(mir_tss_file) as stream:
        asso = {}
        stream.readline()  # Remove header
        regexp5p3p = re.compile(r"3p\(|5p\(")
        regexpmir = re.compile(r"miRNA\|(hsa-.+-)(3p|5p)\(.+\)")
        for line in stream:
            spl = line.split(',')
            mir = spl[0]
            tss = "\t".join(spl[1:5])
            corr = eval(spl[5])
            pval = eval(spl[6])
            pls = eval(spl[10])
            if re.search(regexp5p3p, mir):
                grp = re.match(regexpmir, mir)
                mir_name = grp.group(1)
                mir_arm = grp.group(2)
                if mir_name in asso:
                    if mir_arm in asso[mir_name]:
                        asso[mir_name][mir_arm][tss] = (corr, pval, pls)
                    else:
                        asso[mir_name][mir_arm] = {}
                        asso[mir_name][mir_arm][tss] = (corr, pval, pls)
                else:
                    asso[mir_name] = {}
                    asso[mir_name][mir_arm] = {}
                    asso[mir_name][mir_arm][tss] = (corr, pval, pls)
        return asso


def print_inconsistant_tss(mir_5p_3p_tss_asso):
    for mir in mir_5p_3p_tss_asso:
        #if len(mir_5p_3p_tss_asso[mir]) == 2:
        #    print mir
        if (len(mir_5p_3p_tss_asso[mir]) == 2 and
                mir_5p_3p_tss_asso[mir]["5p"] != mir_5p_3p_tss_asso[mir]["3p"]):
            print "\n\n#######################################################"
            print "Different TSSs for {0}".format(mir.rstrip("-"))
            print "\n{0}5p associated to\n{1}".format(
                mir, '\n'.join(mir_5p_3p_tss_asso[mir]["5p"]))
            print "\n{0}3p associated to\n{1}".format(
                mir, '\n'.join(mir_5p_3p_tss_asso[mir]["3p"]))
            print "#######################################################"


def maximum_correlation(asso):
    maxi = -1
    best_tss = None
    for tss in asso:
        corr, _, _ = asso[tss]
        if corr > maxi:
            maxi = corr
            best_tss = tss
    return maxi, best_tss


def print_inconsistant_specific_tss(mir_5p_3p_tss_asso, min_corr_thr,
                                    max_corr_thr, pval_thr, pls_thr):
    for mir in mir_5p_3p_tss_asso:
        if len(mir_5p_3p_tss_asso[mir]) == 2:
            for tss in mir_5p_3p_tss_asso[mir]["5p"]:
                if tss in mir_5p_3p_tss_asso[mir]["3p"]:
                    corr_5p = mir_5p_3p_tss_asso[mir]["5p"][tss][0]
                    corr_3p = mir_5p_3p_tss_asso[mir]["3p"][tss][0]
                    pval_5p = mir_5p_3p_tss_asso[mir]["5p"][tss][1]
                    pval_3p = mir_5p_3p_tss_asso[mir]["3p"][tss][1]
                    max_corr_5p, max_tss_5p = maximum_correlation(
                        mir_5p_3p_tss_asso[mir]["5p"])
                    max_corr_3p, max_tss_3p = maximum_correlation(
                        mir_5p_3p_tss_asso[mir]["3p"])
                    if ((pval_5p < pval_thr and pval_3p < pval_thr) and
                        ((corr_5p > max_corr_thr and corr_3p < min_corr_thr)
                         or (corr_3p > max_corr_thr and
                             corr_5p < min_corr_thr))):
                        print "\n\n###########################################"
                        print "Different TSSs for {0}".format(mir.rstrip("-"))
                        print "{0}5p - {1} at {2:f},{3:f}".format(mir, tss,
                                                                  corr_5p,
                                                                  pval_5p)
                        print "{0}3p - {1} at {2:f},{3:f}".format(mir, tss,
                                                                  corr_3p,
                                                                  pval_3p)
                        print "Best for {0}5p - {1} with {2:f}".format(mir,
                                                                       max_tss_5p,
                                                                       max_corr_5p)
                        print "Best for {0}3p - {1} with {2:f}".format(mir,
                                                                       max_tss_3p,
                                                                       max_corr_3p)
                        print "###########################################"


def all_analysis(arguments):
    mir_5p_3p_tss_asso = get_all_mir_tss_asso(arguments.mir_tss_file)
    print_inconsistant_specific_tss(mir_5p_3p_tss_asso, arguments.min_corr_thr,
                                    arguments.max_corr_thr, arguments.pval_thr,
                                    arguments.pls_thr)


def all_arg_parsing(subparsers):
    help_str = "Analysis of all the miRNA-TSS associations"
    parser_t = subparsers.add_parser('A', help=help_str)
    parser_t.add_argument("-f", "--file", required=True, type=str,
                          dest="mir_tss_file", action="store",
                          help="miRNA-TSS file association")
    help_str = "Low Pearson correlation coefficient threshold (default:0.3)"
    parser_t.add_argument("-c", "--mincorrelation", required=False, type=float,
                          dest="min_corr_thr", action="store", default=0.3,
                          help=help_str)
    help_str = "High Pearson correlation coefficient threshold (default:0.3)"
    parser_t.add_argument("-C", "--maxcorrelation", required=False, type=float,
                          dest="max_corr_thr", action="store", default=0.7,
                          help=help_str)
    parser_t.add_argument('-p', '--pval', required=False, type=float,
                          dest="pval_thr", action="store", default=0.01,
                          help="P-value threshold (default:0.01)")
    parser_t.add_argument('-l', '--pls', required=False, type=float,
                          dest="pls_thr", action="store", default=None,
                          help="PLS threshold (default:None, no threshold)")
    parser_t.set_defaults(func=all_analysis)


def get_best(asso):
    best = None
    maxi = -1
    for tss in asso:
        corr, _, _ = asso[tss]
        if corr > maxi:
            maxi = corr
            best = tss
    return best, maxi


def print_best_different_tss(mir_5p_3p_tss_asso, diff_thr):
    for mir in mir_5p_3p_tss_asso:
        if "5p" in mir_5p_3p_tss_asso[mir] and "3p" in mir_5p_3p_tss_asso[mir]:
            best_tss_5p, corr_5p = get_best(mir_5p_3p_tss_asso[mir]["5p"])
            best_tss_3p, corr_3p = get_best(mir_5p_3p_tss_asso[mir]["3p"])
            if (best_tss_5p == best_tss_3p and
                    abs(corr_5p - corr_3p) > diff_thr):
                print "\n\n###################################################"
                print "Differences for best mir {0}".format(mir.rstrip('-'))
                print "{0}5p with {1} with {2:f}".format(mir, best_tss_5p,
                                                         corr_5p)
                print "{0}3p with {1} with {2:f}".format(mir, best_tss_3p,
                                                         corr_3p)
                print "#######################################################"


def best_analysis(arguments):
    mir_5p_3p_tss_asso = get_all_mir_tss_asso(arguments.mir_tss_file)
    print_best_different_tss(mir_5p_3p_tss_asso, arguments.diff_thr)


def best_arg_parsing(subparsers):
    help_str = "Analysis of the best miRNA-TSS associations differences"
    parser_t = subparsers.add_parser('B', help=help_str)
    parser_t.add_argument("-f", "--file", required=True, type=str,
                          dest="mir_tss_file", action="store",
                          help="miRNA-TSS file association")
    help_str = "Pearson correlation difference threshold (default:0.1)"
    parser_t.add_argument("-d", "--diff", required=False, type=float,
                          dest="diff_thr", action="store", default=0.1,
                          help=help_str)
    parser_t.set_defaults(func=best_analysis)


def diff_arg_parsing(subparsers):
    help_str = "Analysis of the miRNA-TSS association differences"
    parser_t = subparsers.add_parser('D', help=help_str)
    parser_t.add_argument("-f", "--file", required=True, type=str,
                          dest="mir_tss_file", action="store",
                          help="miRNA-TSS file association")
    help_str = "Pearson correlation difference threshold (default:0.1)"
    parser_t.add_argument("-d", "--diff", required=False, type=float,
                          dest="diff_thr", action="store", default=0.1,
                          help=help_str)
    parser_t.set_defaults(func=diff_analysis)


def print_diff_specific_tss(mir_5p_3p_tss_asso, diff_thr):
    for mir in mir_5p_3p_tss_asso:
        if len(mir_5p_3p_tss_asso[mir]) == 2:
            for tss in mir_5p_3p_tss_asso[mir]["5p"]:
                if tss in mir_5p_3p_tss_asso[mir]["3p"]:
                    corr_5p = mir_5p_3p_tss_asso[mir]["5p"][tss][0]
                    corr_3p = mir_5p_3p_tss_asso[mir]["3p"][tss][0]
                    max_corr_5p, max_tss_5p = maximum_correlation(
                        mir_5p_3p_tss_asso[mir]["5p"])
                    max_corr_3p, max_tss_3p = maximum_correlation(
                        mir_5p_3p_tss_asso[mir]["3p"])
                    if abs(corr_5p - corr_3p) > diff_thr:
                        print "\n\n###########################################"
                        print "TSS difference for {0}".format(mir.rstrip("-"))
                        print "{0}5p - {1} with {2:f}".format(mir, tss,
                                                              corr_5p)
                        print "{0}3p - {1} with {2:f}".format(mir, tss,
                                                              corr_3p)
                        print "Best for {0}5p - {1} with {2:f}".format(mir,
                                                                       max_tss_5p,
                                                                       max_corr_5p)
                        print "Best for {0}3p - {1} with {2:f}".format(mir,
                                                                       max_tss_3p,
                                                                       max_corr_3p)
                        print "###########################################"


def diff_analysis(arguments):
    mir_5p_3p_tss_asso = get_all_mir_tss_asso(arguments.mir_tss_file)
    print_diff_specific_tss(mir_5p_3p_tss_asso, arguments.diff_thr)


def top_analysis(arguments):
    mir_5p_3p_tss_asso = get_top_mir_tss_asso(arguments.mir_tss_file,
                                              arguments.corr_thr,
                                              arguments.pval_thr,
                                              arguments.pls_thr)
    print_inconsistant_tss(mir_5p_3p_tss_asso)


def top_arg_parsing(subparsers):
    help_str = "Analysis of the top miRNA-TSS associations"
    parser_t = subparsers.add_parser('T', help=help_str)
    parser_t.add_argument("-f", "--file", required=True, type=str,
                          dest="mir_tss_file", action="store",
                          help="top miRNA-TSS file association")
    parser_t.add_argument("-c", "--correlation", required=False, type=float,
                          dest="corr_thr", action="store", default=0.3,
                          help="Pearson correlation coefficient (default:0.3)")
    parser_t.add_argument('-p', '--pval', required=False, type=float,
                          dest="pval_thr", action="store", default=0.01,
                          help="P-value threshold (default:0.01)")
    parser_t.add_argument('-l', '--pls', required=False, type=float,
                          dest="pls_thr", action="store", default=None,
                          help="PLS threshold (default:None, no threshold)")
    parser_t.set_defaults(func=top_analysis)


def arg_parsing():
    descr = '''Analysis of the miRNA-TSS associations
       '''
    parser = argparse.ArgumentParser(description=descr,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(help="Choice of the analysis",
                                       title="Subcommands",
                                       description="Valid subcommands")
    top_arg_parsing(subparsers)
    all_arg_parsing(subparsers)
    best_arg_parsing(subparsers)
    diff_arg_parsing(subparsers)
    arguments = parser.parse_args()
    return arguments


###############################################################################
#                               MAIN
###############################################################################
if __name__ == "__main__":
    args = arg_parsing()
    args.func(args)
