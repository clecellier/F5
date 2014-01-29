""" Analysis of the expression difference between 5p and 3p. """

#!/usr/bin/python
#*-* coding: utf-8 *-*

import sys
import getopt
import re
from scipy.stats.stats import pearsonr, spearmanr, ks_2samp


def parse_mirs(fichier):
    """ Parse 5p/3p mirnas and store expression. Return dict of mirs. """

    with open(fichier) as stream:
        stream.readline()  # Header
        mir_dic = {}
        for line in stream:
            grp = re.search("^(.*)(3p|5p)", line)
            if grp:
                name = grp.group(1)
                arm = grp.group(2)
                expr = map(eval, line.split()[1:])
                if name in mir_dic:
                    mir_dic[name][arm] = expr
                else:
                    mir_dic[name] = {}
                    mir_dic[name][arm] = expr
        return mir_dic


def compute_correlations_ks(mirnas):
    """  Compute Pearson, Spearman, and KS for 5p/3p expr data of miRs. """

    print "miRNA\tPearson_r\tPearson_pval\tSpearman_rho\tSpearman_pval\t",
    print "KS_D\tKS_pval"
    for mir in mirnas:
        if len(mirnas[mir]) > 1:
            pears_r, pears_pval = pearsonr(mirnas[mir]["5p"],
                                           mirnas[mir]["3p"])
            spear_rho, spear_pval = spearmanr(mirnas[mir]["5p"],
                                              mirnas[mir]["3p"])
            ks_d, ks_pval = ks_2samp(mirnas[mir]["5p"], mirnas[mir]["3p"])
            mir_name = mir.rstrip('-')
            print "{0}\t{1:f}\t{2:f}\t{3:f}\t{4:f}\t{5:f}\t{6:f}".format(
                mir_name, pears_r, pears_pval, spear_rho, spear_pval, ks_d,
                ks_pval)


###############################################################################
#                               MAIN
###############################################################################
if __name__ == "__main__":
    usage = '''
    %s -f <expr file>
    ''' % (sys.argv[0])

    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:h")
    except getopt.GetoptError:
        sys.exit(str(getopt.GetoptError) + usage)

    expr_file = None
    for o, a in opts:
        if o == "-f":
            expr_file = a
        else:
            sys.exit(usage)
    if not expr_file:
        sys.exit(usage)

    mirs = parse_mirs(expr_file)
    compute_correlations_ks(mirs)
