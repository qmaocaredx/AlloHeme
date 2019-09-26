#!/usr/bin/env python3

import argparse
import collections
import csv
import logging
import math
import pandas
from scipy.special import betaln

#   alloHeme.py is a Python script for AlloHeme's dd-cfDNA fraction
#   calculation.
#
#   This program operates on a single sample at a time, instead of
#   processing multiple samples.
#
# The following command-line arguments are REQUIRED:
#
# --pileup-summary <filename>
#
#   Input file containing pileup summary for a single sample output by
#   secondary analysis.
#
# --config <filename>
#
#   Input file containing pipeline configuration parameters.
#

VERSION = "V1"

# Parses a string as either an int or a float.
def parse_num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)


class DdDNACalculation:

    # Default analysis parameters, used as fallback in case they are
    # missing from the config file.
    default_params = {
        "coverageCutoff": 2000,
        "numSnps": 405
    }

    # Initializes the analysis by extracting the command-line arguments
    def __init__(self, args):
        self.pileup_summary_filename = args.pileup_summary
        self.config_filename = args.config

    # Parses the contents of the configuration file and stores the
    # configuration in a dict, with defaults from the default_params
    # dict.

    def geneticPrior (relationship = "unrelated" , pi = 0.5):
        switcher={
            'parent', 'child': (1-pi)^3, pi * (1-pi)^2, 0, pi * (1-pi)^2, pi * (1-pi), pi^2 * (1-pi), 0, pi^2 * (1-pi), pi^3,
            'sibling': (1-pi)^2*(1-pi/2)^2, pi * (1-pi)^2 * (1-pi/2), pi^2 * (1-pi)^2/4, pi * (1-pi)^2*(1-pi/2), pi * (1-pi) * (1+pi-pi^2), pi^2 * (1-pi) * (1/2+pi/2), 
                            pi^2 * (1-pi)^2/4, pi^2 * (1-pi) * (1/2 + pi/2), pi^2*(1/2+pi/2)^2,
            'uncle','aunt','nephew','niece': (1-pi)^3*(1-pi/2), pi * (1-pi)^2 * (3/2-pi), pi^2 * (1-pi)^2/2,  pi * (1-pi)^2*(3/2-pi), pi * (1-pi) * (1/2+2*pi-2*pi^2), pi^2 * (1-pi) * (1/2+pi), 
                                    pi^2 * (1-pi)^2/2, pi^2 * (1-pi) * (1/2 + pi), pi^3*(1/2+pi/2)
            'unrelated': (1-pi)^4, 2*pi * (1-pi)^3, pi^2 * (1-pi)^2, 2*pi * (1-pi)^3, 4*pi^2 * (1-pi)^2, 2*pi^3 * (1-pi), pi^2 * (1-pi)^2, 2*pi^3 * (1-pi), pi^4
        }
        return switcher.get(i,"Invalid relationship")

    def expectedAF ( GT, lambda, b):
        G = ((1- lambda) * GT + lambda * (2- GT)) / 2
        p = np.dot(G,b)
        q = np.dot(1-G, b)
        return np.concatenate((p, q), axis=1)
        #return(cbind(p = G %*% t(t(b)), q = (1-G) %*% t(t(b)))) ### Why transpose of transpose?

    def PUGT.Prior (contributorNames = contributors, relationships = 'unrelated', pi = 0.5):
        # compute the prior for 3D their corresponding prior probabilities. 
        # Some of the probabilities can be 0
        
        allGT <- allGTID <- c()
        prior1 = c(pi^2, 2*pi*(1-pi), (1-pi)^2)
        allprior = c()

        for (d in 1:length(contributorNames)):
            allGTID = cbind(allGTID, paste(contributorNames[d], 0:2, sep = ':'))
            allGT = cbind(allGT, 0:2)
            #allprior = cbind(allprior, prior1)
            allprior = get.geneticRelatedness.prior(relationship="child",pi=0.5)

        colnames(allGT) = contributorNames
        allGTIDs = apply(as.matrix(expand.grid(as.data.frame(allGTID))), 1, FUN = paste, collapse = '|')
        allpriors = cbind(prior = apply(as.matrix(expand.grid(as.data.frame(allprior))), 1, FUN = prod), as.matrix(expand.grid(as.data.frame(allGT))))
        rownames(allpriors) = allGTIDs;
        return(allpriors)

    def PUGT.expandGT (GTs, prior = NULL):
        # 1) expand partially unknown genotypes
        # 2) compute the posterior probability
        D = np.size(GTs,1) #ncols
        contributorNames = list(GTs) #probably GTs will be a pandas object #shortcut to get column names of df
        GTexp = numpy.zeros(len(GTs.index), 3**D) #GTexp = matrix(0, nrow(GTs), 3^D)
        allGTs = names(prior)
        print("in Expand GT, prior")
        print (prior)
        colnames(GTexp) = allGTs
        rownames(GTexp) = rownames(GTs)
        for (g in 1:nrow(GTs)){
            G1 = GTs[g,] # expand degenerate genotypes for each loci
            allGT.ob = list()
            for (d in 1:D){
            if (is.na(G1[d])){
                Gs = 0:2
            }else{
                Gs = G1[d]
            }
            allGT.ob[[d]] = paste(contributorNames[d], Gs, sep = ':')
            }
            allGT.obs = apply(as.matrix(expand.grid(allGT.ob)), 1, FUN = paste, collapse = '|')
            priorNormalize = prior[allGT.obs];
            priorNormalize = priorNormalize/sum(priorNormalize);
            GTexp[g, allGT.obs] = priorNormalize
        return(GTexp)

    def BB (n1, n2, a, b, r):
        return (betaln(n1 + ni*(1+r)/2*a, n2+ni*(1+r)/2*b) - betaln(ni*(1+r)/2*a, ni*(1+r)/2*b))/log(2)

    def PUGT.LL (beta, r, lambda, expandGT, prior, ni, n1, n2){
        #log 2 transform
        ### compute the expected allele fraction for 3^D possible genotype. Return a vector of length 3^D
        eAF = expectedAF(GT = prior[,2:ncol(prior)], lambda = lambda, b = beta) 
        n = n1 + n2
        mu.n = mean(n)
       
        ### add the constant coefficents to get correct likelihood values
        coeff <- (math.lgamma(n+1)-math.lgamma(n1+1)-math.lgamma(n2+1))/log(2)  
        
        ll.mat = matrix(-Inf, nrow = length(n1), ncol = nrow(eAF))
        for i in 1:nrow(eAF):
            ll.mat[,i] = BB(n1, n2, eAF[i,1], eAF[i,2], r)
        
        colnames(ll.mat) = rownames(eAF)
        rownames(ll.mat) = rownames(expandGT)
        
        ll.mat[expandGT==0] = -Inf
        ### add the genetic prior
        ll.mat[expandGT>0] = log2(expandGT[expandGT>0]) + ll.mat[expandGT>0] 
        ll = rowLogSumExp(ll.mat);
        ll = ll + coeff
        return(sum(ll))
    
    def argmin (x, func):
        minFunc = np.inf
        if !x:
            x = np.array(x, nrow = length(x), ncol = 1)
        for i in x:
            xx = x[i,]
            y = func(xx)
            if y < minFunc:
                minX = xx;
                minFunc = y;
        return(minX)

    def mapping.bias (dat):
        r = dat.Ref_Allele_Counts / (dat.Ref_Allele_Counts + dat.Alt_Allele_Counts)
        cutoffs = [0.05, 0.1, 0.25, 0.5]
        braw <- br <- cutoffs
        names(braw) <- names(br) <- cutoffs
        for c1 in cutoffs:
            braw[paste(c1)] = df.sum(dat.Ref_Allele_Counts[r<c1], skipna = True)/df.sum(dat.Alt_Allele_Counts[r > 1- c1], skipna = True)
            br[paste(c1)] = df.sum(r[r<c1] , skipna = True)/df.sum(1-r[r > 1- c1], skipna = True)

        braw = c(braw, allRef_Alt = df.sum(dat.Ref_Allele_Counts , skipna = True)/df.sum(dat.Alt_Allele_Counts , skipna = True), 
                lambda.Alt2Ref = df.sum(dat.Ref_Allele_Counts[r < 0.25] , skipna = True) /df.sum(dat.Alt_Allele_Counts[r < 0.25] , skipna = True),
                lambda.Ref2Alt = df.sum(dat.Alt_Allele_Counts[r > 0.75] , skipna = True)/df.sum(dat.Ref_Allele_Counts[r > 0.25] , skipna = True))
        return(list(braw = braw, br = br))
    

    # Runs the analysis
    def run(self):
        logging.info("Running tertiary analysis...")
        self.read_config()
        self.read_experiment_name()
        self.read_amplicons()
        self.process_pileup_summary()
        self.get_sample_background()
        self.compute_snp_table()
        self.estimate_cfDNA_fraction()
        self.write_sample_output_file()
        self.write_snp_output_file()

def main():
    logging.basicConfig(level=logging.INFO)
    logging.info("Running tertiary analysis version %s", VERSION)

    parser = argparse.ArgumentParser()
    parser.add_argument("--pileup-summary", required=True, metavar="FILENAME",
                        help="Input file containing pileup summary for a single sample")
    parser.add_argument("--amplicons", required=True, metavar="FILENAME",
                        help="Input file containing annotation of SNPs and amplicons")
    parser.add_argument("--config", required=True, metavar="FILENAME",
                        help="Input file containing pipeline configuration parameters")
    parser.add_argument("--run-name", required=True, metavar="STRING",
                        help="Name of the run which included the sample")
    parser.add_argument("--sequencer", required=True, metavar="STRING",
                        help="Sequencer on which the sample was run")
    parser.add_argument("--output-by-position", required=True, metavar="FILENAME",
                        help="Output file containing data for all alleles at every position of every amplicon")
    parser.add_argument("--output-by-sample", required=True, metavar="FILENAME",
                        help="Output file containing summary data for the sample")
    parser.add_argument("--output-by-snp", required=True, metavar="FILENAME",
                        help="Output file containing data for SNP positions only")

    args = parser.parse_args()
    logging.info("Command-line arguments:")
        logging.info("\tPileup summary file: %s", args.pileup_summary)
    logging.info("\tAmplicons file: %s", args.amplicons)
    logging.info("\tConfig file: %s", args.config)
    logging.info("\tRun name: %s", args.run_name)
    logging.info("\tSequencer: %s", args.sequencer)
    logging.info("\tOutput by position file: %s", args.output_by_position)
    logging.info("\tOutput by sample file: %s", args.output_by_sample)
    logging.info("\tOutput by SNP file: %s", args.output_by_snp)

    pandas.set_option('display.width', 230)
    pandas.set_option('display.max_columns', 30)

    analysis = TertiaryAnalysis(args)
    analysis.run()

if __name__ == "__main__":
    main()