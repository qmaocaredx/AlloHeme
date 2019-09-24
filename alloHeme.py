#!/usr/bin/env python3

import argparse
import collections
import csv
import logging
import math
import pandas

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
        "zeroCutoff": 0.0007848609,
        "upperLimitOfDetection": 0.15,
        "maxBackground": 0.005,
        "refAltCutoff": 0.95,
        "subtractBackground": 1,
        "backgroundMultiplier": 3.345202,
        "minRecipHetCutoff": 0.10,
        "maxRecipHetCutoff": 0.25,
        "minCutoff": 0,
        "maxCutoff": 0.95,
        "numSnps": 266
    }

    # Initializes the analysis by extracting the command-line arguments
    def __init__(self, args):
        self.pileup_summary_filename = args.pileup_summary
        self.amplicons_filename = args.amplicons
        self.config_filename = args.config
        self.run_name = args.run_name
        self.sequencer = args.sequencer
        self.output_by_position_filename = args.output_by_position
        self.output_by_sample_filename = args.output_by_sample
        self.output_by_snp_filename = args.output_by_snp

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

    def PUGT.expandGT (GTs, prior = NULL):
    # 1) expand partially unknown genotypes
    # 2) compute the posterior probability
    D = ncol(GTs)
    contributorNames = colnames(GTs)
    GTexp = matrix(0, nrow(GTs), 3^D)
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