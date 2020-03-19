#!/usr/bin/env python

# Copyright (c) 20120 Tobias Neumann

from __future__ import print_function
import sys

import Bio
import Bio.SeqIO
from Bio.Seq import Seq

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, SUPPRESS

 # Info
usage = "CRISPR2019 - design guides for miRNAs"
version = "1.0"

# Main Parsers
parser = ArgumentParser(description=usage)

parser.add_argument("-f", "--fasta", type=str, required=True, dest="fastaFile", help="Fasta file of hCov MSA")
parser.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", help="Reference fasta file")
parser.add_argument("-l", "--leftMask", type=int, required=False, default=130, dest="leftMask", help="Basepairs masked from left")
parser.add_argument("-b", "--rightMask", type=int, required=False, default=50, dest="rightMask", help="Basepairs masked from right")

args = parser.parse_args()

covidSequencesConservative = dict()
covidSequencesLenient = dict()

referenceID = ""
referenceSeq = ""

for record in Bio.SeqIO.parse(args.referenceFile, 'fasta'):
    referenceID = str(record.id)
    referenceSeq = str(record.seq)

for record in Bio.SeqIO.parse(args.fastaFile, 'fasta'):

    id = str(record.id)
    seq = str(record.seq)

    i = 0
    for base in seq:
        if i not in covidSequencesConservative:
            covidSequencesConservative[i] = dict()
        covidSequencesConservative[i][base] = 0
        if i not in covidSequencesLenient:
            covidSequencesLenient[i] = dict()
        if base != "N":
            covidSequencesLenient[i][base] = 0
        i += 1

conservativeConsensus = ""
lenientConsensus = ""

with open("hCov_conservative_mutated_sites.txt","w") as f:

    for pos in covidSequencesConservative:
        if len(covidSequencesConservative[pos]) == 1:
            conservativeConsensus += next(iter(covidSequencesConservative[pos]))
        else :
            if (pos >= args.leftMask and pos <= (len(referenceSeq) - args.rightMask)):
                print(str(pos) + "\t" + referenceSeq[pos] + ">" + "".join(covidSequencesConservative[pos].keys()),file=f)
            conservativeConsensus += "N"

with open("hCov_lenient_mutated_sites.txt","w") as f:

    for pos in covidSequencesLenient:
        if len(covidSequencesLenient[pos]) == 1:
            lenientConsensus += next(iter(covidSequencesLenient[pos]))
        else :
            if (pos >= args.leftMask and pos <= (len(referenceSeq) - args.rightMask)):
                print(str(pos) + "\t" + referenceSeq[pos] + ">" + "".join(covidSequencesLenient[pos].keys()),file=f)
            lenientConsensus += "N"

with open("hCov_conservative_consensus.fa","w") as f:
    print(">" + referenceID,file=f)
    print(referenceSeq,file=f)
    print(">hCov_conservative_consensus",file=f)
    print(conservativeConsensus,file=f)

with open("hCov_lenient_consensus.fa","w") as f:
    print(">" + referenceID,file=f)
    print(referenceSeq,file=f)
    print(">hCov_lenient_consensus",file=f)
    print(lenientConsensus,file=f)
