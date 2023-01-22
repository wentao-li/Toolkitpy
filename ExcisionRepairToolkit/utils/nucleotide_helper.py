#!/usr/bin/env python3
# coding: utf-8
"""
@author: Yiran Wu  wuyiran55@outlook.com
@last modified by: Yiran Wu
@file:nucleotide_helper.py
@time:1/14/23 11:40 AM
"""
import re


def nucleotide_dis2(input, sampleid, output):
    f = open(input, 'r')
    n_dis = {}
    for LINES in f:
        if LINES[0] != ">" and len(LINES) == 13:
            for i in range(0, len(LINES)):
                if LINES[i] == 'A' or LINES[i] == 'a':
                    key = 'A{0}'.format(i + 1)
                    try:
                        n_dis[key] = n_dis[key] + 1
                    except:
                        n_dis[key] = 1
                elif LINES[i] == 'T' or LINES[i] == 't':
                    key = 'T{0}'.format(i + 1)
                    try:
                        n_dis[key] = n_dis[key] + 1
                    except:
                        n_dis[key] = 1
                elif LINES[i] == 'C' or LINES[i] == 'c':
                    key = 'C{0}'.format(i + 1)
                    try:
                        n_dis[key] = n_dis[key] + 1
                    except:
                        n_dis[key] = 1
                elif LINES[i] == 'G' or LINES[i] =='g':
                    key = 'G{0}'.format(i + 1)
                    try:
                        n_dis[key] = n_dis[key] + 1
                    except:
                        n_dis[key] = 1
        else:
            continue
    f.close()

    out = open(output, 'w')
    for X, V in n_dis.items():
        out.write('{0}\t{1}\t{2}\n'.format(X[0], X[1:3], V))
    out.close()


def get_nucleotide_abundance(input, output, sequence_length=None,  nucleotide_order='ATGC', percentageFlag=True):
    from ..core.fasta import Fasta
    fasta = Fasta(input)
    nucleotide_abundance_dict = fasta.getNucleotideAbundance(sequence_length)
#    out = open(output, "w")
    fasta.writeNucleotideAbundanceTable(nucleotide_abundance_dict, output, nucleotide_order, percentageFlag)
#    out.close()

