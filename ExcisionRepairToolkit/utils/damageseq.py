#!/usr/bin/env python3
# coding: utf-8
"""
@author: Yiran Wu  wuyiran55@outlook.com
@last modified by: Yiran Wu
@file:damageseq.py
@time:1/9/23 11:19 AM
"""


from ExcisionRepairToolkit.core.bed import Bed
from ExcisionRepairToolkit.core.pipeline import Pipeline

def run(sample_bed, gene_list, sampleid):
    dam= Pipeline()
    sample = Bed(sample_bed)
    genelist = Bed(gene_list)

    #get readcount
    genelist.get_readcount(output=f'{sampleid}.readCount.txt')