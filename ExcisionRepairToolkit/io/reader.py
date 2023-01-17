#!/usr/bin/env python3
# coding: utf-8
"""
@author: Yiran Wu  wuyiran55@outlook.com
@last modified by: Yiran Wu
@file:reader.py
@time:12/18/22 5:11 PM
"""

from ExcisionRepairToolkit.core.bed import Bed
from ExcisionRepairToolkit.core.bedline import Bedline
import pandas as pd


def read_bed(file_path,):
    bed = Bed(file_path)
    filein = open(file_path, 'r')
    for line in filein:
        yield (Bedline(line))


def read_samplelist(sample_list):
    from ExcisionRepairToolkit.core.sample import Sample
    f = open(sample_list, "r")
    for line in f.readlines():
        ll = line.strip().split("\t")
        id = ll[0]
        fqs = ll[1]
        sample = Sample(id=id, fqs=fqs)
        return sample