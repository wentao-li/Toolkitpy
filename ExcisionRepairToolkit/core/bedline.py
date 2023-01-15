#!/usr/bin/env python3
# coding: utf-8
"""
@author: Yiran Wu  wuyiran55@outlook.com
@last modified by: Yiran Wu
@file:bedline.py
@time:12/19/22 10:52 AM
"""

import random

class Bedline:

    def __init__(
            self,
            line,
            sep='\t',
            col_idx={'chromosomeCol': 1, 'startCol': 2, 'endCol': 3, 'nameCol': 4, 'scoreCol': 5, 'strandCol': 6}):
        self.line = line.strip()
        self.separator = sep
        self.col_idx = col_idx

    def getLine(self):
        return self.line

    def fields(self):
        return self.line.split(self.separator)

    def chromosome(self,):
        return self.fields()[self.col_idx['chromosomeCol'] - 1]

    def start(self, startCol=2):
        return int(self.fields()[self.col_idx['startCol'] - 1])

    def length(self):
        return self.end() - self.start()

    def end(self, endCol=3):
        return int(self.fields()[self.col_idx['endCol'] - 1])

    def name(self, nameCol=4):
        return self.fields()[self.col_idx['nameCol'] - 1]

    def score(self, scoreCol=5):
        return int(self.fields()[self.col_idx['scoreCol'] - 1])

    def strand(self, strandCol=6):
        return self.fields()[self.col_idx['strandCol'] - 1]

    def midpoint(self, randomness=False):
        if not randomness:
            addition = 0
        else:
            addition = random.sample([0, 1], 1)[0]
        return int((self.start() + self.end() + addition) / 2)

    def newline(self, start, end):
        newList = self.fields()
        newList[1] = str(start)
        newList[2] = str(end)
        return self.separator.join(newList)

