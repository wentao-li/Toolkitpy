#!/usr/bin/env python3
# coding: utf-8
"""
@author: Yiran Wu  wuyiran55@outlook.com
@last modified by: Yiran Wu
@file:bed.py
@time:12/13/22 10:27 AM
"""

import sys
import os
import sys
import unittest



## bed 的每一行是object
## 提取数据操作，长度，分数，列名等

## bed文件是object
## 读取，切窗口等操作

from ExcisionRepairToolkit.core.bedline import Bedline


class Bed:

    def __init__(self, file_path, columns=["chrom", "chromStart", "chromEnd", "name", "score", "strand"]):
        self.file = file_path
        self.columns = columns

    def getTotalRegionLength(self):
        totalLength = 0
        filein = open(self.file, 'r')
        for line in filein:
            ll = line.split('\t')
            beg = int(ll[1])
            end = int(ll[2])
            regionLength = end - beg
            totalLength += regionLength
        return totalLength

    def getColumnNumber(self):
        firstLine = open(self.file, 'r').readline().strip()
        return len(firstLine.split('\t'))

    def getAverageLength(self):
        totalLength = self.getTotalRegionLength()
        hitNum = self.getHitNum()
        return float(totalLength) / hitNum

    def getHitNum(self):
        hitNum = 0
        filein = open(self.file, 'r')
        for line in filein:
            hitNum += 1
        return hitNum

    def fixRange(self, side, length):
        filein = open(self.file, 'r')
        for line in filein:
            print(bedLine2fixedRangedLine(line, side, length))
        return self

    def lengthDistribution(self):
        filein = open(self.file, 'r')
        myDict = {}
        for line in filein:
            ll = line.split('\t')
            start = int(ll[1])
            end = int(ll[2])
            length = end - start
            if length in myDict.keys():
                myDict[length] += 1
            else:
                myDict[length] = 1
        sortedKeys = sorted(myDict.keys())
        for i in range(min(sortedKeys), max(sortedKeys) + 1):
            if not i in myDict.keys():
                count = 0
            else:
                count = myDict[i]
            print(str(i) + '\t' + str(count))

    def printLengths(self):
        filein = open(self.file, 'r')
        for line in filein:
            ll = line.split('\t')
            start = int(ll[1])
            end = int(ll[2])
            length = end - start
            print(length)

    def makeWindowsPerLine(self, interval, windowSize):
        start = interval[0]
        end = interval[1]
        length = end - start
        windowNumber = length / windowSize
        newIntervals = []
        for i in range(windowNumber):
            newIntervals.append([start + i * windowSize, start + i * windowSize + windowSize])
        return newIntervals

    def makeWindows(self, windowSize, noShortFlag=False):
        filein = open(self.file, 'r')
        for line in filein:
            ll = line.strip().split('\t')
            interval = [int(ll[1]), int(ll[2])]
            newIntervals = self.makeWindowsPerLine(interval, windowSize)
            for newInterval in newIntervals:
                newLineList = list(ll)
                newLineList[1] = str(newInterval[0])
                newLineList[2] = str(newInterval[1])
                printFlag = True
                if noShortFlag:
                    newIntervalLength = newInterval[1] - newInterval[0]
                    if newIntervalLength < windowSize:
                        printFlag = False
                if printFlag:
                    print('\t'.join(newLineList))

    def read(self):
        filein = open(self.file, 'r')
        for line in filein:
            yield (Bedline(line))

    def removeNeighbors(self, distance):
        def areNeighbors(previousLine, line, distance):
            bedLine = Bedline(line)
            previousBedLine = Bedline(previousLine)
            if bedLine.chromosome() != previousBedLine.chromosome():
                return False
            start = bedLine.start()
            previousEnd = previousBedLine.end()
            if abs(start - previousEnd) < distance:
                return True
            else:
                return False

        previousLine = "NA\t" + str(-2 * distance) + "\t" + str(-1 * distance)
        filein = open(self.file, 'r')
        start = True
        dontPrintNextline = False
        for line in filein:
            linesAreNeighbors = areNeighbors(previousLine, line, distance)
            if linesAreNeighbors == False:
                if not start:
                    if not previousLineWasNeigbor:
                        print(previousLine.strip())
                dontPrintNextline = False
                previousLineWasNeigbor = False
            else:
                previousLineWasNeigbor = True
            previousLine = line
            start = False
        if linesAreNeighbors == False:
            print(line.strip())

    def divide_transcript(self, sampleid, up_length=6000, up_bin_count=25, tx_bin_count=100):
        from statistics import mean
        up_bin_length = int(up_length / up_bin_count)
        total_number_of_bins = tx_bin_count + up_bin_count * 2
        output = sampleid + "_bin_" + str(total_number_of_bins) + ".bed"
        OUT = open(output, "w")

        from ..io.reader import read_bed
        input_bed = read_bed(self.file)
        for line in input_bed:
            chrom = line.chromosome()
            tx_start = line.start()
            tx_end = line.end()
            strand = line.strand()
            tx_name = line.name()
            tx_length = tx_end - tx_start
            # upstream
            up_starts = []
            up_start = tx_start - up_length
            up_starts.append(up_start)
            for i in range(up_bin_count - 1):
                up_starts.append(up_starts[-1] + up_bin_length)
            up_ends = up_starts[1:]
            up_ends.append(tx_start)

            # transcript
            ideal_bin_size = tx_length / tx_bin_count
            bin_size = int(ideal_bin_size)
            tx_starts = [tx_start]
            bin_sizes = [ideal_bin_size]
            for i in range(1, tx_bin_count):
                if mean(bin_sizes) < ideal_bin_size:
                    tx_starts.append(tx_starts[-1] + int(bin_size) + 1)
                    bin_sizes.append(tx_starts[-1] - tx_starts[-2])
                else:
                    tx_starts.append(tx_starts[-1] + int(bin_size))
                    bin_sizes.append(tx_starts[-1] - tx_starts[-2])
            tx_ends = tx_starts[1:]
            tx_ends.append(tx_end)

            # downstream
            down_starts = [tx_end + x * up_bin_length for x in range(up_bin_count)]
            down_ends = [tx_end + (x + 1) * up_bin_length for x in range(up_bin_count)]

            starts = up_starts + tx_starts + down_starts
            ends = up_ends + tx_ends + down_ends

            if strand == '-':
                starts.reverse()
                ends.reverse()

            for i in range(total_number_of_bins):
                OUT.write('{chrom}\t{tx_start}\t{tx_end}\t{tx_name}\t{bin_no}\t{strand}\n'.format(
                    chrom=chrom,
                    tx_start=starts[i],
                    tx_end=ends[i],
                    bin_no=(i + 1),
                    tx_name=tx_name,
                    strand=strand,
                ))
        OUT.close()
        return output

    def remove_intersection(self, sampleid, gap=6000):
        from ..io.reader import read_bed
        output = sampleid + "_rmintersection.bed"
        op = open(self.file, "r")
        geneDict = {}
        input_bed = read_bed(self.file)
        abandomList = []
        for line in input_bed:
            chrom = line.chromosome()
            start = line.start()
            end = line.end()
            key = chrom + "-" + str(start) + "-" + str(end)
            gene = line.name()
            strand = line.strand()
            length = line.length()
            geneDict[key] = line.fields()
        #for key, value in geneDict.items():
            if key in abandomList:
                continue
            #chrom, start, end, gene, length, strand=value
            #chrom, start, end, strand, length, gene = value
            start = int(start) - gap if int(start) - gap >= 0 else 0
            end = int(end) + gap
            tag = False

            input_bed2 = read_bed(self.file)
            for line2 in input_bed2:
                chrom2 = line2.chromosome()
                start2 = line2.start()
                end2 = line2.end()
                key2 = chrom2 + "-" + str(start2) + "-" + str(end2)
                gene2 = line2.name()
                strand2 = line2.strand()
                length2 = line2.length()
            #for key2, value2 in geneDict.items():
                if key2 == key:
                    continue
                #chrom2, start2, end2, gene2, length2, strand2=value2
                #chrom2, start2, end2, strand2, length2, gene2=value2
                #start2, end2 = int(start2),int(end2)
                if chrom2 == chrom:
                    if end2 < start or start2 > end:
                        pass
                    else:
                        tag = True
                        abandomList.append(key2)
            if tag:
                abandomList.append(key)
        ot = open(output, "w")
        for k, v in geneDict.items():
            if k not in abandomList:
                start = v[1]
                v[1] = str(v[1])
                end = v[2]
                v[2] = str(v[2])
                ot.write('{v}\n'.format(v='\t'.join(v)))
        op.close()
        ot.close()
        return output

    def get_readcount(self, sampleid):
        output = sampleid + "_readCount.txt"
        command = 'grep -c "^" ' + f'{self.file} > {output}'
        os.system(command)
        return output

def bedpeLine2bedLine(bedpeLine):
    # Converts paired bed line to single fragment by getting the maximum range between to genomic locations.
    ll = bedpeLine.strip().split('\t')
    chr1 = ll[0]
    beg1 = int(ll[1])
    end1 = int(ll[2])
    chr2 = ll[3]
    beg2 = int(ll[4])
    end2 = int(ll[5])
    name = ll[6]
    score = ll[7]
    strand1 = ll[8]
    strand2 = ll[9]

    if chr1 != chr2:
        sys.exit('Chromosome is not identical between mates! Exiting...')

    if strand1 == '+':
        beg = beg1
        end = end2
    elif strand1 == '-':
        beg = beg2
        end = end1

    newLine = chr1 + \
              '\t' + str(beg) + \
              '\t' + str(end) + \
              '\t' + name + \
              '\t' + score + \
              '\t' + strand1

    if len(ll) > 10:
        newLine += '\t' + '\t'.join(ll[10:])

    return newLine


def bedLine2fixedRangedLine(bedLine, side, length, strandFlag=True):
    ll = bedLine.strip().split('\t')
    beg = int(ll[1])
    end = int(ll[2])
    strand = ll[5]
    if (strandFlag == True and strand == '+') or (strandFlag == False):
        if side == 'l':
            newEnd = min(beg + length, end)
            newBeg = beg
        elif side == 'r':
            newBeg = max(end - length, beg)
            newEnd = end
    elif strandFlag == True and strand == '-':
        if side == 'l':
            newBeg = max(end - length, beg)
            newEnd = end
        elif side == 'r':
            newEnd = min(beg + length, end)
            newBeg = beg
    if side != 'l' and side != 'r':
        sys.exit('Error: side must be either l or r. Exiting...')
    ll[1] = str(newBeg)
    ll[2] = str(newEnd)
    return '\t'.join(ll)
