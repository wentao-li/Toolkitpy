#!/usr/bin/env python3
# coding: utf-8
"""
@author: Yiran Wu  wuyiran55@outlook.com
@last modified by: Yiran Wu
@file:pipeline.py
@time:12/13/22 10:23 AM
"""
import os

from ExcisionRepairToolkit.core.bedline import Bedline
from ExcisionRepairToolkit.core.bed import Bed

class XRseqPipeline(object):
    def __init__(self, sample_bed, gene_list, sampleid):
        self.sampleid = sampleid
        self.sample_bed = Bed(sample_bed)
        self.gene_list = Bed(gene_list)
        self.result = dict()
        self.process = []

    def run_commands(self, command, output, step):
        os.system(command)
        self.process.append(step)
        self.result[step] = output

    def cut_bylength(self, bed_file, sampleid):
        from ..io.reader import read_bed
        input_bed = read_bed(bed_file)
        output = f"{sampleid}_filtered.bed"
        command = "awk '($5>300 && ($3-$2)>3000){print $0}' " + bed_file + ">" + output
        self.run_commands(command, output=output, step="cut_bylength")
        return output


    def get_readcount(read_count_file): ## 重新计算最后一列得到 rpkm
        '''

        :param read_count_file: "11173882"
        :return:
        '''
        with open(read_count_file) as f:
            read_count = int(f.readline().strip())
            nf = 10 ** 9 / read_count
        return nf

    def get_TSandNTS(
            self,
            sample_bed,
            gene_list,
            output,
            module_load="module load BEDTools/2.30.0-GCC-8.3.0",
            param="-wa -c -S -F 0.5",

    ):
        ## generate *.zip
        #output = sampleid + "_" + type + ".bed"
        command = f"bedtools intersect -c -a {gene_list} -b {sample_bed} {param} > {output}"
        if module_load is not None:
            command = module_load + "\n" + command
        self.run_commands(command, output, "get_TSandNTS")
        return output

    @staticmethod
    def get_rpkm(read_count_file, ts_bed, nts_bed, sampleid): ## 重新计算最后一列得到 rpkm
        '''

        :param read_count_file: "11173882"
        :param ts_bed: "# chr1    303262  303502  AP006222.2  1   -   0"
        :param nts_bed: "# chr1    303262  303502  AP006222.2  1   -   0"
        :return:
        '''
        with open(read_count_file) as f:
            read_count = int(f.readline().strip())
            nf = 10 ** 9 / read_count
        for intersect in [ts_bed, nts_bed]:
            output = sampleid + "rpkm_" + intersect
            out = open(output, "w")
            with open(intersect) as f:
                for line in f:
                    bed_line = Bedline(line)
                    #count = bed_line.fields()[-1]
                    bed_line.fields()[-1] = (bed_line.fields()[-1] * nf) / (bed_line.end()-bed_line.start())
                    #with open(output, 'a') as d:
                    out.write(bed_line.newline(start=bed_line.start(), end=bed_line.end()))
                    #out.write(line.bed_line())
        return output

    @staticmethod
    def get_rpgcm(fasta, read_count_file, intersect_count_file, sampleid):
        # read count file
        with open(read_count_file) as f:
            read_count = int(f.readline().strip())
            nf = 10 ** 9 / read_count
        ## read tad count

        output = sampleid + "_rpkgcm.txt"
        gc_dict = {}
        pre_region = None
        f = open(fasta, "r")
        for i in f.readlines():
            i = i.strip()
            if i.startswith(">"):
                pre_region = i
                gc_dict[pre_region] = 0
                continue
            count_sub = {}
            for sub in i:
                count_sub[sub] = i.count(sub)
            for letter in ['G', 'C', 'g', 'c']:
                if letter in count_sub.keys():
                    gc_dict[pre_region] += count_sub[letter]
        out = open(output, "w")
        with open(intersect_count_file) as f:
            for line in f:
                region = f"{bed_line.chromosome()}:{bed_line.start()}-{bed_line.end()}"
                gc_count = gc_dict[region]
                bed_line = Bedline(line)
                bed_line.fields()[-1] = (bed_line.fields()[-1] * nf) / gc_count
                out.write(bed_line.newline(start=bed_line.start(), end=bed_line.end()))
        out.close()


    @staticmethod
    def get_rpkm_mean(bed_file, output):
        '''
        :param bed_file: "chr1	35313	35553	WASH7P	1	-	0.0"
        :return:
        '''

        import numpy as np
        import pandas as pd
        #output = sampleid + "_mean.txt"
        rp = pd.read_csv(bed_file, sep="\t", header=None)
        rp.columns = ['chr', 'start', 'end', 'gene', 'bin', "strand", 'rpkm']
        rp.groupby(by='bin')['rpkm'].agg([np.mean]).to_csv(output,sep=" ", header=None)
        return output
