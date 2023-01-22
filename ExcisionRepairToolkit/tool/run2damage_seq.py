#!/usr/bin/env python3
# coding: utf-8
"""
@author: Yiran Wu  wuyiran55@outlook.com
@last modified by: Yiran Wu
@file:run2damage_seq.py
@time:1/14/23 10:02 AM
"""

import sys
sys.path.append("/work/wllab/yiran/01.scripts/Toolkitpy/")
from ExcisionRepairToolkit.core.sample import Sample
from ExcisionRepairToolkit.core.bed import Bed
from ExcisionRepairToolkit.core.pipeline import XRseqPipeline
from ExcisionRepairToolkit.plots.ts_nts import plot_ts_nts
from ExcisionRepairToolkit.io.writer import check_dir


def runner(sample_list, work_dir):
    #allsample = read_samlplelist(sample_list)
    from ExcisionRepairToolkit.core.sample import Sample
    f = open(sample_list, "r")
    for line in f.readlines():
        ll = line.strip().split("\t")
        sampleid = ll[0]
        fqs = ll[1]
        outdir = work_dir + "/" + sampleid
        sample = Sample(sampleid=sampleid, fqs=fqs, outdir=outdir)
        check_dir(outdir)
        out_merged = sample.merge_gz()
        sample.get_fastqc(input=out_merged)
        [out_adapter,out_adapter_ungz]= sample.cut_adaptor(input=out_merged)
        out_rmdup = sample.rmdup(input=out_adapter_ungz)  # 4. Remove duplicate reads:
        out_sam = sample.bowtie1_mapping(input=out_adapter, bt_ref="/work/wllab/wentao/ref/hg38_female/female_bt1.hg38") # 5. bowtie1 mapping
        out_bam = sample.sam2bam(input=out_sam)
        out_bed = sample.bamtobed(input=out_bam)
        out_slop = sample.slop_bed(input=out_bed,ref_size="/work/wllab/wentao/afb1/hg38.chrom.sizes",)
        out_fixed = sample.fixed_slop(input=out_slop)
        out_upper = sample.low2uppercase(input=out_fixed)
        out_dis = sample.nucleotide_dis(input=out_upper)
        out_abun = sample.nucleotide_abundance(input=out_dis,sequence_length=None,nucleotide_order='ATGC',percentageFlag=True)


if __name__=='__main__':
    sample_list = sys.argv[1]
    outdir = sys.argv[2]
    runner(sample_list,outdir)
