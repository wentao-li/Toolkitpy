#from ExcisionRepairToolkit.core.sample import Sample
#from ExcisionRepairToolkit.core.bedline import Bedline
from ExcisionRepairToolkit.core.bed import Bed
from ExcisionRepairToolkit.core.pipeline import XRseqPipeline
from ExcisionRepairToolkit.plots.ts_nts import plot_ts_nts


def runner(sample_bed, gene_list, sampleid):
    '''
    sl = open(file, "r")
    for lines in sl:
        [sample_id, fqs] = lines.strip().split('\t')
        outdir = "./"
        sample = Sample(sample_id, fqs, outdir)
        output = sample.check_gz()
        sample.get_fastqc(input=output)
        sample.writer()
    '''
    # getreadcount
    sample_bed = Bed(sample_bed)
    out_readcount = sample_bed.get_readcount(sampleid=sampleid)

    gene_list = Bed(gene_list)
    # expand region to 6kb and split into 150 bins.
    out_divide = gene_list.divide_transcript(sampleid=sampleid)
    # rm intersection
    out_intersection = Bed(out_divide).remove_intersection(sampleid=sampleid, gap=6000)
    #
    pipe = XRseqPipeline(sample_bed, gene_list, sampleid)
    #
    out_filter = pipe.cut_bylength(out_intersection, sampleid=sampleid)
    out_ts = pipe.get_TSandNTS(sample_bed=sample_bed, gene_list=out_filter, output=f"{sampleid}_TS.bed")
    out_nts = pipe.get_TSandNTS(sample_bed=sample_bed, gene_list=out_filter, output=f"{sampleid}_NTS.bed", param="-wa -c -s -F 0.5")
    out_rpkm = pipe.get_rpkm(out_readcount, out_ts, out_nts, sampleid=sampleid)
    out_tsmean = pipe.get_rpkm_mean(out_rpkm,output=f"{sampleid}_ts_mean.txt")
    out_ntsmean = pipe.get_rpkm_mean(out_rpkm, output=f"{sampleid}_nts_mean.txt")
    plot_ts_nts(ts_file=out_tsmean,nts_file=out_ntsmean,output=f"{sampleid}_ts_nts.png")


if __name__ == "__main__":
    runner("rep1.bed", "./test.bed", "test")
