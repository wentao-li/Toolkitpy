import os
from ..utils.cluster_module import check_module
import sys

class Sample(object):
    def __init__(self, sampleid, fqs, outdir,):
        self._sampleid = sampleid
        self._fqs = fqs
        self._outdir = outdir
        self._commands = []
        self._process = []
        self._result = dict()

    @property
    def sampleid(self):
        return self._sampleid

    @property
    def fqs(self):
        return self._fqs

    @property
    def outdir(self):
        return self._outdir

    @property
    def commands(self):
        return self._commands

    @property
    def process(self):
        return self._process

    def run_commands(self, command, output, step, module_load, ):
        command = check_module(module_load, command)
        os.system(command)
        self.commands.append(command)
        self.process.append(step)
        print(step)
        #self.result[step] = output

    def merge_gz(self,sep=";"):
        fq_list = self._fqs.split(sep)
        output = self._outdir + "/" + self._sampleid + "_merged.fq.gz"
        if len(fq_list) > 1:
            command = "cat %s > %s" %(" ".join(fq_list), output)
            self.run_commands(
                command=command,
                output="zip and html",
                step=type(self).merge_gz.__name__,
                module_load=None)

        return output

    def get_fastqc(self, input, module_load="module load FastQC/0.11.9-Java-11\nmodule list", ):
        ## generate *.zip and *.html
        command= f"fastqc {input}"
        self.run_commands(
            command=command,
            output="zip and html",
            step=type(self).get_fastqc.__name__,
            module_load=module_load)

    def cut_adaptor(
            self,
            input,
            param="--discard-trimmed -g GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT",
            module_load="module load cutadapt"):
        output = self._outdir + "/" + self._sampleid + "_cutadapt.fq.gz"
        output2 = self._outdir + "/" + self._sampleid + "_cutadapt.fq"
        command = f"cutadapt {param} -o {output} {input}\nzcat {output} > {output2}\n"
        self.run_commands(
            command=command,
            output=output,
            step=type(self).cut_adaptor.__name__,
            module_load=module_load)
        return [output,output2]

    def rmdup(
            self,
            input,
            module_load="module load FASTX-Toolkit/0.0.14-GCCcore-8.3.0",
            param="-v -Q33",):
        output = self._outdir + "/" + self._sampleid + "_rmdup.fq.gz"
        command = f"fastx_collapser -i {input} -o {output} {param}\n"
        self.run_commands(
            command=command,
            output=output,
            step=type(self).rmdup.__name__,
            module_load=module_load)
        return output

    def bowtie1_mapping(
            self,
            input,
            bt_ref, ## "/work/wllab/wentao/ref/hg38_female/female_bt1.hg38"
            module_load="module load Bowtie",
            param="--nomaqround --phred33-quals --chunkmbs 2000 -m 4 -n 2 -e 70 -l 20 --strata --best -p 4 --seed=123"):
        output = self._outdir + "/" + self._sampleid + "_bt.sam"
        command = f"bowtie {param} -x {bt_ref} -q {input} -S {output}"
        self.run_commands(
            command=command,
            output=output,
            step=type(self).bowtie1_mapping.__name__,
            module_load=module_load)
        return output

    def sam2bam(
            self,
            input,
            module_load="module load SAMtools",
            param="-b -S"):
        output = self._outdir + "/" + self._sampleid + "_bt.bam"
        command = f"samtools view {param} -o {output} {input}"
        self.run_commands(
            command=command,
            output=output,
            step=type(self).sam2bam.__name__,
            module_load=module_load)
        return output

    def bamtobed(
            self,
            input,
            module_load="module load  BEDTools",):
        output = self._outdir + "/" + self._sampleid + "_bam2bed.sorted.bed"
        command = f"bedtools bamtobed -i {input} | sort -k1,1 -k2,2n > {output}"
        self.run_commands(
            command=command,
            output=output,
            step=type(self).bamtobed.__name__,
            module_load=module_load)
        return output

    def slop_bed(
            self,
            input,
            ref_size,
            module_load="module load  BEDTools",
            param="-b 6 -s"):
        output = self._outdir + "/" + self._sampleid + "_slop.sorted.bed"
        command = f"bedtools slop -i {input} -g {ref_size} {param} > {output}"
        self.run_commands(
            command,
            output=output,
            step=type(self).slop_bed.__name__,
            module_load=module_load)
        return output

    def fixed_slop(self, input, side="l", length=12, strandFlag=True):
        output = self._outdir + "/" + self._sampleid + f"_slop_{length}.sorted.bed"
        from ..io.reader import read_bed
        input_bed = read_bed(input)
        out = open("tmp.bed", "w")
        for ll in input_bed:
        #ll = bedLine.strip().split('\t')
            beg = ll.start()
            end = ll.end()
            strand = ll.strand()
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
            newline = ll.newline(start=str(newBeg), end=str(newEnd))
            out.write(newline + "\n")
            # ll[1] = str(newBeg)
            # ll[2] = str(newEnd)
            # return '\t'.join(ll)
        out.close()
        program = "bedline2fixedRangedLine.py"
        #command = f"python {program} {input} {side} {length} tmp.bed\n"
        command = f"sort -k1,1 -k2,2n tmp.bed > {output}"

        self.run_commands(
            command,
            output=output,
            step=type(self).fixed_slop.__name__,
            module_load=None)
        return output

    def getfasta(
            self,
            input,
            ref, # "/work/wllab/yiran/01.scripts/Damage_seq/refined.hg38.female.fa"
            param,
            module_load="module load  BEDTools"):
        output = self._outdir + "/" + self._sampleid + f".fa"
        command = f"bedtools getfasta -fi {ref} -bed {input} -fo {output}\n"
        command = check_module(module_load, command)
        self.run_commands(
            command,
            output=output,
            step=type(self).getfasta.__name__,
            module_load=module_load)
        return output

    def low2uppercase(self, input,):
        output = self._outdir + "/" + self._sampleid + f"_uppercase.fa"
        command = f"cat {input} |tr 'a-z' 'A-Z' > {output}\n"
        self.run_commands(
            command,
            output=output,
            step=type(self).low2uppercase.__name__,
            module_load=None)
        return output

    def nucleotide_dis(self, input):
        #program = "/work/wllab/yiran/01.scripts/Damage_seq/Nucleotide_Dis2_V1.py"
        from ..utils.nucleotide_helper import nucleotide_dis2
        output = f"{self._sampleid}_nt_dis.txt"
        nucleotide_dis2(input=input, sampleid=self._sampleid, output=output)
        return output
        '''
        command = f"python {program} {input}"
        self.run_commands(
            command,
            output=output,
            step=type(self).nucleotide_dis.__name__,
            module_load=None)
        '''

    def nucleotide_abundance(
            self,
            input,
            sequence_length=None,
            nucleotide_order='ATGC',
            percentageFlag=True,
            module_load="module load Python/3.8.2-GCCcore-8.3.0",
            ):
        from ..utils.nucleotide_helper import get_nucleotide_abundance
        program = "/work/wllab/yiran/01.scripts/Damage_seq/fastatonucleotideAbundanceTable.py"
        output = self._outdir + "/" + self._sampleid + f"_getNuAbTa.csv"
        get_nucleotide_abundance(input=input, output=output,sequence_length=sequence_length, nucleotide_order=nucleotide_order,percentageFlag=percentageFlag)
        return output
        '''
        self.run_commands(
            command,
            output=output,
            step=type(self).nucleotide_abundance.__name__,
            module_load=None)
        '''

    def plot_nucleotide_abundance(self,module_load="module load R"):
        command = "Rscript {program} {input} toplot\n"
        command = check_module(module_load, command)
        output = self._outdir + "/" + self._sampleid + f"_getNuAbTa.pdf"
        self.run_commands(
            command,
            output=output,
            step=type(self).plot_nucleotide_abundance.__name__,
            module_load=None)
        return output

    def fa2kmerAbundanceTable(
            self,
            input,
            kmer=2,
            percentage=True,
            first_n_letters=None,):
        from .fasta import Fasta
        fasta = Fasta.fasta(input)
        output = f"{self._sampleid}_get2kmer.csv"
        kmerAbundanceDict = fasta.getKmerAbundance(int(kmer), firstNletters=first_n_letters)
        fasta.writeKmerAbundanceTable(kmerAbundanceDict, output, percentage)
        return output
        '''
        program = "/work/wllab/yiran/01.scripts/Damage_seq/fa2kmerAbundanceTable.py"
        output = self._outdir + "/" + self._sampleid + f"_get2kmer.csv"
        command = f"python {program} -i {input} -o {output} {param}\n"
        command = check_module(module_load, command)
        self.run_commands(
            command,
            output=output,
            step=type(self).fa2kmerAbundanceTable.__name__,
            module_load=None)
        '''

    def writer(self):
        output = self._outdir + "/" + self._sampleid + f".sh"
        with open(output, "w") as f:
            f.write("\n".join(self.commands))
        f.close()
