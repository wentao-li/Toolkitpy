import os
from ..utils.cluster_module import check_module


class Sample(object):
    def __init__(self, id, fqs, outdir,):
        self._id = id
        self._fqs = fqs
        self._outdir = outdir
        self._commands = []

    @property
    def id(self):
        return self.id

    @property
    def fqs(self):
        return self._fqs

    @property
    def outdir(self):
        return self._outdir

    @property
    def commands(self):
        return self._commands

    def run_commands(self, command, output, step, module_load, ):
        check_module(module_load)
        os.system(command)
        self.commands.append(command)
        self.process.append(step)
        self.result[step] = output

    def check_gz(self,):
        fq_list = self._fqs.split(';')
        output = self._outdir + "/" + self._id + "_merged.fq.gz"
        if len(fq_list) > 1:
            command = "cat %s > %s" %(" ".join(fq_list), output)
            self.run_commands(command)
        return output

    def get_fastqc(self, input, module_load="module load FastQC/0.11.9-Java-11", ):
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
        output = self._outdir + "/" + self._id + "_cutadapt.fq.gz"
        command = f"cutadapt {param} -o {output} {input}\n"
        self.run_commands(
            command=command,
            output=output,
            step=type(self).cut_adaptor.__name__,
            module_load=module_load)

    def rmdup(
            self,
            module_load="module load FASTX-Toolkit/0.0.14-GCCcore-8.3.0",
            param="-v -Q33",
            input=None):
        output = self._outdir + "/" + self._id + "_rmdup.fq.gz"
        command = f"fastx_collapser -i {input} -o {output} {param}\n"
        self.run_commands(
            command=command,
            output=output,
            step=type(self).rmdup.__name__,
            module_load=module_load)

    def bowtie1_mapping(
            self,
            input,
            bt_ref, ## "/work/wllab/wentao/ref/hg38_female/female_bt1.hg38"
            module_load="module load Bowtie",
            param="--nomaqround --phred33-quals --chunkmbs 2000 -m 4 -n 2 -e 70 -l 20 --strata --best -p 4 --seed=123"):
        output = self._outdir + "/" + self._id + "_bt.sam"
        command = f"bowtie {param} -x {bt_ref} -q {input} -S {output}"
        self.run_commands(
            command=command,
            output=output,
            step=type(self).bowtie1_mapping.__name__,
            module_load=module_load)

    def sam2bam(
            self,
            input,
            module_load="module load SAMtools",
            param="-b -S"):
        output = self._outdir + "/" + self._id + "_bt.bam"
        command = f"samtools view {param} -o {output} {input}"
        self.run_commands(
            command=command,
            output=output,
            step=type(self).sam2bam.__name__,
            module_load=module_load)

    def bamtobed(
            self,
            input,
            param,
            module_load="module load  BEDTools",):
        output = self._outdir + "/" + self._id + "_bam2bed.sorted.bed"
        command = f"bedtools bamtobed -i {input} | sort -k1,1 -k2,2n > {output}"
        self.run_commands(
            command=command,
            output=output,
            step=type(self).bamtobed.__name__,
            module_load=module_load)

    def slop_bed(
            self,
            input,
            module_load="module load  BEDTools",
            param="-b 6 -s"):
        output = self._outdir + "/" + self._id + "_slop.sorted.bed"
        ref_size="/work/wllab/wentao/afb1/hg38.chrom.sizes"
        command = f"bedtools slop -i {input} -g {ref_size} {param} > {output}"
        self.run_commands(
            command,
            output=output,
            step=type(self).slop_bed.__name__,
            module_load=module_load)

    def fixed_slop(self, input, side="l", length=12,):
        output = self._outdir + "/" + self._id + f"_slop_{length}.sorted.bed"
        program = "bedline2fixedRangedLine.py"
        command = f"python {program} {input} {side} {length} tmp.bed\nsort -k1,1 -k2,2n tmp.bed > {output}\n"
        self.run_commands(
            command,
            output=output,
            step=type(self).fixed_slop.__name__,
            module_load=None)

    def getfasta(
            self,
            input,
            ref, # "/work/wllab/yiran/01.scripts/Damage_seq/refined.hg38.female.fa"
            param,
            module_load="module load  BEDTools"):
        output = self._outdir + "/" + self._id + f".fa"
        command = f"bedtools getfasta -fi {ref} -bed {input} -fo {output}\n"
        command = check_module(module_load, command)
        self.run_commands(
            command,
            output=output,
            step=type(self).getfasta.__name__,
            module_load=module_load)

    def low2uppercase(self, input, param):
        output = self._outdir + "/" + self._id + f"_uppercase.fa"
        command = f"cat {input} |tr 'a-z' 'A-Z' > {output}\n"
        self.run_commands(
            command,
            output=output,
            step=type(self).low2uppercase.__name__,
            module_load=None)

    def nucleotide_dis(self, input, param):
        #program = "/work/wllab/yiran/01.scripts/Damage_seq/Nucleotide_Dis2_V1.py"
        from ..utils.nucleotide_helper import nucleotide_dis2
        nucleotide_dis2(input=input, sampleid=self.id)
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
            module_load="module load Python/3.8.2-GCCcore-8.3.0",
            param="-n GACT --percentage"):
        program = "/work/wllab/yiran/01.scripts/Damage_seq/fastatonucleotideAbundanceTable.py"
        output = self._outdir + "/" + self._id + f"_getNuAbTa.csv"
        command = f"python {program} -i {input} -o {output} {param}\n"
        self.run_commands(
            command,
            output=output,
            step=type(self).nucleotide_abundance.__name__,
            module_load=None)

    def plot_nucleotide_abundance(self,module_load="module load R"):
        command = "Rscript {program} {input} toplot\n"
        command = check_module(module_load, command)
        self.run_commands(
            command,
            output=output,
            step=type(self).plot_nucleotide_abundance.__name__,
            module_load=None)

    def fa2kmerAbundanceTable(
            self,
            input,
            kmer=2,
            percentage=True,
            first_n_letters=None,):
        from .fasta import Fasta
        fasta = Fasta.fasta(input)
        output = f"{self.id}_get2kmer.csv"
        kmerAbundanceDict = fasta.getKmerAbundance(int(kmer), firstNletters=first_n_letters)
        fasta.writeKmerAbundanceTable(kmerAbundanceDict, output, percentage)
        '''
        program = "/work/wllab/yiran/01.scripts/Damage_seq/fa2kmerAbundanceTable.py"
        output = self._outdir + "/" + self._id + f"_get2kmer.csv"
        command = f"python {program} -i {input} -o {output} {param}\n"
        command = check_module(module_load, command)
        self.run_commands(
            command,
            output=output,
            step=type(self).fa2kmerAbundanceTable.__name__,
            module_load=None)
        '''

    def writer(self):
        output = self._outdir + "/" + self._id + f".sh"
        with open(output, "w") as f:
            f.write("\n".join(self.commands))
        f.close()
