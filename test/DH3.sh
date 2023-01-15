cat fq_A_1.gz  > .//DH3_merged.fq.gz
module load FastQC/0.11.9-Java-11
fastqc .//DH3_merged.fq.gz