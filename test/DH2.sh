cat fq_A_1.gz fq_B_1.gz fq_c_1.gz > .//DH2_merged.fq.gz
module load FastQC/0.11.9-Java-11
fastqc .//DH2_merged.fq.gz