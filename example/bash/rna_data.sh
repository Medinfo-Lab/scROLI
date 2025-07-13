files=($(find trim_fastq/scRNA_seq -maxdepth 1 -name "*_R1.fastq.gz" -printf "%f\n"))

for file in "${files[@]}"
do
	#echo ${file}
	name=${file%_*}
	echo ${name}
	bowtie2 -p 10 -x ~/genomes/mm10_bowtie2/mm10 -1 trim_fastq/scRNA_seq/${name}_R1.fastq.gz -2 trim_fastq/scRNA_seq/${name}_R2.fastq.gz | samtools view -b -@ 10 | samtools sort -@ 10 -m 2G -T ${name} -o rna_bam/${name}.sorted.bam && samtools index -@ 10 rna_bam/${name}.sorted.bam
done

/work/tools/featureCounts/subread-2.0.1-Linux-x86_64/bin/featureCounts -T 10 -t exon -g gene_id -a /work/genomes/Gencode/GRCm38/gencode.vM25.annotation.gtf -o SRP151137_RNA_counts.tsv rna_bam/*.sorted.bam
