files=($(find ../trim_fastq/scBS_seq -maxdepth 1 -name "*_R1.fastq.gz" -printf "%f\n"))

for file in "${files[@]}"
do
	name=${file%_*}
	mkdir bam/${name}
	echo ${name}
	bismark --genome ~/genomes/mm10_bismark/mm10 --bowtie2 --parallel 5 -p 4 --bam -1 ../trim_fastq/scBS_seq/${name}_R1.fastq.gz -2 ../trim_fastq/scBS_seq/${name}_R2.fastq.gz --output_dir bam/${name}
done

files=($(find bam -maxdepth 1 -printf "%f\n"))

for file in "${files[@]}"
do
	echo ${file}
	mkdir dedup_bam/${bam}
	deduplicate_bismark --bam bam/${file}/${file}_R1_bismark_bt2_pe.bam --output_dir dedup_bam/${file}
done

files=($(find dedup_bam -maxdepth 1 -printf "%f\n"))

for file in "${files[@]}"
do
	echo ${file}
	mkdir methylation_extractor/${file}
	bismark_methylation_extractor -p --cytosine_report --parallel 10 --gzip --genome_folder ~/genomes/mm10_bismark/mm10 dedup_bam/${file}/${file}_R1_bismark_bt2_pe.deduplicated.bam --output_dir methylation_extractor/${file}
done

files=($(find methylation_extractor/one -maxdepth 1 -printf "%f\n"))

for file in "${files[@]}"
do
	echo ${file}
	mkdir coverage_cytosine/${file}
	coverage2cytosine --nome-seq --gzip --genome_folder ~/genomes/mm10_bismark/mm10 --dir coverage_cytosine/${file} -o ${file} methylation_extractor/one/${file}/${file}_R1_bismark_bt2_pe.deduplicated.bismark.cov.gz
done