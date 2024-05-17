#! /usr/bin/sh
output_dir="/panfs/compbio/users/zzhu244/volumes/GAP3/qc_output/outcome"
mkdir -p $output_dir
input_dir="/panfs/compbio/users/zzhu244/volumes/GAP3/qc_output"
fastp_output_dir="/panfs/compbio/users/zzhu244/volumes/GAP3/qc_output/QC"
mkdir -p $fastp_output_dir


for fq in $input_dir/*.fastq
do
  output_prefix=$(basename -s .fastq -- "$fq")
  fastp -i "$fq" -o "$fastp_output_dir/${output_prefix}_clean.fastq"
done

for fq in $fastp_output_dir/*.fastq
do
  output_prefix=$(basename -s _clean.fastq -- "$fq")
  STAR --readFilesIn "$fq" --outSAMtype BAM SortedByCoordinate --runThreadN 15 --genomeDir /projects/compbio/users/sshar28/refData/ucsc/hg38/star_index_ucsc_refgene_hg38/ --outFileNamePrefix "${output_dir}/${output_prefix}" --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --quantMode GeneCounts --twopassMode Basic
done
