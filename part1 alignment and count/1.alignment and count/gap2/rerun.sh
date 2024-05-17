#! /usr/bin/sh
output_dir="/panfs/compbio/users/zzhu244/volumes/GAP2/qc_output/outcome"
mkdir -p $output_dir
fastp_output_dir="/panfs/compbio/users/zzhu244/volumes/GAP2/qc_output/QC"
mkdir -p $fastp_output_dir

while read -r filename; do
  STAR --readFilesIn "${fastp_output_dir}/${filename}_clean.fastq" --outSAMtype BAM SortedByCoordinate --runThreadN 15 --genomeDir /projects/compbio/users/sshar28/refData/ucsc/hg38/star_index_ucsc_refgene_hg38/ --outFileNamePrefix "${output_dir}/${filename}" --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --quantMode GeneCounts --twopassMode Basic
done < fail_sample.txt
