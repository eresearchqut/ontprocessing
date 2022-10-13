#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process MERGE {
  publishDir "${params.outdir}/merge", pattern: '*.fastq.gz', mode: 'link'
  tag "${sampleid}"
  label 'process_very_low'

  input:
    tuple val(sampleid), path(lanes), path(reference)
  output:
    tuple val(sampleid), path("${sampleid}.fastq.gz"), path(reference), emit: merged
  script:
  """
  cat ${lanes} > ${sampleid}.fastq.gz
  """

}

process FLYE {
  publishDir "${params.outdir}/flye", mode: 'link'
  tag "${sampleid}"
  label 'process_medium'

  container "quay.io/biocontainers/flye:2.9.1--py310h590eda1_0"

  input:
    tuple val(sampleid), path(sample), path(reference)
  output:
    path 'outdir/*'
    tuple val(sampleid), path('outdir/assembly.fasta'), path(reference), emit: assembly
  script:
  """
  flye  --out-dir outdir --threads ${task.cpus} --read-error 0.03 --nano-hq ${sample}
  """
}

process BLASTN {
  publishDir "${params.outdir}/blastn", mode: 'link'
  tag "${sampleid}"
  label 'process_low'

  container 'quay.io/biocontainers/blast:2.13.0--hf3cf87c_0'

  input:
    tuple val(sampleid), path(aligned_sample), path(reference)
  output:
    path "BLASTN_${reference}_vs_${sampleid}_assembly.txt"
  script:
  """
  blastn -query ${aligned_sample} -subject ${reference} -evalue 1e-5 -out blastn_${reference}_vs_${sampleid}_assembly.txt \
    -outfmt '6 qseqid sacc length pident mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovhsp qcovs'
  
  echo "qseqid sacc length pident mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovhsp qcovs" > header
  
  cat header blastn_${reference}_vs_${sampleid}_assembly.txt > BLASTN_${reference}_vs_${sampleid}_assembly.txt
  """
}

process MINIMAP2 {
  publishDir "${params.outdir}/minimap2", mode: 'link'
  tag "${sampleid}"
  label 'process_low'

  container 'quay.io/biocontainers/minimap2:2.24--h7132678_1'

  input:
    tuple val(sampleid), path(sample),path(reference)
  output:
    tuple val(sampleid), file("${sampleid}_aln.sam"), path(reference), emit: aligned_sample
  script:
  """
  minimap2 -a --MD ${reference} ${sample} > ${sampleid}_aln.sam
  """
}

process INFOSEQ {
  publishDir "${params.outdir}/infoseq", mode: 'link'
  tag "${sampleid}"
  label 'process_low'

  container "quay.io/biocontainers/emboss:6.6.0--h1b6f16a_5"

  input:
    tuple val(sampleid), path(sample), path(reference)
  output:
    tuple val(sampleid), path(sample), path("${reference}_list.txt"), emit: infoseq_ref
  script:
  """
  infoseq ${reference} -only -name -length | sed 1d > ${reference}_list.txt
  """

}

process SAMTOOLS {
  publishDir "${params.outdir}/samtools", mode: 'link'
  tag "${sampleid}"
  label 'process_low'

  container 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'

  input:
    tuple val(sampleid), path(sample), path(reference_list)
  output:
    tuple val(sampleid), path("${sampleid}_aln.sorted.bam"), path("${sampleid}_aln.sorted.bam.bai"), emit: sorted_sample
  script:
  """
  samtools view -bt ${reference_list} -o ${sampleid}_aln.bam ${sample}
  samtools sort -T /tmp/aln.sorted -o ${sampleid}_aln.sorted.bam ${sampleid}_aln.bam
  samtools index ${sampleid}_aln.sorted.bam
  """
}

process NANOQ {
  publishDir "${params.outdir}/nano-q", mode: 'link'
  tag "${sampleid}"
  label 'process_medium'
  
  container 'ghcr.io/eresearchqut/nano-q:1.0.0'

  input:
    tuple val(sampleid), path(sorted_sample), path(sorted_sample_index)
  output:
    path 'Results/*'

  script:
  """
  nano-q.py -b ${sorted_sample} -c 1 -l 9000 -nr 1 -q 5 -j 10
  """
}

workflow {
// if (params.samplesheet) {
//     Channel
//       .fromPath(params.samplesheet, checkIfExists: true)
//       .splitCsv(header:true)
//       .map{ row-> tuple(row.sampleid), file("row.sample_folder"), file(row.reference) }
//       .set{ ch_sample }
//   } else { exit 1, "Input samplesheet file not specified!" }
  if (params.samplesheet) {
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> tuple(row.sampleid), file(row.sample_files), file(row.reference) }
      .set{ ch_sample }
  } else { exit 1, "Input samplesheet file not specified!" }

  MERGE ( ch_sample )
  FLYE ( MERGE.out.merged )
  BLASTN ( FLYE.out.assembly )
  MINIMAP2 ( MERGE.out.merged )
  INFOSEQ ( MINIMAP2.out.aligned_sample )
  SAMTOOLS ( INFOSEQ.out.infoseq_ref )
  NANOQ ( SAMTOOLS.out.sorted_sample )
}

