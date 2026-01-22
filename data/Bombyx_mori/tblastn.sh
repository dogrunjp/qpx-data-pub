#!/bin/bash
# tblastn実行
tblastn -query idmapping_2026_01_21.fasta \
  -db ICPK01_TSA_nucl \
  -evalue 1e-10 -num_threads 8 \
  -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore sframe" \
  -out tblastn-out.tsv
# 遺伝子ごとに最小e-value（０が複数あれば全て）を抽出
awk -F'\t' 'BEGIN{OFS="\t"; ec=12}
FNR==NR {
  q=$1; e=$ec+0
  if(!(q in min) || e < min[q]) min[q]=e
  next
}
{
  q=$1; e=$ec+0
  if(e == min[q]) print
}' tblastn-out.tsv tblastn-out.tsv > best_by_query.tsv
# カイコの遺伝子IDだけを抽出
awk -F'\t' 'BEGIN{
  OFS="\t";
  print "qseqid","sseqid","stitle","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","sframe","B_mori_gene"
}
{
  st=$3
  gene="NA"

  # KWMTBOMO05295.mrna1 など → KWMTBOMO05295
  if (match(st, /KWMTBOMO[0-9]+/)) {
    gene = substr(st, RSTART, RLENGTH)

  # MSTRG.21962.8 → そのまま
  } else if (match(st, /MSTRG\.[0-9]+\.[0-9]+/)) {
    gene = substr(st, RSTART, RLENGTH)
  }

  print $0, gene
}' best_by_query.tsv > best_by_query_withBmGene.tsv
# ヒトタンパク質配列ID（UniProtKB/Swiss-Prot ID）を抽出
awk -F'\t' 'BEGIN{OFS="\t"}
NR==1{
  print $0, "UniProtKB/Swiss-Prot ID"
  next
}
{
  q=$1
  id="NA"

  # 例: sp|P14618|KPYM_HUMAN → P14618
  n = split(q, a, /\|/)
  if(n >= 3) id = a[2]

  print $0, id
}' best_by_query_withBmGene.tsv > best_by_query_withBmGene_withUniProt.tsv