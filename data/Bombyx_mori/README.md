# カイコ　WikiPathways リフトオーバー作業  

対応表になかった場合の対応  
- Ensembl bioMartでEnsembl Gene IDをクエリとしてUniProt IDを取得(MANE Selectでフィルタ)  
- Uniprot IDをもとにUniProtKBでヒトの対応するアミノ酸配列を取得(SilkBaseでカイコの全タンパク質配列データ(Protein sequences.Gene models based on the genome assembly (Nov.2016))を取得)  
-  →カイコのリファレンスをTSAの[ICPK00000000.1](https://www.ncbi.nlm.nih.gov/Traces/wgs?val=ICPK01)に変更  
- (blastpで確認(クエリ：ヒト、データベース：カイコ))  
-  →tblastnに変更

実行コマンド
```
# makeblastdb -in Bomo_gene_models_prot.fa -dbtype prot -hash_index -parse_seqids
# blastp -query idmapping_2025_12_11.fasta -db Bomo_gene_models_prot.fa -evalue 1e-10 -num_threads 8 -outfmt 6 -out WP534blastp-out.tsv

実行コマンド（WP534のもの）：
# 1) gzipを展開
gunzip -c ICPK01.1.fsa_nt.gz > ICPK01.1.fsa_nt.fa
# 2) nucl DBを作る（dbtypeをnuclに）
makeblastdb -in ICPK01.1.fsa_nt.fa -dbtype nucl -hash_index -parse_seqids -out ICPK01_TSA_nucl
# 3) blastn実行
tblastn -query idmapping_2025_12_11.fasta \
  -db ICPK01_TSA_nucl \
  -evalue 1e-10 -num_threads 8 \
  -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore sframe" \
  -out WP534tblastn-out.tsv
# 4) 遺伝子ごとに最小e-value（０が複数あれば全て）を抽出
awk -F'\t' 'BEGIN{OFS="\t"; ec=12}
FNR==NR {
  q=$1; e=$ec+0
  if(!(q in min) || e < min[q]) min[q]=e
  next
}
{
  q=$1; e=$ec+0
  if(e == min[q]) print
}' WP534tblastn-out.tsv WP534tblastn-out.tsv > best_by_query.tsv
# 5) カイコの遺伝子IDだけを抽出
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
# 6) ヒトタンパク質配列ID（UniProtKB/Swiss-Prot ID）を抽出
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
```


WP534
対応表になしのblastp確認結果(qseqid	sseqid	stitle	pident	evalue)
```
LDHA	
LDHC	
LDHAL6B	
 →対応済みLDHBと同一
G6PC	KWMTBOMO00773		sp|P35575|G6PC1_HUMAN	dbj|ICPK01023394.1|	TSA: Bombyx mori mRNA, KWMTBOMO00773.mrna1, mRNA sequence	29.444	8.57E-35
PFKL	11hits
 →対応済みPFKM,PFKPと同一
ENO3	4hits
ENO2	4hits
 →対応済みENO1と同一
PKM1	
PKM2	
 →PKMに統一 3hits
 PKM	MSTRG.21962.8		sp|P14618|KPYM_HUMAN	dbj|ICPK01051677.1|	TSA: Bombyx mori mRNA, MSTRG.21962.8, mRNA sequence	62.948	0
 PKM	MSTRG.21962.7		sp|P14618|KPYM_HUMAN	dbj|ICPK01051676.1|	TSA: Bombyx mori mRNA, MSTRG.21962.7, mRNA sequence	62.948	0
 PKM	KWMTBOMO05295		sp|P14618|KPYM_HUMAN	dbj|ICPK01051673.1|	TSA: Bombyx mori mRNA, KWMTBOMO05295.mrna1, mRNA sequence	62.948	0
PGK2	
→対応済みPGK1と同一
ALDOC	
ALDOB	
 →対応済みALDOAと同一
PGAM2	
 →対応済みPGAM1と同一
PCK1	7hits
 MSTRG.10528.3		sp|P35558|PCKGC_HUMAN	dbj|ICPK01024919.1|	TSA: Bombyx mori mRNA, MSTRG.10528.3, mRNA sequence	63.621	0
 KWMTBOMO12094		sp|P35558|PCKGC_HUMAN	dbj|ICPK01024920.1|	TSA: Bombyx mori mRNA, KWMTBOMO12094.mrna1, mRNA sequence	63.621	0
 KWMTBOMO12093		sp|P35558|PCKGC_HUMAN	dbj|ICPK01024916.1|	TSA: Bombyx mori mRNA, KWMTBOMO12093.mrna1, mRNA sequence	62.292	0
 MSTRG.10528.2		sp|P35558|PCKGC_HUMAN	dbj|ICPK01024918.1|	TSA: Bombyx mori mRNA, MSTRG.10528.2, mRNA sequence	66.186	0
 MSTRG.10528.1		sp|P35558|PCKGC_HUMAN	dbj|ICPK01024917.1|	TSA: Bombyx mori mRNA, MSTRG.10528.1, mRNA sequence	66.186	0
 MSTRG.10526.2		sp|P35558|PCKGC_HUMAN	dbj|ICPK01024915.1|	TSA: Bombyx mori mRNA, MSTRG.10526.2, mRNA sequence	65.155	0
 MSTRG.10526.1		sp|P35558|PCKGC_HUMAN	dbj|ICPK01024914.1|	TSA: Bombyx mori mRNA, MSTRG.10526.1, mRNA sequence	67.268	0
FBP2	
 →対応済みFBP1と同一
HK1	
HK2	
HK3
→対応済みGCKと同一
PGI　IdentiferがEnzyme Nomenclatureのため参照不可
```

