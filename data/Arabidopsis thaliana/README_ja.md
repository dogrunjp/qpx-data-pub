## シロイヌナズナ 光合成炭素還元経路（WP1461）× QPX 対応データセット

シロイヌナズナ（*Arabidopsis thaliana*）の地上部器官発生RNA-seqデータ（[ArrayExpress E-MTAB-7978](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-7978)）と、WikiPathwaysの光合成炭素還元経路図（[WP1461](https://sandbox.wikipathways.org/pathways/WP1461.html)）をもとに、さらに文献情報を用いて編集した経路図です

この編集した経路データ上に、[QPX](https://github.com/bonohu/qpx)を用いてRNA-seqデータを表示できるようにした一式です。

## Pathway制作に使用した引用文献

* Evans SE et al. Rubisco supplies pyruvate for the 2-C-methyl-D-erythritol-4-phosphate pathway. *Nat Plants.* 2024;10:1445–1455. doi:10.1038/s41477-024-01791-z.

* Wang R et al. Structural insights into the functions of Raf1 and Bsd2 in hexadecameric Rubisco assembly. *Mol Plant.* 2023;16(12):1927–1936. doi:10.1016/j.molp.2023.10.011.

* Hanke GT, Hase T. Variable photosynthetic roles of two leaf-type ferredoxins in Arabidopsis, as revealed by RNA interference. *Photochem Photobiol.* 2008;84:1302–1309. doi:10.1111/j.1751-1097.2008.00411.x.

* Lehtimäki N et al. Posttranslational modifications of FERREDOXIN-NADP+ OXIDOREDUCTASE in Arabidopsis chloroplasts. *Plant Physiol.* 2014;166(4):1764–1776. doi:10.1104/pp.114.249094.

* Tikkanen M et al. Electron flow from PSII to PSI under high light is controlled by PGR5 but not by PSBS. *Front Plant Sci.* 2015;6:521. doi:10.3389/fpls.2015.00521.


## その他変更点

1. **GPMLのXrefを全部「接尾辞なしのプレーンなAGIロケルス」に統一**
   - 入手した発現データが全てAGIコードで紐つけられていたため、Grameneの`-TAIR-G`を除去し、Entrez番号3件をNCBI Geneで調べてAGIに変換(`821114`→`AT3G01550`=ATPPT2、`833635`→`AT5G36700`=ATPGLP1、`819611`→`AT3G04550`=Raf1)
   その後、Database="TAIR"に統一
   

2. **発現データ(TSV・TXT)の編集**
   - `transcript_id`はRefSeqに個別レコードが存在しないため空欄のまま。`xref_id`はATPPT2/ATPGLP1/Raf1の3件のみ判明した番号を記入、残りは空欄

## 成果物

- `RuBisco_Arabidopsis_thaliana.gpml` ── 編集した経路図
- `RuBisco_Arabidopsis_thaliana.tsv` ── 編集した発現データ(TSV形式)
- `RuBisco_Arabidopsis_thaliana.txt` ── 同上(TXT形式)

## 使い方

コードを以下の形に変更し、使用する

visualizer = qpx_widgets.GpmlD3Visualizer(
    expression_data_path=expression_data_path,
    expression_columns_index=[4, 4],
    filter_key="gene",
)
visualizer.show()
```
発現データの用語の意味については`E-MTAB-7978-experiment-design.tsv`を参照
