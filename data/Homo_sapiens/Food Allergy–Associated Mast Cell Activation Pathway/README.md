# Human Food Allergy–Associated Mast Cell Activation Pathway

## 概要

ヒト食物アレルギーにおけるIgE依存的な肥満細胞活性化を中心とした、試験的なGPMLパスウェイです。

論文・既存の生物学的知識・データベース等を参考に作成しており、PathVisioでの可視化やQPXによる遺伝子発現データのマッピングを目的としています。

本パスウェイは、完成された学術的パスウェイではなく、**学習・可視化・データマッピングのための初期モデル**として位置づけています。

## QPXテストデータ

QPXでの動作確認用として、NCBI GEOから取得した以下のRNA-seq発現データを使用しています。

* `GSE189149_norm_counts_TPM_GRCh38.p13_NCBI.tsv`
* `GSE201154_norm_counts_TPM_GRCh38.p13_NCBI.tsv`

いずれも、**NCBIが生成した「Series RNA-seq normalized counts matrix」のTPMデータ**を使用しています。

### GSE189149

食物アレルギー患者および非アレルギー者のnCD4 T細胞を対象としたRNA-seqデータです。48サンプルが登録されています。

[NCBI GEO: GSE189149](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189149)

### GSE201154

ヒト食道上皮のbulk RNA-seqデータです。PDPN（podoplanin）とE-cadherinをもとに分離したbasal層とsuprabasal層を対象としており、8サンプルが登録されています。

[NCBI GEO: GSE201154](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201154)

これらのデータは、主に**QPXによる発現データのマッピングおよび動作確認を目的としたテストデータ**として使用しています。

## 注意点

使用しているTPMデータはNCBIによって生成された正規化データであり、本パスウェイにおける可視化・動作確認を目的としています。

そのため、これらのTSVを用いた結果を正式な差次的発現解析や研究間の定量比較として解釈することは想定していません。

## Identifier

遺伝子ノードには可能な限りヒトEntrez Gene IDを使用しています。一方、代謝物や一部のノードにはChEBI、HMDB、ChemSpider等のidentifierを使用しています。

現時点では一部のidentifierが未整理のため、今後確認・修正する予定です。

## 今後の改善

今後は、

* 文献に基づくinteractionの再検証
* 省略した経路の詳細化
* Gene ID・isoform等の整理
* identifierの再確認
* 既存のキュレーション済みpathway databaseとの比較
* 細胞種や実験条件の整理

などを進め、より再現性の高いパスウェイへ発展させる予定です。

## 位置づけ

本GPMLは、**「正しいパスウェイを完成させたもの」ではなく、「正しいパスウェイを理解・整理するために作成した試作モデル」**です。

学術的に完全なモデルを提供することを目的とせず、今後の学習・検証・改良のための試案として公開しています。

## 参考文献

* https://doi.org/10.18452/21632
* https://doi.org/10.1038/ni.f.216
* https://doi.org/10.1016/j.febslet.2010.08.006

## データソース

* NCBI GEO: [GSE189149](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189149)
* NCBI GEO: [GSE201154](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201154)
