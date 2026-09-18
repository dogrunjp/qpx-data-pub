# Human Food Allergy–Associated Mast Cell Activation Pathway

## 概要

ヒト食物アレルギーにおけるIgE依存的な肥満細胞活性化を中心とした、試験的なGPMLパスウェイです。

論文・既存の生物学的知識・データベース等を参考に作成しており、PathVisioでの可視化やQPXによる遺伝子発現データのマッピングを目的としています。

本パスウェイは、完成された学術的パスウェイではなく、**学習・可視化・データマッピングのための初期モデル**として位置づけています。

## QPXテストデータ

QPXでの動作確認用として、NCBI GEOから取得した以下の2つのRNA-seq発現データを使用しています。

* `GSE201154_norm_counts_TPM_GRCh38.p13_NCBI.tsv`
* `GSE189149_norm_counts_TPM_GRCh38.p13_NCBI.tsv`

いずれも **NCBIが生成した「Series RNA-seq normalized counts matrix」のTPMデータ**を使用しています。NCBIのRNA-seq count pipelineにより処理されたデータであり、各GEO研究の投稿者が提供したprocessed dataとは異なる場合があります。

GSE201154はヒト食道上皮のbulk RNA-seqデータ、GSE189149は食物アレルギー関連のRNA-seqデータです。これらは主として、GPMLへの遺伝子発現値のマッピングおよびQPXの動作確認を目的として使用しています。


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

また、付属のRNA-seqデータは研究解析用の最終データではなく、**QPXでの可視化・動作確認を目的としたテストデータ**として収録しています。
