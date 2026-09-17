# ザゼンソウ(Symplocarpus renifolius) QPX表示用データ

ザゼンソウ(*Symplocarpus renifolius*)の発熱組織 Hot_F / Hot_P の遺伝子発現・メタボロームを、
[QPX](https://github.com/bonohu/qpx) でパスウェイマップ上に表示するための最小限の再現セット。

## 内容

```
data/Symplocarpus_renifolius/
├── README.md
├── skunk_cabbage.ipynb              表示用ノートブック
└── skunk_cabbage/
    ├── gpml/                        マップ(PathVisioで作成)
    │   ├── FigS10A_qpx.gpml         ペントースリン酸経路
    │   ├── FigS10B_qpx.gpml         ピリミジン合成
    │   └── FigS10C_qpx.gpml         プリン・His合成
    └── data/
        ├── annotation_Hot_qpx.tsv   遺伝子発現(15,904行、xref_id = Araport11 AGIコード)
        └── metabolome_Hot_qpx.tsv   メタボローム(93行、xref_id = CHEBI:xxxxx)
```

元データ: Tanimoto et al., *Plant Physiol.*, 2024 の補足データ(https://doi.org/10.1093/plphys/kiae059)

## 再現手順

動作確認に使ったQPXのバージョン: `bonohu/qpx` の commit `83d27bf`。

```bash
# 1. QPX を取得し、動作確認したバージョンに合わせる
git clone https://github.com/bonohu/qpx.git
cd qpx
git checkout 83d27bf

# 2. qpx-data-pub の data/Symplocarpus_renifolius/ にある2つを qpx/ の直下にコピーする
#    (ノートブック内のパスは qpx/ からの相対パス)
cp -r /path/to/qpx-data-pub/data/Symplocarpus_renifolius/skunk_cabbage .
cp /path/to/qpx-data-pub/data/Symplocarpus_renifolius/skunk_cabbage.ipynb .

# 3. ビルドして起動(初回のビルドには時間がかかる)
docker compose build
docker compose up -d
```

1. ブラウザで http://localhost:8888/lab/tree/skunk_cabbage.ipynb を開く(トークン不要)
2. セルを上から順に実行する
3. ドロップダウンでマップを選び、ノードをクリックするとヒートマップがその遺伝子・代謝物に絞り込まれる

停止: `docker compose down -v`

## ファイルの書式

- マップの遺伝子ノードは `Database="TAIR"` + 大文字のAGIコード(`AT3G60750`)、代謝物ノードは `Database="ChEBI"` + `CHEBI:` 接頭辞付きID(`CHEBI:61548`)
- TSVはタブ区切り・ヘッダーあり。0始まりで3番目の列が `xref_id`、4番目以降の列が発現値(Hot_F1–F4, Hot_P1–P4)。ノートブックの `expression_columns_index=[4, 4]` はこの配置を指す
- マップの `ID` とTSVの `xref_id` は、大文字・小文字の違いも含めて完全に一致するものだけが紐付く

マップ上のノードのうち、TSVに対応する行があるもの:

| マップ | 遺伝子(TAIR) | 代謝物(ChEBI) |
|---|---|---|
| FigS10A | 19 / 33 | 4 / 13 |
| FigS10B | 7 / 7 | 2 / 8 |
| FigS10C | 13 / 19 | 1 / 15 |

TSVに行がないノードも表示はされる。クリックしてもヒートマップに行が出ないだけ。

## 既知の制約

- QPXのPython側(`visualizers.py`)は `xref_id` を整数として読み込む。このため、AGIコードや `CHEBI:` 付きIDのような文字列IDは null になり、`visualizer.selected_expression_data` が0行になる。一方、マップの表示・ノード選択・ヒートマップの絞り込みはブラウザ側で処理するので影響を受けない
- ベースイメージは x86_64 用。Apple Silicon ではエミュレーションで動くので遅い
- ポート 8888 を使う
- 同じPCに別の qpx のコンテナが既にあって名前がぶつかる場合は、`docker compose -p <任意の名前> ...` でプロジェクト名を変える
