# Skunk cabbage (*Symplocarpus renifolius*) data for QPX

A minimal set for reproducing the display of gene expression and metabolome data from the thermogenic tissues Hot_F / Hot_P of the Asian skunk cabbage (*Symplocarpus renifolius*) on pathway maps with [QPX](https://github.com/bonohu/qpx).

## Contents

```
data/Symplocarpus_renifolius/
├── README.md
├── skunk_cabbage.ipynb              Notebook for display
└── skunk_cabbage/
    ├── gpml/                        Maps (drawn with PathVisio)
    │   ├── FigS10A_qpx.gpml         Pentose phosphate pathway
    │   ├── FigS10B_qpx.gpml         Pyrimidine biosynthesis
    │   └── FigS10C_qpx.gpml         Purine and histidine biosynthesis
    └── data/
        ├── annotation_Hot_qpx.tsv   Gene expression (15,904 rows; xref_id = Araport11 AGI code)
        └── metabolome_Hot_qpx.tsv   Metabolome (93 rows; xref_id = CHEBI:xxxxx)
```

Source data: supplementary data of Tanimoto et al., *Plant Physiol.*, 2024 (https://doi.org/10.1093/plphys/kiae059)

**How the maps were made:** The three maps (FigS10A–C) were digitized and built manually in PathVisio 3.3.0, using the pathway figures in Tanimoto et al., *Plant Physiol.*, 2024 as reference. Automation of this step is currently in progress.

## How to reproduce

Tested with `bonohu/qpx` commit `83d27bf`.

```bash
# 1. Get QPX and check out the tested version
git clone https://github.com/bonohu/qpx.git
cd qpx
git checkout 83d27bf

# 2. Copy the two items in qpx-data-pub/data/Symplocarpus_renifolius/ to the top level of qpx/
#    (paths in the notebook are relative to qpx/)
cp -r /path/to/qpx-data-pub/data/Symplocarpus_renifolius/skunk_cabbage .
cp /path/to/qpx-data-pub/data/Symplocarpus_renifolius/skunk_cabbage.ipynb .

# 3. Build and start (the first build takes a while)
docker compose build
docker compose up -d
```

1. Open http://localhost:8888/lab/tree/skunk_cabbage.ipynb in a browser (no token required)
2. Run the cells from top to bottom
3. Select a map from the dropdown. Clicking a node filters the heatmaps to that gene or metabolite

To stop: `docker compose down -v`

## File format

- In the maps, gene nodes use `Database="TAIR"` with upper-case AGI codes (`AT3G60750`). Metabolite nodes use `Database="ChEBI"` with IDs that include the `CHEBI:` prefix (`CHEBI:61548`)
- The TSV files are tab-separated and have a header row. Counting from 0, column 3 is `xref_id` and columns 4 onward hold the expression values (Hot_F1–F4, Hot_P1–P4). `expression_columns_index=[4, 4]` in the notebook refers to this layout
- A map `ID` is linked to a TSV `xref_id` only when the two match exactly, including case

Map nodes that have a matching row in the TSV files:

| Map | Genes (TAIR) | Metabolites (ChEBI) |
|---|---|---|
| FigS10A | 19 / 33 | 4 / 13 |
| FigS10B | 7 / 7 | 2 / 8 |
| FigS10C | 13 / 19 | 1 / 15 |

Nodes without a matching row are still drawn. Clicking one just shows no rows in the heatmaps.

## Known limitations

- The Python side of QPX (`visualizers.py`) reads `xref_id` as an integer. String IDs such as AGI codes and `CHEBI:`-prefixed IDs therefore become null, and `visualizer.selected_expression_data` has 0 rows. Map display, node selection and heatmap filtering are not affected, because the browser side handles them
- The base image is built for x86_64. On Apple Silicon it runs under emulation and is slow
- Port 8888 is used
- If another qpx container on the same machine causes a name conflict, set a different project name with `docker compose -p <any name> ...`
