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

## Work log (BioHackathon 2026, from 2026-09-14)

Contributor: Haruka Tanimoto (@haru-tanimoto)

This section records where each file in this directory came from and how it was changed.
The ChEBI candidate/recommendation tool is out of scope here.

### Source data

| File used | Source | Notes |
|---|---|---|
| Supplementary Data Set 1 (`ds1.xlsx`, sheet `TPM`): *Transcripts per million (TPM) normalized expression values of genes identified in S. renifolius* | Tanimoto et al., *Plant Physiol.*, 2024, supplementary data (https://doi.org/10.1093/plphys/kiae059) | 15,904 transcripts × 24 samples (Pre/Hot/Post × Florets/Pith) |
| Supplementary Data Set 7 (`ds7.xlsx`, sheet `MetaboAnalyst 5.0`): *Metabolite concentrations across various thermogenic stages in S. renifolius* | Same as above | 93 metabolites, formatted as a MetaboAnalyst 5.0 input table |
| Supplemental Figure S10A–C | Same as above | Pathway figures used as drawing references for the maps |
| `SrEurofins_1000.fasta_Araport11_unique.tsv` | BLAST of the *S. renifolius* transcripts against Araport11 proteins | One best hit per query (59,810 queries, all E-value ≤ 1e-5) |
| `IDmapping.xlsx` | Made in this work from Data Set 7 (see below) | Data Set 7 as a standard data frame, with HMDB / PubChem / KEGG / ChEBI IDs |

### Gene expression table (`skunk_cabbage/data/annotation_Hot_qpx.tsv`)

1. **2026-09-15** Saved sheet `TPM` of Data Set 1 as tab-separated text (`ds1.txt`). No values were changed.
2. **2026-09-15** Added an Arabidopsis gene ID (AGI code) to each transcript.
   - The Araport11 ID was assigned from the BLAST result (`add_araport_ids.py`): `Gene ID` of Data Set 1 was joined to `qseqid` of `SrEurofins_1000.fasta_Araport11_unique.tsv`, and `sseqid` was taken with the transcript version removed (`AT4G26910.1` → `AT4G26910`).
   - Result: 11,438 / 15,904 transcripts got an AGI code (7,724 distinct AGI codes; one AGI code can be shared by several transcripts). No transcript had more than one candidate. The other 4,466 rows have an empty ID.
3. **2026-09-15** Reshaped for QPX.
   - Kept `Gene ID`, `Description (SWISS PROT)`, `Description (Araport11)` unchanged.
   - Inserted the AGI code as column 3 with header `xref_id`.
   - Kept only the 8 Hot columns (`Hot_F1`–`Hot_F4`, `Hot_P1`–`Hot_P4`) and dropped the 16 Pre/Post columns. The values are unchanged from Data Set 1.
   - Row count and row order are the same as Data Set 1 (15,904 rows); rows without an AGI code were kept.

### Metabolome table (`skunk_cabbage/data/metabolome_Hot_qpx.tsv`)

1. **2026-09-15** Converted Data Set 7 into a standard data frame (`IDmapping.xlsx`, saved as `metabolome.txt`).
   - Data Set 7 is laid out for MetaboAnalyst 5.0: the first row holds sample names and the second row holds group labels (`Label`, `Hot_F`, ...). It was changed to one header row and one row per metabolite.
   - `Hit`, `HMDB`, `PubChem` and `KEGG` were obtained by submitting the compound names to the MetaboAnalyst 6.0 *Metabolite ID Conversion* tool.
   - ChEBI accessions were assigned by hand for all 93 compounds, choosing the accession more commonly used in UniProt (Rhea) reaction annotations. We did not rely on automatic conversion because:
     - MetaboAnalyst does not return ChEBI IDs;
     - MetaboAnalyst could not assign IDs to some compounds (`Hit` is empty for 2 compounds);
     - automatic conversion from HMDB and similar IDs does not give a unique ChEBI ID, and the ID it gives may not be the form commonly used in reaction annotations.
2. **2026-09-15** Reshaped for QPX.
   - `Metabolite` = the original Data Set 7 compound name (`Query`). The matched name (`Hit`) was not used; it differs from `Query` in 16 rows (e.g. `L-Glutamine` → `Glutamine`) and is empty for 2 rows.
   - Kept `HMDB` and `KEGG`; dropped `PubChem` and `Label`.
   - The ChEBI ID (`Uniprot-Chebi`) became column 3 `xref_id`. All 93 rows have a ChEBI ID.
   - Kept only the 8 Hot columns. The duplicated headers (`Hot_F` ×4, `Hot_P` ×4) were renamed to `Hot_F1`–`Hot_F4`, `Hot_P1`–`Hot_P4`, because QPX needs unique column names. The values are unchanged.
   - Row count and order are the same as Data Set 7 (93 rows).
3. **2026-09-16** Added the `CHEBI:` prefix to all `xref_id` values (`32966` → `CHEBI:32966`), to match the GPML maps, WikiPathways, the QPX sample files and the PathVisio metabolite database (`metabolites_20260102.bridge`).

No ChEBI ID was replaced with a different ID after step 1.

> **Open issue:** the ChEBI assignment in step 1 is manual curation. Making this step reproducible is one of the improvements planned for future work and is under discussion by the team.

### Maps (`skunk_cabbage/gpml/FigS10A–C_qpx.gpml`)

1. **2026-09-14** Drew the three maps by hand in PathVisio 3.3.0, using Supplemental Figure S10A–C as reference.
   - Gene nodes: `Database="TAIR"`, AGI codes.
   - Metabolite nodes: `Database="ChEBI"`.
   - Node counts: S10A 33 genes + 13 metabolites; S10B 7 + 8; S10C 19 + 15 (+ 1 pathway node).
2. **2026-09-15** Edited the maps for QPX (saved as `*_qpx.gpml`).
   - Changed AGI codes to upper case (`At3g60750` → `AT3G60750`), because QPX links a node to a table row only when the IDs match exactly, including case.
   - Set `Organism="Symplocarpus renifolius"` for S10B and S10C (they had `Arabidopsis thaliana`).
   - Fixed three metabolite nodes:
     - S10A `E4P`: added the missing ChEBI ID `16897`
     - S10A `PRPP`: `15378` → `58017`
     - S10B `UMP`: `58017` → `57865`, and removed a line break from the label
     - These errors were made during manual drawing and were found by LLM review (see *Use of an LLM* below).
   - S10A only: redrew the reaction edges. The 7 curved (`ConnectorType="Curved"`) `Arrow` edges between metabolites were replaced with straight `mim-conversion` edges, the catalysis edges from the enzymes were re-attached to the new edges, two enzyme nodes were moved slightly, and the canvas size changed.
     - Reason: QPX does not draw S-shaped curves; it handles straight lines and elbow connectors only (`qpx_widgets/src/widget.ts`). Straight and elbow lines are enough to lay out these maps. Limiting maps to a few simple elements also matches the QPX design idea: keep maps easy for computers and AI to handle, while keeping them easy for people to read.
3. **2026-09-16** Added the `CHEBI:` prefix to all metabolite IDs (`61548` → `CHEBI:61548`), for the same reason as the metabolome table.

No nodes were added or removed. S10B and S10C have no layout or edge changes. S10A edges were redrawn (see step 2 above).

### Use of an LLM

This work was done together with an LLM coding agent: Claude Code (model: Claude Opus 5, Claude Max plan). The agent was used to:

- review the hand-made maps and tables against the source data and the QPX rules. This review found the three wrong or missing ChEBI IDs in the maps (see *Maps*, step 2);
- read the QPX source code to find the file format rules (for example the `xref_id` column, 0-based `expression_columns_index`, exact-case ID matching, supported line types);
- write and check the scripts for joining the BLAST result, and reshape the tables;
- write the documentation, including this record.

All decisions about the data, such as which ChEBI ID to use, were made by the contributor.

### Display setup and publication

- **2026-09-15/16** Wrote `skunk_cabbage.ipynb` (both tables, `expression_columns_index=[4, 4]`) and checked in QPX (`bonohu/qpx` commit `83d27bf`) that the maps, node selection and heatmap filtering work.
- **2026-09-17** Published the minimal set to https://github.com/haru-tanimoto/skunk-cabbage-pathway-analysis and to this directory.
- Found that the Python side of QPX reads `xref_id` as an integer (see *Known limitations*).

## Known limitations

- The Python side of QPX (`visualizers.py`) reads `xref_id` as an integer. String IDs such as AGI codes and `CHEBI:`-prefixed IDs therefore become null, and `visualizer.selected_expression_data` has 0 rows. Map display, node selection and heatmap filtering are not affected, because the browser side handles them
- The base image is built for x86_64. On Apple Silicon it runs under emulation and is slow
- Port 8888 is used
- If another qpx container on the same machine causes a name conflict, set a different project name with `docker compose -p <any name> ...`
