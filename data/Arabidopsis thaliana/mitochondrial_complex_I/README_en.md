# GPML pathway diagram of Arabidopsis mitochondrial complex I

[日本語](README.md) | English

This dataset contains a pathway diagram of mitochondrial complex I (NADH:ubiquinone oxidoreductase) in *Arabidopsis thaliana*. [At_mitochondrial_complex_I.gpml](At_mitochondrial_complex_I.gpml) records component names and identifiers, layout, lines representing reactions, groups, and literature references.

**Currently, only the GPML pathway diagram is included.** No transcriptome or proteome data, gene expression tables, or protein abundance tables are provided. The current use is to display the pathway diagram and inspect its annotations.

## Included files

| File | Description |
| --- | --- |
| [At_mitochondrial_complex_I.gpml](At_mitochondrial_complex_I.gpml) | Mitochondrial complex I pathway diagram in GPML 2013a format |
| [README.md](README.md) | Japanese documentation: information extracted from GPML, viewing instructions, and references |
| [README_en.md](README_en.md) | English documentation |

## GPML metadata

The following values are recorded in the GPML file.

| Field | Value |
| --- | --- |
| Pathway name (`Name`) | `mitochondrial complex I` |
| Organism (`Organism`) | `Arabidopsis thaliana` |
| Version attribute (`Version`) | `20260918` |
| GPML namespace | `http://pathvisio.org/GPML/2013a` |
| Drawing area | Approximately 2383.33 × 1228 in GPML coordinates |
| Pathway-level literature references | 3 (`b23`, `ec5`, and `e66`) |

`Version` is reported exactly as stored in the file. This GPML does not record a WikiPathways pathway ID (WP number), an author, or a distribution URL.

## Information in the pathway diagram

### Diagram elements

| GPML element | Count | Description |
| --- | ---: | --- |
| `DataNode` / `GeneProduct` | 61 | Nodes recorded as gene products |
| `DataNode` / `Protein` | 1 | Node recorded as a protein |
| `DataNode` / `Metabolite` | 5 | NADH, NAD+, and H+ (H+ appears in three locations) |
| `DataNode` / `Complex` | 6 | Headings for modules, domains, and candidate subunits |
| `Interaction` | 3 | Lines representing conversion, catalysis, and H+ movement |
| `Shape` | 5 | Shapes enclosing diagram regions |
| `Group` | 8 | Groups of nodes |

There are 73 `DataNode` elements in total. These counts describe diagram nodes and do not represent the number of confirmed complex I subunits.

### Modules and domains

The GPML contains the following six headings:

- `N module`
- `Q module`
- `PP module`
- `PD module`
- `CA domain`
- `Candidate plant-specific subunits`

Literature on plant complex I composition and assembly, the domain containing carbonic anhydrases, and internal architecture is listed in the References section below.

### Selected gene and protein nodes

The following display names and cross-references are extracted from the GPML. Names follow the annotations in the file.

| Display name | Description in GPML | `Xref Database` | `Xref ID` |
| --- | --- | --- | --- |
| CI51 | 51 kDa subunit of complex I | Ensembl Plants | `AT5G08530` |
| EMB1467 | NADH-ubiquinone dehydrogenase | Ensembl | `AT5G37510` |
| FRO1 | NADH-ubiquinone oxidoreductase-like protein | Ensembl | `AT5G67590` |
| nad7 | NADH dehydrogenase subunit 7 | Ensembl | `ATMG00510` |
| nad9 | NADH dehydrogenase subunit 9 | Ensembl Plants | `ATMG00070` |
| nad4 | NADH dehydrogenase subunit 4 | Ensembl | `ATMG00580` |
| nad4L | NADH dehydrogenase subunit 4L | Ensembl | `ATMG00650` |
| nad6 | NADH dehydrogenase subunit 6 | Ensembl | `ATMG00270` |
| GAMMA CA1 | gamma carbonic anhydrase 1 | Ensembl | `AT1G19580` |
| GAMMA CA2 | gamma carbonic anhydrase 2 | Ensembl | `AT1G47260` |
| GAMMA CA3 | gamma carbonic anhydrase 3 | Ensembl | `AT5G66510` |
| GAMMA CAL1 | gamma carbonic anhydrase like 1 | Ensembl | `AT5G63510` |
| GAMMA CAL2 | gamma carbonic anhydrase-like 2 | Ensembl | `AT3G48680` |
| AT2G42310 | ESSS subunit of NADH:ubiquinone oxidoreductase (complex I) protein | Uniprot-TrEMBL | `A0A178VYI3` |
| TIM23-2 | translocase inner membrane subunit 23-2 | Ensembl | `AT1G72750` |

### Metabolites and reaction lines

| Metabolite label | `Xref Database` | `Xref ID` | Node count |
| --- | --- | --- | ---: |
| NADH | ChEBI | `CHEBI:16908` | 1 |
| Nadide (NAD+) | ChEBI | `CHEBI:15846` | 1 |
| hydron (H+) | ChEBI | `CHEBI:15378` | 3 |

The three `Interaction` elements describe conversion from NADH to NAD+ (`mim-conversion`), a catalytic connection from the shape enclosing the N module to the conversion line (`mim-catalysis`), and an arrow between two H+ nodes. The H+ arrow is dashed. The GPML does not include ubiquinone/ubiquinol nodes or the stoichiometry of the complete reaction.

### Identifier considerations

There are 66 `DataNode` elements with cross-references: 54 use `Ensembl`, 6 use `Ensembl Plants`, 1 uses `Uniprot-TrEMBL`, and 5 use `ChEBI`. Seven nodes have empty cross-references: the six `Complex` headings and one nad2 node.

When mapping omics data in the future, check the following points:

- The nad2 label includes `AtMg00285/AtMg01320`, but both `Xref Database` and `Xref ID` are empty.
- The nad3 label is `AT2G07751 (nad3)`, while its `Xref ID` is `ATMG00990`. Distinguish the display label from the cross-reference ID when mapping data.
- The protein node labeled `AT2G42310` uses the UniProt accession `A0A178VYI3` as its cross-reference rather than an AGI code.
- Multiple loci for genes such as nad1 and nad5 are recorded as separate nodes. Decide whether to consolidate these nodes into a single gene according to the identifiers in the data being used.

## Viewing the pathway diagram

### Viewing with QPX (configuration example for `qpx.ipynb`)

1. Prepare a Jupyter environment with [QPX](https://github.com/bonohu/qpx) available.
2. Copy `At_mitochondrial_complex_I.gpml` into the `gpml/` directory of the QPX project.
3. Replace the visualization cell in `qpx.ipynb` with the configuration below and run it.
4. Select `At_mitochondrial_complex_I.gpml` in the GPML file selector.
5. Check labels such as `CI51` and `GAMMA CA1`, the module layout, and the line from NADH to NAD+.

```python
import qpx_widgets

# Use an empty list of expression data paths because no omics data are available.
visualizer = qpx_widgets.GpmlD3Visualizer(
    expression_data_path=[],
    expression_columns_index=[],
)
visualizer.show()
```

This example is based on the [current QPX implementation](https://github.com/bonohu/qpx/blob/main/qpx_widgets/qpx_widgets/visualizers.py), which constructs the pathway view without loading tables when the expression data path list is empty. Rendering in a QPX runtime environment has not been verified. Expression tables, heatmaps, HN-scores, and DAP flags are not part of this viewing check.

### Viewing with PathVisio

You can also view GPML pathway diagrams in [PathVisio](https://pathvisio.org/documentation/GPML). Start PathVisio and select `At_mitochondrial_complex_I.gpml` through `File → Open`. See the [official tutorial](https://pathvisio.org/tutorials/HandsOnTutorial.pdf) for instructions on opening a pathway file.

## Adding omics data in the future

Prepare separate transcriptome or proteome TSV files with a column corresponding to the GPML cross-reference IDs. Because AGI codes and UniProt accessions are both used, an identifier mapping table or identifier standardization may be needed. Specify whether mapping uses `TextLabel` or `Xref ID`, and check the key column and data types expected by QPX.

After adding TSV files, set `expression_data_path` to the actual file paths and configure `expression_columns_index` and `filter_key` according to the table structure. No TSV column structure or expression values are defined at this stage.

## References

The following three publications are recorded in the GPML's `Biopax` element and referenced by pathway-level `BiopaxRef` elements. Authors, titles, years, and PubMed IDs were extracted from the GPML; volume, issue, page, and DOI information was supplemented from PubMed records. These publications are cited as pathway references, not as sources of omics data included in this dataset.

1. **Braun HP et al. (2014).** The life of plant mitochondrial complex I. *Mitochondrion.* 19 Pt B:295–313. [DOI: 10.1016/j.mito.2014.02.006](https://doi.org/10.1016/j.mito.2014.02.006). [PubMed: 24561573](https://pubmed.ncbi.nlm.nih.gov/24561573/). GPML reference ID: `b23`.
   A review of plant mitochondrial complex I composition, transcript maturation, assembly, and the carbonic anhydrase domain.

2. **Subrahmanian N, Remacle C, Hamel PP (2016).** Plant mitochondrial Complex I composition and assembly: A review. *Biochim Biophys Acta.* 1857(7):1001–1014. [DOI: 10.1016/j.bbabio.2016.01.009](https://doi.org/10.1016/j.bbabio.2016.01.009). [PubMed: 26801215](https://pubmed.ncbi.nlm.nih.gov/26801215/). GPML reference ID: `ec5`.
   A review of plant complex I subunit composition, plant-specific components, assembly factors, and assembly processes.

3. **Klodmann J, Sunderhaus S, Nimtz M, Jänsch L, Braun HP (2010).** Internal architecture of mitochondrial complex I from Arabidopsis thaliana. *Plant Cell.* 22(3):797–810. [DOI: 10.1105/tpc.109.073726](https://doi.org/10.1105/tpc.109.073726). [PubMed: 20197505](https://pubmed.ncbi.nlm.nih.gov/20197505/).
   GPML reference ID: `e66`. A study of Arabidopsis complex I internal architecture and subunit arrangement based on subcomplex separation and mass spectrometry.
