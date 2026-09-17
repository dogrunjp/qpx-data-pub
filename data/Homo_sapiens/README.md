# Test set for transcriptome and proteome data integration

- Use case for the integration of multiple omics data 
  - Transcriptome: hypoxic stress RNA-seq data (HN-score) collected in "Multi-Omic Meta-Analysis of Transcriptomes and the Bibliome Uncovers Novel Hypoxia-Inducible Genes. [DOI: 10.3390/biomedicines9050582](https://doi.org/10.3390/biomedicines9050582)"
  - Proteome: differentially abundant proteins (DAPs) flags in "Proteomic-Based Analysis of Hypoxia- and Physioxia-Responsive Proteins and Pathways in Diffuse Large B-Cell Lymphoma. [DOI: 10.3390/cells10082025](https://doi.org/10.3390/cells10082025)"

1. Place `WP534.gpml` file (Glycolysis and gluconeogenesis [`WP534`](https://www.wikipathways.org/pathways/WP534.html)) in `gpml` directory.
2. Modify cells (test code below).
3. Choose `WP534.gpml` for GPML file.
4. Click `ENO2` to verify that it is set to DAP and corresponding HN-score.

## test code for `qpx.ipynb`
```
# Specify the path of the TSV file in which the expression information is described.
expression_data_path = ['data/TP.tsv']
# Running this cell will visualise the pathways described in the GPML; you can also visualise different pathways by modifying the GPML file.
# The columns after the ‘expression_columns_index’ parameter of the GpmlD3Visualizer described in the cell are coloured in the table as a heatmap of expression levels.
import qpx_widgets
visualizer = qpx_widgets.GpmlD3Visualizer(expression_data_path=expression_data_path, expression_columns_index=[3])
visualizer.show()
```
