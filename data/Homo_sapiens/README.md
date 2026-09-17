# demo set for transcriptome and proteome data integration

Place `WP534.gpml` file in `gpml` directory.

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
