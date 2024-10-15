Code for analyze the single cell transcriptome data of endometriosis

1. alldata2024.ipynb, is the scanpy workflow used to analyze our single cell sequencing data, and get the cell clusters and cell types.

2.1 readdata_E-MTAB-10287_NG.R.ipynb, is the R workflow to load the endometrium scRNA data of E-MTAB-10287, and convert to h5ad format.
2.2 Integration_ourEU_with_E-MTAB-10287_NG.ipynb, is the scanpy workflow to integrated analysis of our endometriuma scRNA data with E-MTAB-10287, and compare the endoemtrial stromal cell subtypes.

3. ourOMA_OvayControl_integration.ipynb, is the scanpy workflow to integrated analysis of our endometrioma scRNA data with public ovary scRNA datasets (GSE118127 and GSE260685) , and compare the ovary stromal cells.

4.0 Fonseca_rds_to_h5ad.R.ipynb, is the R workflow to load the Fonseca's dataset, and convert to h5ad format. Fonseca's data was well analyzed and downloaded from https://cedars.app.box.com/s/1ks3eyzlpnjbrseefw3j4k7nx6p2ut02
4.1 Fonseca_NGdata_scanpy_OMA.ipynb, is the scanpy workflow used to subset the endometrioma data from Fonseca's dataset and reanalyzed.
4.2 Fonseca_NGdata_PEM_EU_OV_scanpy.ipynb, is the workflow to analyze the Fonseca's data subset according to the derived tissue. Specifically, the peritoneal endometriosis, eutopic endometrium, and un-affected ovaries.
4.3 Fonseca_NGdata_4groups_BBKNN.ipynb, is the scanpy workflow to integrated analysis of the scRNA dataset derived from different lesions, and compare the tissuse specific stromals.

5.1 Menstrual_dynamics_ourdata.R.ipynb, is the workflow to analyze the gene expression dynamic along the menstrual cycle for eutopic endometrial stromal cells, thereby, get the gene list to define the menstrual features.
5.2 NMdata_menstrual_score.R.ipynb, is the workflow to evaluate the menstrual feature gene list with the downloaded control endometrium data from healthy donors.
5.3 E-MTAB-10287_NGEU_reanalysis_menstrual.R.ipynb, is the workflow to evaluate the expression of our identified menstrual feature gene list in E-MTAB-10287 dataset.
5.4a Fonseca_NGdata_all_EnS.ipynb, is the workflow to subset all the endometrial stromal cells from Fonseca's dataset, prepared for the menstrual featrue gene expression score analysis.
5.4b Fonseca_EnS_menstrual_score.ipynb, is the workflow to evaluate the gene expression score with Fonseca's dataset using our defined menstrual feature gene list.
5.5 myEnS_menstrual_score_diffgene.R.ipynb, is the workflow to evaluate the menstrual feature genes expression score and identify the overexpressed genes for the ectopic endometrial stromal cells across the menstrual cycle.
5.6 Menstrual_gene_dotplot.ipynb, is the workflow to plot the selected gene expression for our data and Fonseca's data.

6.Fonseca_NG_OSC_reanalyze.R.ipynb, is the workflow to reanalyze the Ovarain Stromal Cells derived from the Fonseca's data, revealing the characteristics of the different subgroups of OSC.
7.multsamle_ST_seurat_SCT.R.ipybn,
