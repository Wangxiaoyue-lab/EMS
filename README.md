Code for analyze the single cell transcriptome data of endometrisois

1. The file alldata2024.ipynb, is the scanpy workflow used to analyze our data, and get the cell clusters and cell types.
2. OVcancer_scanpy.ipynb, is the workflow to reanalyze the genes specifically expressed in stromal cells derived from different types of ovarain cancer.
3. rds2h5ad.r.ipynb, used to convert the downloaded Fonseca's data (https://cedars.app.box.com/s/1ks3eyzlpnjbrseefw3j4k7nx6p2ut02) into h5ad format.
4. NGdata_scanpy_OMA.ipynb, is the workflow to analyze the ovarian endometriosis data, a subset of Fonseca' data.
5. NGdata_PEM_EU_OV_scanpy.ipynb, is the workflow to analyze the Fonseca's data subset according to the derived tissue.
6. NGdata_all_EnS.ipynb, is the workflow to extract all the endometrial stromal cells from eutopic endometrium,peritoneal endometriosis, and endometrioma.
7. Menstrual_dynamics_ourdata.R.ipynb, is the workflow to analyze the gene expression dynamic along the menstrual cycle for eutopic endometrial stromal cells, thereby, get the gene list to define the menstrual features.
8. myEnS_menstrual_score_diffene.R.ipynb, is the workflow to evaluate the menstrual feature genes score and discover the overexpressed genes for the ectopic endometrial stromal cells across the menstrual cycle.
9. NMdata_menstrual_score.R.ipynb, is the workflow to evaluate the menstrual feature gene list with the downloaded control endometrium data from healthy donors.
10. Menstrual_gene_dotplot.ipynb, is the workflow to evalute the selected menstrual gene expression for our data and Fonseca's data.
