## scRNA-seq data analysis: DOT1L activity affects cell lineage progression in the developing brain by controlling metabolic programs (Appiah et al)
Documentation of workflow used for analysing scRNA-seq data from labelled cells acquired following microinjection or electroporation in embryonic mouse brain (E14.5). Once published, the original data will be available for download at GEO via accession number: GSE176323.

## Brief description of data acquisition procedure
Our aim was to understand the role of DOT1L in cortical development. DOT1L functions as both an enzyme and a scaffolding protein. In these experiments, effects of DOT1L inhibition on cell fate decisions were assessed.

Apical progenitors in embryonic mouse brain (E14.5) were labelled either via microinjection of Dextran coupled to Alexa fluorophore or via ex utero electroporation. Subsequently, the DOT1L inhibitor pinometostat (EPZ5676) or vehicle (DMSO) was applied to organotypic slices or electroporated hemispheres in culture for 24 h. At the end of culture, tissues were dissociated, labelled cells were sorted into 384 well plates and processed for scRNA-seq using the mCEL-Seq2 protocol (Herman et al., 2018).

## Reproducing the analysis - important notes
To run the scripts in the markdown files smoothly, define the working directory by using knitr:opts_knit$set(root.dir = 'INSERT PATH TO DIRECTORY').
Specific versions of packages used must be installed properly. Environment in which analyses were performed is provided (both as .yml and text files) and commented on below:

sc3.yml: Data pre-processing, clustering, and downstream analysis and SCENIC anaysis were performed in this conda environment.

Overview of the main seurat analysis has been provided in the folder: Overview_of_Main_Seurat_Analysis_Appiah_et_al. To view how the figures in the markdown files should look like, download the files as html and view them in your browser.

### Additional information

Platform: x86_64-conda_cos6-linux-gnu (64-bit)

About software:

Run under: Ubuntu 18.10 

R version 3.6.1 (2019-07-05)

Python 3.7.6

anaconda-client (version 1.7.2)

Conda 4.8.3

Seurat v3

RaceID3





