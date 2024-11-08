# Colibactin-driven colon cancer requires adhesin-mediated epithelial binding

# About
Various bacteria are suggested to contribute to colorectal cancer (CRC) development including pks+ Escherichia coli, which produces the genotoxin colibactin that induces characteristic mutational signatures in host epithelial cells. However, it remains unclear how the highly unstable colibactin molecule is able to access host epithelial cells to cause harm. Here, using the microbiota-dependent ZEB2-transgenic mouse model of invasive CRC, we demonstrate that the oncogenic potential of pks+ E. coli critically depends on bacterial adhesion to host epithelial cells, mediated by the type 1 pilus adhesin FimH and the F9 pilus adhesin FmlH. Blocking bacterial adhesion using a pharmacological FimH inhibitor attenuates colibactin-mediated genotoxicity and CRC exacerbation. We also show that allelic switching of FimH strongly influences the genotoxic potential of pks+ E. coli and can induce a genotoxic gain-of-function in the probiotic strain Nissle 1917. Adhesin-mediated epithelial binding subsequently allows the production of the genotoxin colibactin in close proximity to host epithelial cells, which promotes DNA damage and drives CRC development. These findings present promising therapeutic routes for the development of anti-adhesive therapies aimed at mitigating colibactin-induced DNA damage and inhibiting the initiation and progression of CRC, particularly in individuals at risk for developing CRC. 

# Content of this repository
In this repository you can find all the code used for the Bulk RNA-seq data analysis. Below you can find an overview:

- 1.Script_DESeq2_Pre_processing.R: Allows you to get the original DESeq2 Object used for all the analysis.
- 2.Script_Functions.R: Contains all functions used for the analysis.
- 3.Script_Figure_2_Supplementary_Figure_1G.R: Script used for generating Figure 2 and Supplementary Figure 1G.
- 4.Script_Supplementary_Figure_1_C_D_I_J_K.R: Script used for generating Supplementary Figure 1(C,D,I,J,K).
- 5.session_info.txt: Contains all packages used and their respective versions for the analysis.

Notes:
- To be able to run the provided R scripts you need the raw_counts.csv file provided on [GSE241851](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241851) and additional files found in the data folder.

# Citation

Published Manuscript: [http://img.shields.io/badge/DOI-10.1038/s41586-024-08135-z.svg](https://doi.org/10.1038/s41586-024-08135-z)
Preprint Manuscript: [![DOI:10.1101/2023.08.16.553526](http://img.shields.io/badge/DOI-10.1101/2023.08.16.553526-B31B1B.svg)](https://www.biorxiv.org/content/10.1101/2023.08.16.553526v1)

Code: [![DOI](https://zenodo.org/badge/682061416.svg)](https://zenodo.org/doi/10.5281/zenodo.10046232)
