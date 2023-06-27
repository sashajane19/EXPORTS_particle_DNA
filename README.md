# EXPORTS_particle_DNA
Processing code for analysis of DNA sequences collected from surface water samples, in bulk sediment traps, and from individual sinking particles as part of the EXport Processes in the Ocean from RemoTe Sensing (EXPORTS) field campaigns. DNA extraction and 18S rRNA gene amplification details can be found in Durkin et al. (2022).

Processing of 18S rRNA gene sequence files was conducted following Durkin et al. (2022) using QIIME2 (Bolyen et al., 2019) and the DADA2 (Callahan et al., 2016) workflow. Taxonomic identities were assigned to amplicon sequence variants (ASVs) based on comparison to the PR2 database (Guillou et al., 2013; v 4.14.0, downloaded September 2022) with a naïve Bayes classifier (Bokulich et al., 2018). Taxa were assigned to photosynthetic or heterotrophic following Supplementary Table 2 in Durkin et al. (2022) and adding photosynthetic taxa found in this dataset that were not found in that dataset (photosynthetic_taxa.csv). Given the large variability in sequence counts across different sample types, rarefaction of the data was also performed using QIIME2 using only the photosynthetic taxa (Weiss et al., 2017). 

Please reference the following work if you use this code to inform your own work:

Bokulich, N. A., Kaehler, B. D., Rideout, J. R., Dillon, M., Bolyen, E., Knight, R., et al. (2018). Optimizing taxonomic classification of marker-gene amplicon sequences with QIIME 2’s q2-feature-classifier plugin. Microbiome, 6(1), 90. https://doi.org/10.1186/s40168-018-0470-z

Bolyen, E., Rideout, J. R., Dillon, M. R., Bokulich, N. A., Abnet, C. C., Al-Ghalith, G. A., et al. (2019). Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology, 37(8), 852–857. https://doi.org/10.1038/s41587-019-0209-9

Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: high-resolution sample inference from Illumina amplicon data. Nature Methods, 13(581), 581–583. https://doi.org/10.1038/nMeth.3869

Durkin, C., Cetinić, I., Estapa, M. L., Ljubešić, Z., Mucko, M., Neeley, A., & Omand, M. M. (2022). Tracing the path of carbon export in the ocean though DNA sequencing of individual sinking particles. The ISME Journal, 1–11. https://doi.org/10.1038/s41396-022-01239-2

Guillou, L., Bachar, D., Audic, S., Bass, D., Berney, C., Bittner, L., et al. (2013). The Protist Ribosomal Reference database (PR2): a catalog of unicellular eukaryote Small Sub-Unit rRNA sequences with curated taxonomy. Nucleic Acids Research, 41(D1), D597–D604. https://doi.org/10.1093/nar/gks1160

Weiss, S., Xu, Z. Z., Peddada, S., Amir, A., Bittinger, K., Gonzalez, A., et al. (2017). Normalization and microbial differential abundance strategies depend upon data characteristics. Microbiome, 5(1), 27. https://doi.org/10.1186/s40168-017-0237-y
