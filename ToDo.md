# ToDo

+ fix supplementary table 3 with all model performances (order predictor type column)
- What are the parameters of elastic net? 
    - Visual justification of the selected parameters?
    * plot of all alpha by lambda values creates to cluttered of a plot; discuss further in meeting

## Main Ideas

+ What type of precictor has the highest variable importance. Start with CTCF and down the list.
    + Is distance universally good for all genomic annotations? Test each annotation using OC, OP, and distance as an individual model, check
    
+ How model performance depend on the type of genomic annotations? Evaluate the performance when using 1) histone marks (+ DNAse I hypersensitive sites) only, 2) TFs only, 3) chromatin states only.
    + start with tfbs

+ Distance is the most important predictor. Thus, distance distributions for the most predictive features should be significantly different (KS-test?) from random
    + test if different compared to random bins; stratify by annotation type
    
+ Train the best model in one cell type, apply to other cell types
    + rerun models on clustered TFBS
    + compare performance
  
- Combine TAD boundaries from all cell lines (Intersection, or union) and predict them using features common to all cell lines to estimate performance and identify universal set of features
    - Use cluster specific TFBS
    - wait with chromstates and histone modifications
    - use common flanked TADs between two cell lines; how many are common?
    * Question: how to handle non uniform overlaps? what size bins to use for non TADs?
  

## Secondary Ideas

- Look into using grouped lasso to see if a particular group of annotations get regularized 

- Even the best model cannot predict all TADs accurately. Why some TADs are failing to be predicted?
    - look at profiles of annotations for tads that cant be predicted; wilcoxon tests
    - cross section with tad lengths

- Consider filtering out TAD boundaries that form TADs greater than 2mb bp in length prior to modelling and compare

- Use stacking/ensemble model building techniques and compare performances with RF, GBM, etc
    - "Comparison of Bagging, Boosting and Stacking Ensembles Applied to Real Estate Appraisal" shows that stacking outperforms other aggregation techniques like additive regression and bagging

- Prediction of differential TADs. Idea from Crow, Megan, Nathaniel Lim, Sara Ballouz, Paul Pavlidis, and Jesse Gillis. “Predictability of Human Differential Gene Expression.” Proceedings of the National Academy of Sciences 116, no. 13 (March 26, 2019): 6491–6500. https://doi.org/10.1073/pnas.1802973116.

- Check if we can adapt F-racing strategy: Dal Pozzolo, Andrea, Olivier Caelen, Serge Waterschoot, and Gianluca Bontempi. “Racing for Unbalanced Methods Selection.” In International Conference on Intelligent Data Engineering and Automated Learning, 24–31. Springer, 2013. - F-racing strategy to select best performing method to deal with class imbalance. Overview of class imbalance techniques, including SMOTE, Ensemble methods. The F-Racing approach tests in parallel a set of alternatives and uses Friedman test to determine if an alternative is significantly worse than others.Random Forest and SMOTEnsemble generally perform best. race R package to perform F-racing algorithm https://cran.r-project.org/web/packages/race/index.html

- Add more classifiers: GBM, SVM, Naive Bayes, neural network, MLR, Elastic-Net 

- Explore https://github.com/scikit-learn-contrib/imbalanced-learn - A Python Package to Tackle the Curse of Imbalanced Datasets in Machine Learning http://imbalanced-learn.org

- Use SpectralTAD to detect hierarchical TAD boundaries. Goal - predict level 1, 2, and 3 hierarchical boundaries, identify common and hierarchy-specific predictive features. 
    - Before predicting hierarchical boundaries separately, use combined set of TADs



## TBD 

- directly compare mourad using percent predictor types
- directly compare hong using count predictor types
- Compare proteins used by Mourad with our tfbs 
- investigate what features mourad used
   * We simulated data that were similar to real ChIP-seq data (see Subsection Materials and
     Methods, Data simulation, first paragraph). Both genomic coordinate data (e.g., ChIP-seq peak
     coordinates) and quantitative data (e.g., ChIP-seq signal intensity log ChIP
     Input) were generated.
   * We first simulated 100 genomic coordinate and 100 quantitative datasets
     that comprised 6 proteins and learned models without considering any interaction terms.
   * We also compared MLR with ET and RF using real data in human. For this purpose, we analyzed new 3D domains detected from           recent high resolution Hi-C data at 1 kb for GM12878 cells for which 69 ChIP-seq data were available
   * For human analysis, we used publicly available ChIP-seq peaks of 69 chromatin proteins (ATF2, ATF3, BATF, BCL11A, BCL3,            BCLAF1, BHLHE40,       BRCA1, CEBPB, CHD1, CHD2, CTCF, E2F4, EBF1, EGR1, ELF1, ELK1, ETS1, EZH2, FOS, FOXM1, IKZF1, IRF3,          IRF4, JUND, MAFK, MAX, MAZ, MEF2A,          MEF2C, MTA3, MXI1, MYC, NFATC1, NFE2, NFIC, NFYA, NFYB, NRF1, P300, PAX5, PBX3,        PIGG, PML, POU2F2, RAD21, REST, RFX5, RUNX3, RXRA, SIN3A,      SIX5, SP1, SRF, STAT1, STAT3, STAT5A, TAF1, TCF12, TCF3, USF1,      USF2, YY1, ZBTB33, ZEB1, ZNF143, ZNF274, ZNF384 and ZZZ3) of GM12878 cells      from ENCODE 


# Papers to describe and cite (bold highlights relevant point)

- Dixon, Jesse R., Inkyung Jung, Siddarth Selvaraj, Yin Shen, Jessica E. Antosiewicz-Bourget, Ah Young Lee, Zhen Ye, et al. “Chromatin Architecture Reorganization during Stem Cell Differentiation.” Nature 518, no. 7539 (February 19, 2015): 331–36. https://doi.org/10.1038/nature14222. - Reorganization of the 3D genome during human embryonic stem cell differentiation. ~36% of active/inactive A/B compartments change. Fig. 1 - how to compare PC1 vectors. Gene expression changes in compartments that switch states (up/down in B to A vs. A to B switching), small but significant. Changes in _interaction frequency_ (increase/decrease - switch from B to A vs. A to B), correlated positively with activating epigenomic marks H3K27ac, DNAse hypersensitive sites, CTCF, negatively with repressive marks H3K27me3 and H3K9me3. **Random Forest to select epigenomic features classifying changes in interaction frequency - H3K4me1 enhancer mark**. Haplotype-resolved Hi-C maps are similar.



