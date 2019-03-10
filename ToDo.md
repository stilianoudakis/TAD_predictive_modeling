
- Check if we can adapt F-racing strategy: Dal Pozzolo, Andrea, Olivier Caelen, Serge Waterschoot, and Gianluca Bontempi. “Racing for Unbalanced Methods Selection.” In International Conference on Intelligent Data Engineering and Automated Learning, 24–31. Springer, 2013. - F-racing strategy to select best performing method to deal with class imbalance. Overview of class imbalance techniques, including SMOTE, Ensemble methods. The F-Racing approach tests in parallel a set of alternatives and uses Friedman test to determine if an alternative is significantly worse than others.Random Forest and SMOTEnsemble generally perform best. race R package to perform F-racing algorithm https://cran.r-project.org/web/packages/race/index.html



Backburner: 

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

- Predict new TAD boundaries and investigate (visually, analytically) why they were not detected by the arrowhead method
- TAD boundaries can be hierarchical. Goal - use your framework to model level 1, 2, and 3 hierarchical boundaries, compare features. To discuss together with Kellen.
- Apply model to other cell lines in order to identify shared and unique features
- Combine TAD boundaries from all cell lines and predict them using features common to all cell lines to estimate performance and identify universal set of features
- Add more classifiers: GBM, SVM, Naive Bayes, neural network, MLR, Elastic-Net 
- Explore https://github.com/scikit-learn-contrib/imbalanced-learn - A Python Package to Tackle the Curse of Imbalanced Datasets in Machine Learning http://imbalanced-learn.org
 