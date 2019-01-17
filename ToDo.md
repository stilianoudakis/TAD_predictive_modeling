
+ clearly define framework
- redo figures/tables 
- directly compare mourad using percent predictor types
- directly compare hong using count predictor types
- Use the numerical quantity "imbalance ratio" (IR), when first using, cite A. Orriols-Puig, E. Bernad-Mansilla, "Evolutionary rule-based systems for imbalanced data sets", Soft Comp., vol. 13, pp. 213-225, 2009.
- Compare proteins used by Mourad with our tfbs 
+ figure or table to validate log2 transform distances 
   * have figure; is figure enough? what type of table to include? supplementary?
- investigate what features mourad used
   * We simulated data that were similar to real ChIP-seq data (see Subsection Materials and
     Methods, Data simulation, first paragraph). Both genomic coordinate data (e.g., ChIP-seq peak
     coordinates) and quantitative data (e.g., ChIP-seq signal intensity log ChIP
     Input) were generated.
   * We first simulated 100 genomic coordinate and 100 quantitative datasets
     that comprised 6 proteins and learned models without considering any interaction terms.
   * We also compared MLR with ET and RF using real data in human. For this purpose, we analyzed new 3D domains detected from           recent high resolution Hi-C data at 1 kb for GM12878 cells for which 69 ChIP-seq data were available
   * For human analysis, we used publicly available ChIP-seq peaks of 69 chromatin proteins (ATF2, ATF3, BATF, BCL11A, BCL3,            BCLAF1, BHLHE40,       BRCA1, CEBPB, CHD1, CHD2, CTCF, E2F4, EBF1, EGR1, ELF1, ELK1, ETS1, EZH2, FOS, FOXM1, IKZF1, IRF3,          IRF4, JUND, MAFK, MAX, MAZ, MEF2A,          MEF2C, MTA3, MXI1, MYC, NFATC1, NFE2, NFIC, NFYA, NFYB, NRF1, P300, PAX5, PBX3,        PIGG, PML, POU2F2, RAD21, REST, RFX5, RUNX3, RXRA, SIN3A,      SIX5, SP1, SRF, STAT1, STAT3, STAT5A, TAF1, TCF12, TCF3, USF1,      USF2, YY1, ZBTB33, ZEB1, ZNF143, ZNF274, ZNF384 and ZZZ3) of GM12878 cells      from ENCODE 
- learn GIMP




- TBD - how iterations and cross-validation play together?

- `figure14_5kb_rus_results.png` - Add legend
- `figure15_other_resolutions_rus_results.png` - Add legend
- `figure16_5kb_dist_results.png` - Add legend
- `figure10_allresolutions_varimps.png` - To be redone with the correct reatures.

- `table2_class_distributions_5kb.csv` - Add column description to the legend.
- `table8_class_distributions_other_resolutions.csv` - Add column header. Add legend.
- `table9_5kb_none_performances.csv` - Add legend. Add other metrics (TP, FP etc.)
- `table10_other_resolutions_none_performances.csv` - Add legend. Add other metrics (TP, FP etc.)
- `table11_5kb_rus_performances.csv` - Add legend. Add other metrics (TP, FP etc.)
- `table12_other_resolutions_rus_performances.csv` - Add legend. Add other metrics (TP, FP etc.)
- `table13_5kb_dist_performances.csv` - Add legend. Add other metrics (TP, FP etc.)



- Make supplementary figures of varimps for OC and OP predictor types
- Justify in text why not using ROC curves
- Understand Figure 6 and 7 from https://www.nature.com/articles/ng.3539 - predictive importance of features, will be OC, OP, and distance in our case
- Most important features - are they the same at different resolutions?
- Revisit comparing results with mourad (take his finding and compare; in discussion discuss differences)
- Add discussion; discuss literature: previous results, recommendations
- Explore https://github.com/scikit-learn-contrib/imbalanced-learn - A Python Package to Tackle the Curse of Imbalanced Datasets in Machine Learning http://imbalanced-learn.org

Backburner: 

- Use tfbs encode clustered in separate model (/home/sequencing/data/ExtData/db_5.00_07-22-2015/grsnp_db/hg19/ENCODE/). Goal - to compare the results with Hong et. al 2017 paper, https://doi.org/10.1093/nar/gkx738
- Predict new TAD boundaries and investigate (visually, analytically) why they were not detected by the arrowhead method
- TAD boundaries can be hierarchical. Goal - use your framework to model level 1, 2, and 3 hierarchical boundaries, compare features. To discuss together with Kellen.
- Apply model to other cell lines in order to identify shared and unique features
- Combine TAD boundaries from all cell lines and predict them using features common to all cell lines to estimate performance and identify universal set of features
- Add more classifiers: GBM, SVM, Naive Bayes, neural network, MLR, Elastic-Net 
 