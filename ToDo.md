# Legend

- To do items
+ Completed to do items (awaiting discussion)
* Outstanding questions

## Deliverables

- Class imbalance
- OC, OP, Distance
- Metric
- Resolution
- Collinearity
- TBD - Feature selection
- Prediction


## ToDo

- reconfigure sup fig 1

- fix ???

- add distance boxplots as last sup figs

+ put legends at bottom of manuscript

+ make heatmap larger in fig4/sup fig 8

+ fix lines on heatmaps fig4/sup fig 8

+ play around with jaccard heatmaps; consider making diagonal na's

+ change to chromhmm sup fig 7; put A and B; reorder chromhmm and hm

+ combine sup fig 4 and 5; fix xaxis

+ cut sup fig 3 to .2; put A and B

- use arial font

- reference kellens bioarxiv paper

- look at instructions for journals

- make figures in final form

- for predicted regions: provide summaries of widths

- look at distances between called and predicted regions

- look into developing a way to enumerate the strength of each predicted region

- write paragraphs of all previous methods; refer them to supplementary
   - Comparison with Mourad? Hong?

- Relevant supplementary figures should be referred consequentially.

### For later

10. Use SpectralTAD to call TADs. Read the preprint, the package vignette, talk with Kellen, if needed.
    - Extract Hi-C matrices in text format. Replicates should be in individual matrices.
    - Call consensus TADs. Ask Kellen how. This is analogous what Arrowhead does now.
    - Call TADs in individual replicates.

### Main Ideas

- Write each paragraph using the structure (https://github.com/mdozmorov/manuscript_template):
    - What is the question?
    - What data/methods are used?
    - Any relevant details/results?
    - What is the answer (take home message)?
    
    
### Secondary Ideas

- https://sites.google.com/site/raphaelmouradeng/programs

- Prediction of hierarchical TADs

- Prediction using signal, like gene expression/methylation level. Wang,Y. et al. (2016) Predicting DNA methylation state of CpG dinucleo- tide using genome topological features and deep networks.

- Prediction of differential TADs. Idea from Crow, Megan, Nathaniel Lim, Sara Ballouz, Paul Pavlidis, and Jesse Gillis. “Predictability of Human Differential Gene Expression.” Proceedings of the National Academy of Sciences 116, no. 13 (March 26, 2019): 6491–6500. https://doi.org/10.1073/pnas.1802973116.

- Add more classifiers: GBM, SVM, Naive Bayes, neural network, MLR, Elastic-Net 

- Explore https://github.com/scikit-learn-contrib/imbalanced-learn - A Python Package to Tackle the Curse of Imbalanced Datasets in Machine Learning http://imbalanced-learn.org

- TAD boundaries are not created equal: CTCF/Cohesin (Smc1) mark long-range interactions, Mediator (Med12)/Cohesin mark short-range interactions within and between larger subdomains [@Phillips-Cremins2013:aa]. TAD boundary strength correlates with occupancy of architectural proteins [@Van-Bortle:2014aa]. Perform PCA on TAD annotations to see if they separate into groups. Use A/B compartments as a proxy to separate TADs into two groups, predict then. Goal - to find differential signatures.

- Can size of TAD boundaries be an important predictor?


### TBD 

- investigate what features mourad used
   * We simulated data that were similar to real ChIP-seq data (see Subsection Materials and
     Methods, Data simulation, first paragraph). Both genomic coordinate data (e.g., ChIP-seq peak
     coordinates) and quantitative data (e.g., ChIP-seq signal intensity log ChIP
     Input) were generated.
   * We first simulated 100 genomic coordinate and 100 quantitative datasets
     that comprised 6 proteins and learned models without considering any interaction terms.
   * We also compared MLR with ET and RF using real data in human. For this purpose, we analyzed new 3D domains detected from           recent high resolution Hi-C data at 1 kb for GM12878 cells for which 69 ChIP-seq data were available
   * For human analysis, we used publicly available ChIP-seq peaks of 69 chromatin proteins (ATF2, ATF3, BATF, BCL11A, BCL3,            BCLAF1, BHLHE40,       BRCA1, CEBPB, CHD1, CHD2, CTCF, E2F4, EBF1, EGR1, ELF1, ELK1, ETS1, EZH2, FOS, FOXM1, IKZF1, IRF3,          IRF4, JUND, MAFK, MAX, MAZ, MEF2A,          MEF2C, MTA3, MXI1, MYC, NFATC1, NFE2, NFIC, NFYA, NFYB, NRF1, P300, PAX5, PBX3,        PIGG, PML, POU2F2, RAD21, REST, RFX5, RUNX3, RXRA, SIN3A,      SIX5, SP1, SRF, STAT1, STAT3, STAT5A, TAF1, TCF12, TCF3, USF1,      USF2, YY1, ZBTB33, ZEB1, ZNF143, ZNF274, ZNF384 and ZZZ3) of GM12878 cells      from ENCODE 


### Papers to describe and cite (bold highlights relevant point)

- Dixon, Jesse R., Inkyung Jung, Siddarth Selvaraj, Yin Shen, Jessica E. Antosiewicz-Bourget, Ah Young Lee, Zhen Ye, et al. “Chromatin Architecture Reorganization during Stem Cell Differentiation.” Nature 518, no. 7539 (February 19, 2015): 331–36. https://doi.org/10.1038/nature14222. - Reorganization of the 3D genome during human embryonic stem cell differentiation. ~36% of active/inactive A/B compartments change. Fig. 1 - how to compare PC1 vectors. Gene expression changes in compartments that switch states (up/down in B to A vs. A to B switching), small but significant. Changes in _interaction frequency_ (increase/decrease - switch from B to A vs. A to B), correlated positively with activating epigenomic marks H3K27ac, DNAse hypersensitive sites, CTCF, negatively with repressive marks H3K27me3 and H3K9me3. **Random Forest to select epigenomic features classifying changes in interaction frequency - H3K4me1 enhancer mark**. Haplotype-resolved Hi-C maps are similar.



