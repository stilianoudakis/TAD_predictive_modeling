
## For Package

- When running the last example in the README file (line 120) and the cell-type-specific predictions in the vignette (line 318), the code breaks with error:
```
Error in makePSOCKcluster(names = spec, ...) : 
  Cluster setup failed. 2 of 2 workers failed to connect.
```
The error is in `preciseTAD()` function, when the `parallel` argument is a number (`=2`, `=1`). The error occurs on Mac. Setting `parallel = NULL` allows error-free execution. The error does not occur on Windows, but multiple CPUs seem not to be utilized. Something wrong with parallelization, investigate.


## For publication

- Prepare submission material folder, [Genome Biology](https://genomebiology.biomedcentral.com/submission-guidelines/preparing-your-manuscript/software). Check:
    - Submitting to special issue: https://genomebiology.biomedcentral.com/regulatory-elements
    - Suggested reviewers:
        + Mattia Forcato mattia.forcato@unimore.it
        + Kadir Akdemir, kcakedemir@mdanderson.org
        + Geoff Fudenberg, geoff.fudenberg@gladstone.ucsf.edu
        - Kazuhiro Maeshima kmaeshim@nig.ac.jp
        - Biola M. Javierre bmjavierre@carrerasresearch.org
        - Anton Goloborodko, goloborodko.anton@gmail.com
    
- New experiment, priority on par with ExpermentalHub
    - Gene expression for boundary prediction. 
        - Start of a gene +/-500bp defines feature
        - Gene expression defines signal
        - Same questions - whether overlap count, overlap percent, signal, and distance are predictive of TAD/loop boundaries

- Look into ExperimentalHub

- look up “calibration curve machine learning” and think if/how we can apply it

- Need either repo or all scripts generating all figures, like https://www.biorxiv.org/content/10.1101/2020.08.17.254615v1.supplementary-material. For preprint, Additional file will suffice
    - Add README.md. For each script, specify input, what it does, what is the output

- Look into permutation test as replacement for normalized overlap. Quantify permutation enrichment p-values for two things:
    - Supplementary Figure 9, when described. Also, Discussion
    - Figure 5 B,C (and supplementary 5) - to statistically show that CTCF et al. are more significantly enriched.





### Writing tips

- Write using the structure (https://github.com/mdozmorov/manuscript_template):
    - What is the question?
    - What data/methods are used?
    - Any relevant details/results?
    - What is the answer (take home message)?

## Goals

* TAD boundaries are detected at low resolution, thousands of kilobases
* Genome was annotated with vast amount of cell type-specific annotations
    * Annotated regions are small, on the order of tens-hundreds bases
* Precisely predict TAD boundaries using genomic annotations using a machine learning model

## Problems to be addressed

- Class imbalance - many more ordinary regions than TAD boundaries
- Feature engineering - different types of association between TAD boundaties and genomic annotations
- Metric - Which metric to use to evaluate model performance in class imbalance settings?
- Resolution - How the model performs across resolutions
- Collinearity - Many genomic annotations are colocalize (i.e., collinear). How to properly address it?

## Current methods

- Class imbalance - random undersampling, random oversampling, SMOTE
- Feature engineering - OC, OP, Distance
- Model - currently, Random Forest as a balanced combination of performance and interpretability
- Metric - currently, Matthew Correlation Coefficient
- Resolution - currently, 10/25,50/100kb
- Collinearity - currently, addressed by recursive feature elimination
- Using optimally trained model, forget about TAD boundary locations used for training and predict the location of TAD boundaries only from genome annotation data

## ToDo, Kellen

- Orient within Methods
    - Data, Supplementary Table S1 - Add column specifying what type of replicate (technical/biological) the sample is. 
    - Find how the data was processed. Use non-normalized and the combined file (Need to ask Spiro exact)
    - Where the processed data is located? We need 5/10/25,50/100kb data

- Call TADs using consensus score across replicates
    - Check the number of TADs, Spiro obtained very few TADs using consensus score
    - Need to update Table 1

- Orient in genome annotation data, Supplementary Table S2
    - Where the data is located?

- Orient within code Spiro created
    - Feature engineering code - where it is?
    - Class imbalance code - where it is?
    - Make README.md annotating each file

### New results

- CTCF and other selected genomic annotations were closer to the predicted than to the Hi-C identified TAD boundaries. Test the significance of it using GenometriCorr.

- For predicted regions: provide summaries of widths

- Predict precise TAD boundary locations for other tissue/cell types using the framework developed above
    - Release them online, as in http://promoter.bx.psu.edu/hi-c/publications.html

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



