## For Grubert results 

### Spiro

- Tables with all performance metrics. Do we have them?

- Table with performance metrics of other machine learning methods (SVM, XGBoost, etc.)

- Per chromosome model performance. Do we have them?

- **Additional File 11: Figure S6. Comparing enrichment levels between TAD/chromatin loop calling tools.** Signal profile plots comparing the binding strength of top TFBS around flanked Arrowhead called TAD boundaries (blue), Peakachu chromatin loop boundaries, Grubert chromatin loop boundaries (purple), and SpectralTAD called TAD boundaries (green) on (A) GM12878 and (B) on K562.
    - Are you using midpoints of boundaries? Important because it is a very strong statement.

- **Additional File 17: Figure S12. *preciseTAD* trained on Arrowhead accurately predicts boundaries on cell lines using annotation data only.** Venn diagrams and signal profile plots comparing flanked predicted boundaries using Arrowhead trained models. (A) Models trained on GM12878 and predicted on GM12878 (red, GM on GM) vs. models trained on K562 and predicted on GM12878 (blue, K on GM). (B) Models trained on K562 and predicted on K562 (red, K on K) vs. models trained on GM12878 and predicted on K562 (blue, GM on K). Boundaries were flanked by 5 kb.
    - Need BED files for PTBRs from "K on G" and "G on K" predictions.
    - Same for **Additional File 18: Figure S13. *preciseTAD* trained on Peakachu accurately predicts boundaries on cell lines using annotation data only.**

- Figure 6, Additional_File_10_figureS5 - get rid of panel C; instead, add Grubert signal (original midpoints of Grubert anchors) to each plot on panels A and B

- Figure 5, panel A - make more editable, remove boundaries from the Hi-C triangle, ungroup vertical lines on the base-level probability panel.

- Figure that compares model performance using combinations of annotations. Needed for predicting boundaries in all cell types. `26_dnase_vs_topTFBS.R` - keep transcription factor only, all combinations of them.

### MD

- Figure 1 - make new (MD)
- What is the resolution (width) of Grubert peaks? All genomic annotation data? BroadHMM in particular (MD/Maggie, Spiro?)
    - How Grubert peaks overlap with CTCF, RAD21 (native Chia-PET factor), SMC3? Jaccard matrix.
- Jaccard overlap between Arrowhead/Peakachu boundaries and Grubert peaks? Venn diagram. (MD/Maggie, Spiro?)
- To Be Decided. Train the model using Grubert boundaries and Grubert as ground truth, using CTCF/RAD21/SMC3/ZNF143.
    - What's the order of feature importance for those?
- Figure 5B - colocalization with CTCF etc. peaks is not intuitive.

- Update Spiro on 1) Avocado predictions (Maggie) and 2) base-level predictions

- Remove:
    - Additional_File_12_figureS7.png
    - Additional_File_13_figureS8.png
    - Additional_File_14_figureS9.png
    


### MD notes

Grubert et al. generated cell line-specific ChIA-PET chromatin interaction maps at single-cohesin peak resolution (about 2kb??). They chose the RAD21 subunit of the cohesin complex, which facilitates physical contacts between genes and enhancers, and CTCF binding sites [@Heidari:2014wr]. Their loops are similar in size to those identified by Rao et al. 2014. 

The increasingly complex definitions of TADs, compartments, and loops indicate that distinct classes of chromatin domains overlap and share some functional, mechanistic, and structural properties (see Beagan & Phillips-Cremins [@Beagan:2020wz] for an excellent review of the current knowledge about chromatin compartments).

Rao et al. 2014 detected "contact domains" of median length 185kb [@Rao:2014aa], 3-5 times smaller than the 1Mb "Topologically Associated Domains" detected by Dixon et al. 2012. They are similar in size to "physical domains" detected in Drosophila [@Sexton:2012aa] and "topologically constrained domains" in structural studies of human chromatin [@Cook:1975ve][@Vogelstein:1980uf]. Boundaries of those contact domains were enriched in loops, suggesting they are "loop domains" following the definition of Beagan & Phillips-Cremins (A contact domain formed via loop-extrusion mechanisms and often but not always having a dot at the corner (corner-dot domain), owing to the transient nature of extrusion (non-compartment or non-corner-dot domain)) [@Beagan:2020wz].

Most peak loci encompass a unique DNA site containing a CTCF-binding motif, to which all three proteins (CTCF, SMC3, and RAD21) were bound (5-fold enrichment) [@Rao:2014aa].

Although CTCF and cohesin binding tend to colocalize, cohesin peaks are slightly shifted to the 3' ends of convergently oriented motifs [@Tang:2015aa][@Fudenberg:2016aa].

Numerous studies demonstrated the association of domain boundaries with genomic annotations [@Dixon:2015aa; @Whalen:2016aa; @Dekker:2016aa; @Bonev:2016aa; @Wang:2016ab].

Heidari et al. identified cohesin, CTCF, and ZNF143 as key components of the 3D chromatin structure [@Heidari:2014wr]. ZNF143 [@Bailey2015]. Boundaries of TADs were found to be enriched in architectural factors such as CTCF, RAD21, SMC3, YY1, and ZNF143 [@Dixon:2012aa; @Rao:2014aa; @Phillips-Cremins2013:aa; @Corces:2016aa]

subTAD boundaries exhibit weaker insulation strength, as evidenced by their relatively lower capacity to attenuate long-range contacts between domains, and they are also significantly more likely than TADs to exhibit cell-type-dynamic folding properties [@Dixon:2012aa][@Phillips-Cremins2013:aa][@Norton:2018aa].

Cohesin translocation extrudes DnA in an ATP-dependent manner into long-range looping interactions that form the topological basis for TAD and subTAD loop domains.

Boundaries of TADs were found to be enriched in architectural factors such as CTCF, RAD21, SMC3, YY1, and ZNF143 [@Dixon:2012aa; @Rao:2014aa; @Phillips-Cremins2013:aa; @Corces:2016aa], and boundary strength correlates with their occupancy [@Van-Bortle:2014aa; @Fudenberg:2019ab]. Strong and weak boundary concept [@Chang:2020vw].


##### Reviewer's suggestion

- Doug Phanstiel
- Yun Li
- Ferhat Ay
- Sushmita Roy
- Sunduz Keles
- Anthony Schmidt

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



