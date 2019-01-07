- how to select testing set in class balance papers when validating models; do they balance test set also?
- investigate further why using HIC001-HIC029
- figure or table to validate log2 transform distances
- include all predictor types, use lasso and report which type is best
- include feature engineering step in model flowchart figure
+ `table1_genomic_elements.csv` - Make supplemental. Make up-to-date (some features you don't use). Add column "Cell type". Add column description to the legend.
+ `figure2_predictor_types.png` - Good.
- `figure3_ensemble_framework.png` - Where is resampling shown?
+ 50 iterations were done for ROS and RUS. Wasn't it done for SMOTE?
- TBD - how iterations and cross-validation play together?
- `figure4_class_imbalance_plots.png` - in "Over-Sampling", is the majority class the same as in "Complete" case? Looks smaller.
  * it is the same as the majority class for the training set
+ `figure18_model_construction_V2.png` - Add legend
+ `figure13_other_resolutions_none_results.png` - Add legend
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


- Make Figures and Tables in order
    - Rename files. 
    - Combine some figures into one, with A, B, C panels
    - Adjust figure/table legends
    - Adjust the text

- Use the numerical quantity "imbalance ratio" (IR), when first using, cite A. Orriols-Puig, E. Bernad-Mansilla, "Evolutionary rule-based systems for imbalanced data sets", Soft Comp., vol. 13, pp. 213-225, 2009.
- Add ROSE R package to the class-balance techniques
- Make supplementary figures of varimps for OC and OP predictor types
- Justify in text why not using ROC curves
+ Understand the paper and Figure 5C (justification of the number of features kept in the model) from https://www.nature.com/articles/ng.3539 - can (should) we generate similar?
+ Understand Figure 6 and 7 from https://www.nature.com/articles/ng.3539 - predictive importance of features, will be OC, OP, and distance in our case
+ We should investigate these observations: Blagus and Lusa, “Class Prediction for High-Dimensional Class-Imbalanced Data.” - Variable selection introduces additional bias towards classification into the majority class. Undersampling helps, oversampling does not. Variable normalization (centering) can worsen the performance.
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
 