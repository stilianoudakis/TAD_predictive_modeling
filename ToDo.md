- `table1_genomic_elements.csv` - Make supplemental. Make up-to-date (some features you don't use). Add column "Cell type". Add column description to the legend.
- `figure2_predictor_types.png` - Good.
- `figure3_ensemble_framework.png` - Where is resampling shown?
- 50 iterations were done for ROS and RUS. Wasn't it done for SMOTE?
- TBD - how iterations and cross-validation play together?

- `figure18_model_construction_V2.png` - Write legend
- `table2_class_distributions_5kb.csv` - Add column description to the legend.
- `table8_class_distributions_other_resolutions.csv` - Add column header. Add legend.

- Make Figures and Tables in order
    - Rename files. 
    - Combine some figures into one, with A, B, C panels
    - Adjust figure/table legends
    - Adjust the text

- Add ROSE R package to the class-balance techniques
- Make supplementary figures of varimps for OC and OP predictor types
- Justify in text why not using ROC curves
- Understand the paper and Figure 5C (justification of the number of features kept in the model) from https://www.nature.com/articles/ng.3539 - can (should) we generate similar?
- Understand Figure 6 and 7 from https://www.nature.com/articles/ng.3539 - predictive importance of features, will be OC, OP, and distance in our case
- We should investigate these observations: Blagus and Lusa, “Class Prediction for High-Dimensional Class-Imbalanced Data.” - Variable selection introduces additional bias towards classification into the majority class. Undersampling helps, oversampling does not. Variable normalization (centering) can worsen the performance.
- Most important features - are they the same at different resolutions?
- Revisit comparing results with mourad (take his finding and compare; in discussion discuss differences)
- Add discussion; discuss literature: previous results, recommendations


Backburner: 

- Use tfbs encode clustered in separate model (/home/sequencing/data/ExtData/db_5.00_07-22-2015/grsnp_db/hg19/ENCODE/). Goal - to compare the results with Hong et. al 2017 paper, https://doi.org/10.1093/nar/gkx738
- Predict new TAD boundaries and investigate (visually, analytically) why they were not detected by the arrowhead method
- TAD boundaries can be hierarchical. Goal - use your framework to model level 1, 2, and 3 hierarchical boundaries, compare features. To discuss together with Kellen.
- Apply model to other cell lines in order to identify shared and unique features
- Combine TAD boundaries from all cell lines and predict them using features common to all cell lines to estimate performance and identify universal set of features
- Add more classifiers: GBM, SVM, Naive Bayes, neural network, MLR, Elastic-Net 
 