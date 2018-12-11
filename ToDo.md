
- fix figure 18 y axis
- clearly specify binning process; mention that we use contact matrix coordinates to obtain bins
- don't break up figures 18, 3 into parts
- for class imbalance make figure that tells story about how severe class imbalance is; currently figure 12 is evaluating predictors; comparisons need to be re-sampling 
- update table and figure legends
- fix csl path
- look at varimps for oc and op predictor types

##11/28/18


- add discussion; discuss literature: previous results, recommendations
+ fix figure 10; make 4 panels specifying resolution
+ make figure 10 into a table
- look at models with all predictor types 
- look at correlations amoung predictors for each type
- Read papers!
- Plot ROC curves
- Select features using Boruta
- Predict new TAD boundaries and investigate (visually, analytically)
- We should investigate these observations: Blagus and Lusa, “Class Prediction for High-Dimensional Class-Imbalanced Data.” - Variable selection introduces additional bias towards classification into the majority class. Undersampling helps, oversampling does not. Variable normalization (centering) can worsen the performance.

Backburner: 
- use tfbs encode clustered in separate model (/home/sequencing/data/ExtData/db_5.00_07-22-2015/grsnp_db/hg19/ENCODE/). Goal - to compare the results with Hong et. al 2017 paper, https://doi.org/10.1093/nar/gkx738
- look at chromosome specific models and compare
- revisit comparing results with mourad (take his finding and compare; in discussion discuss differences)
- Apply model to other cell lines in order to identify shared and unique features
- Combine TAD boundaries from all cell lines and predict them using features common to all cell lines to estimate performance and identify universal set of features
- GBM
- SVM
- Naive Bayes
- neural network
- MLR
- Elastic-Net 
- Understand the paper and Figure 5C (justification of the number of features kept in the model) from https://www.nature.com/articles/ng.3539 - can (should) we generate similar?
- Understand Figure 6 and 7 from https://www.nature.com/articles/ng.3539 - predictive importance of features, will be OC, OP, and distance in our case
 






