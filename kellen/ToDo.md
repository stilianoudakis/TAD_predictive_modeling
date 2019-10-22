# ToDo

- Train the model on resolution-level response vector (detected TAD boundaries)
- Predict on basepair-level data - get global maximum and convert to resolution-level response vector (predicted TAD boundaries)
- Goal - to get rid of multiple false positives. Instead, we expect detected and predicted response vectors be very similar (test with flanking). Consequently, MCC or F1 may show better, more distinguishable performance.

- Dismal performance of random forest TBD
    - People use Linear, Logistic, GBM, RF. See Ramírez et al., “High-Resolution TADs Reveal DNA Sequences Underlying Genome Organization in Flies”, Methods section, how to use

- Describe variable importance in Methods. After polishing penalized regression and/or RF.

- Need legends for tables

- Can we predict boundary STRENGTH? That is, response variable is continuous boundary score.

- Outcome: An R package to predict TAD boundaries at basepair resolution
    - Input: Coordinates of TAD boundaries (resolution-level), and the corresponding genomic annotations
    - Modeling: Using any method we tested, train a model with optimal parameters we identified, with cross-validation. 
    - Output: Predict basepair-level coordinates of TAD boundaries
    - Case scenario: K562 cell line
    