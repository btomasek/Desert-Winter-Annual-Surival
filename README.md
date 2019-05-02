# Desert-Winter-Annual-Survival

This repository is for archiving code used for survival data analysis of Sonoran Desert Winter Annual plants. This model is intended to simultaneously account for three primary issues common to many ecological survival data sets.

1. Allow for the interpretation of the effects of higher-frequency environmental data with lower-frequency survival data from periodic census-taking.

2. Account for uncertainty arising from interval censoring.

3. Allow for the partitioning of the effects of predictor variables on survival coming from acute (single day) versus chronic (multiple days) conditions. Chronic conditions effect survival through an autoregressive structure in the model, i.e. AR(1). 

The underlying dependent variable for this analysis is binary. This binary condition is assessed at discrete time points which may lead to interval censoring. Predictor variables can categorical or continuous and may be measured at a higher frequency than the dependent variable.

Within subfolders, there is data and example code of how to use custom functions can replicate the analysis and visualizations. We are continuing to improve the efficiency and user-friendliness of the code.