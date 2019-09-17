# Scripts to extract and view model results

## Folders

### `functions_for_R_package` - lists of functions for predicting TAD boundaries at bp resolution to be included in an R package

   * **extract_boundaries_func** - extract unique TAD boundaries from domain data created from a TAD caller
   * **distance_predictor_types_func** - calculating the distance from particular genomic annotation to genomic bins
   * **percent_predictor_types_func** - calculating percentage of overlap between regions defined by genomic annotations and genomic bins
   * **count_predictor_types_func** - calculating the total number of regions defined by genomic annotations that are inside genomic bins
   * **binary_predictor_types_func** - creating indicator variable for whether or not regions defined by genomic annotations overlap with genomic bins
   * **annots_to_granges_func** - converting a folder of .bed files defining genomic annotations into a GenomicRangesList object
   * **annots_to_datamatrix_func** - converting annotations in GenomicRangesList form into a datamatrix (with choice of predictor type) to be used as the feature space for downstream models
   * **predict_tad_bounds_func** - performing cell line and chromosome specific random forest classification model
   * predict_at_bp_resolution_func - predicting TAD boundaries using model derived from predict_tad_bounds_func at base pair resolution
   * **all_function** - all combined functions
   
### `function_for_tad_prediction` - different functions to run predictive models to predict TAD boundaries on the linear genome (outdated)

## Results Scripts

   * **evaluating_model_performances_rf** - plotting average model performance (MCC) across chromosomes of all random forest classifiers in the ensemble framework
   * **comparing_predicted_boundaries_vs_called_boundaries** - plotting enrichment heatmaps of signal strengths of top genomic annotations (CTCF, RAD21, SMC3, ZNF143) around 10 kb flank predicted TAD boundaries and called TAD boundaries
   * **visualizing_example_predicted_tad_bounds** - visualizing an example of clustered predicted TAD boundary regions on the 17390000-18640000 section of chr22 (first called TAD on chr 22; 2 boundaries)
   * **evaluating_jaccard_indices** - visualizing the pairwise jaccard indices of genomic annotations using heatmaps
   * **proportion_of_importance_per_annotation** - plotting the proportion of variable importance specific to sets of genomic annotations (histone modifications, chromhmm, and transcription factor binding sites) from random forest classifiers
   * **evaluating_rfe_accross_chromosomes** - plotting where performances of random forest models stabilized as the number of genomic annotations included increases by powers of 2
   * **evaluating_preprocessed_tad_summaries** - evaluating how many TADs were removed due to exceeding 2 mb in width, and the subsequent number of unique TAD boundaries extracted
   * **evaluating_imbalance_ratios** - summaries of class imbalance created from modeling framework
   * **comparing_tad_regions_across_cl_using_vd** - venn diagrams comparing the overlap of TAD boundary regions between GM12878 and K562
   * **comparing_tfbs_specific_models_with_all_annotations** - comparing annotation (histone modifications, chromhmm, tfbs) specific models
   
   