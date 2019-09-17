#Scripts to extract and view model results

## Folders

### `functions_for_R_package` - lists of functions for predicting TAD boundaries at bp resolution to be included in an R package

   * extract_boundaries_func - extract unique TAD boundaries from domain data created from a TAD caller
   * distance_predictor_types_func - calculating the distance from particular genomic annotation to genomic bins
   * percent_predictor_types_func - calculating percentage of overlap between regions defined by genomic annotations and genomic bins
   * count_predictor_types_func - calculating the total number of regions defined by genomic annotations that are inside genomic bins
   * binary_predictor_types_func - creating indicator variable for whether or not regions defined by genomic annotations overlap with genomic bins
   * annots_to_granges_func - converting a folder of .bed files defining genomic annotations into a GenomicRangesList object
   * annots_to_datamatrix_func - converting annotations in GenomicRangesList form into a datamatrix (with choice of predictor type) to be used as the feature space for downstream models
   * predict_tad_bounds_func - performing cell line and chromosome specific random forest classification model
   * predict_at_bp_resolution_func - predicting TAD boundaries using model derived from predict_tad_bounds_func at base pair resolution
   * all_function - all combined functions
   
### `function_for_tad_prediction` - different functions to run predictive models to predict TAD boundaries on the linear genome (outdated)

## Results Scripts


   