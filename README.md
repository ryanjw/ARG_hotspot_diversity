# ARG_hotspot_diversity
This is the repo for all files for the AMR paper.  Scripts were run interactively and contain many functions used to visualize data, generate summary statistics, and make decisions regarding further analyses and figures to generate.  The included dataset also inlcudes all MG-RAST ids used in this study.  Note that AMR and ARG (antibiotic resistance gene) are used interchangeably and as shorthand throughout these scripts.  Several scripts here may reflect analyses not included in the study but reflect additional information and data exploration.  Below are brief summaries for each script located in the scripts/ directory:

## row_column_correction.R
This contains a function for managing issues where headers and columns did not line up correctly when reading in data.

## column_splitter_func.R
This script contains a function used to split strings in a column by a common delimiter.

## libraries.R
This script contains all of the appropriate libraries used across these scripts.

## merging_data_w_metadata.R
This script was used to manage linking metadata with data.

## mixed_model_rsq.R
This script was used to generate summary statistics for mixed models.

## arg_hotspots.R
This script was used to find samples that belong to the 95th percentile based on AMR (ARG) abundance.  Figure 2 was generated within this script using the ggtern() function.  Due to issues with generating this figure, other legends were generated using ggplot2() and added into the image post analysis.   

## finding-reca-int-arg-agreements.R
This script was used to find co-occurrence of different genes in the study within genomes as part of the RefSoil database.  Similar analysis was run on HMP data (not included in this paper), and results are discussed in the text along with Supplemental Table 1.

## making_distances.R
This script was used to make distances between recA and AMR samples along with Figure 1.

## final-co-occurrence.R
This is the final co-occurrence script that was used to make Figure 3B (which was also edited for labeling post-analysis).  Note that other co-occurrence scripts do not produce figures here but may contain important exploratory analyses.

## picking_out_indicator_args.R
This script was used to pick out indicator ARGs and do produce Figure 3A.
