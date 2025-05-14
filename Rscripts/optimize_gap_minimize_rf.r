#!/usr/bin/env Rscript

# --- Bayesian Optimization to Minimize RF Distance via Gap Threshold ---
#
# Description:
# This script uses Bayesian Optimization to find the optimal gap threshold for
# sequence filtering. The objective is to minimize the Robinson-Foulds (RF)
# distance between a UPGMA tree built from the filtered subset (with singleton
# rescue) and a reference tree.
#
# Usage:
# Rscript optimize_gap_minimize_rf.r <bayes_output_csv> <original_dist_matrix_rds> \
#                                    <reference_tree_file> \
#                                    <gap_min> <gap_max> <gap_init_pts> <gap_n_iter> \
#                                    <output_prefix>
#
# Arguments:
#   bayes_output_csv         : Path to CSV from a previous run (e.g., bayes_kmedoids.r)
#                                containing original cluster assignments and sequence info.
#                                Needs columns: tip_label, cluster, specimen_id, sequence.
#   original_dist_matrix_rds : Path to the original distance matrix RDS file.
#   reference_tree_file      : Path to the reference tree file (Newick).
#   gap_min                  : Minimum gap threshold (%) to search (0-100).
#   gap_max                  : Maximum gap threshold (%) to search (0-100).
#   gap_init_pts             : Number of initial random points (gap values) for BO.
#   gap_n_iter               : Number of Bayesian Optimization iterations for gap.
#   output_prefix            : Prefix for output files (e.g., "optimized_gap_minRF").
#
# Output Files:
#   <output_prefix>_bayesopt_results.rds : BO results object.
#   <output_prefix>_bayesopt_summary.csv : BO summary table (Score = 1/(1+RF)).
#   <output_prefix>_bayesopt_plots.png   : BO diagnostic plots.
#   <output_prefix>_best_dist_matrix.rds : Subsetted distance matrix for the best gap threshold.
#   <output_prefix>_best_info.csv        : Info CSV for the best subset.
#   <output_prefix>_best_hclust.tre      : Hclust tree for the best subset.
#   <output_prefix>_log.txt              : Log file detailing the process.
#
# Example:
# Rscript optimize_gap_minimize_rf.r kmedoids_run_clusters.csv dist_obj.rds ref.tre \
#                                    10 90 10 30 \
#                                    optimized_gap_minRF
#

# --- Dependencies ---
cat("Loading required packages...\n")
suppressPackageStartupMessages({
    library(ParBayesianOptimization)
    library(ape)        # Tree functions, RF.dist
    library(dplyr)      # Data manipulation
    library(readr)      # read_csv
    library(argparse)   # Command-line arguments
    library(stats)      # hclust, as.dist
    library(purrr)      # map functions
    library(tibble)     # tibble
    library(phangorn)   # RF.dist alternative
    library(Biostrings) # Sequence handling (AAStringSet)
    library(stringr)    # str_count
    library(ggplot2)    # Plotting
    library(patchwork)  # Plotting
    library(ggtree)     # For tree visualization
    library(scales)     # For hue_pal
})

# --- Argument Parsing ---
parser <- ArgumentParser(description="Optimize gap threshold to minimize RF distance to reference, with singleton rescue.")
parser$add_argument("bayes_output_csv", help="Path to the input CSV with cluster assignments and sequence info.")
parser$add_argument("original_dist_matrix_rds", help="Path to the original input distance matrix RDS file.")
parser$add_argument("reference_tree_file", help="Path to the reference tree file.")
parser$add_argument("gap_min", type="double", help="Minimum gap threshold (%) (0-100).")
parser$add_argument("gap_max", type="double", help="Maximum gap threshold (%) (0-100).")
parser$add_argument("gap_init_pts", type="integer", help="Number of initial random points for BO.")
parser$add_argument("gap_n_iter", type="integer", help="Number of Bayesian Optimization iterations.")
parser$add_argument("output_prefix", help="Prefix for output files.")
parser$add_argument("num_cores", nargs='?', type="integer", default=1, help="[Optional] Number of CPU cores for parallel BO evaluation (default: 1).")

# Check if enough arguments are provided
min_mandatory_args <- 8
if (length(commandArgs(trailingOnly = TRUE)) < min_mandatory_args) {
    parser$print_help()
    stop("Error: Not enough required arguments provided.", call. = FALSE)
}

args <- parser$parse_args()

# Validate arguments
if (args$gap_min < 0 || args$gap_max > 100 || args$gap_min >= args$gap_max) {
    stop("Error: Invalid gap_min/gap_max range. Must be 0 <= gap_min < gap_max <= 100.")
}
if (args$gap_init_pts < 1) stop("Error: gap_init_pts must be >= 1.")
if (args$gap_n_iter < 1) stop("Error: gap_n_iter must be >= 1.")
if (!file.exists(args$bayes_output_csv)) stop("Error: Input CSV file not found: ", args$bayes_output_csv)
if (!file.exists(args$original_dist_matrix_rds)) stop("Error: Original distance matrix RDS file not found: ", args$original_dist_matrix_rds)
if (!file.exists(args$reference_tree_file)) stop("Error: Reference tree file not found: ", args$reference_tree_file)

# --- Setup Output Redirection ---
log_file <- paste0(args$output_prefix, "_log.txt")
cat("Redirecting console output to:", log_file, "\n")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
    cat("Closing log file...\n")
    sink(type = "message")
    sink(type = "output")
    close(log_con)
}, add = TRUE)

# Re-print configuration to the log file
cat("--- Configuration ---\n")
cat("Input Cluster CSV       :", args$bayes_output_csv, "\n")
cat("Input Original Dist RDS :", args$original_dist_matrix_rds, "\n")
cat("Reference Tree File     :", args$reference_tree_file, "\n")
cat("Gap Search Range        : [", args$gap_min, ",", args$gap_max, "]\n")
cat("Gap BO Init Points      :", args$gap_init_pts, "\n")
cat("Gap BO Iterations       :", args$gap_n_iter, "\n")
cat("Output Prefix           :", args$output_prefix, "\n")
cat("Num Cores Requested     :", args$num_cores, "\n")
cat("---------------------\n")

# --- Setup Parallel Backend ---
num_cores <- args$num_cores
if (is.na(num_cores) || num_cores < 1) {
    num_cores <- 1
    cat("Invalid num_cores provided, defaulting to 1 (serial execution).\n")
} else if (num_cores > parallel::detectCores(logical = FALSE)) {
    warning("Requested more cores than physically available. Using maximum available physical cores.")
    num_cores <- parallel::detectCores(logical = FALSE)
}
cat("Using", num_cores, "core(s) for optimization evaluation.\n")

cl <- NULL
if (num_cores > 1) {
    cat("Setting up parallel cluster...\n")
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)
    # Ensure cluster is stopped even if script errors. Add to existing on.exit for log_con.
    on.exit({ if(!is.null(cl)) { cat("Stopping parallel cluster...\n"); parallel::stopCluster(cl); doParallel::registerDoSEQ(); cl <- NULL }}, add = TRUE)
    cat("Parallel cluster setup complete.\n")
} else {
    cat("Running in serial mode.\n")
}

# --- Load Base Data (Load once outside the function) ---
cat("Loading base data...\n")
# Cluster assignments and sequence info from CSV
cluster_df_raw <- tryCatch(read_csv(args$bayes_output_csv, col_types = cols(.default = "c")),
                     error = function(e) stop("Error reading input CSV: ", e$message))
required_cols <- c("cluster", "tip_label", "specimen_id", "sequence")
missing_cols <- setdiff(required_cols, names(cluster_df_raw))
if (length(missing_cols) > 0) {
    stop("Error: Input CSV is missing required columns: ", paste(missing_cols, collapse=", "))
}
# Pre-process
cluster_df_raw <- cluster_df_raw %>%
    mutate(cluster = as.factor(cluster),
           sequence = as.character(sequence))
if(!is.character(cluster_df_raw$sequence)) {
    stop("Failed to convert 'sequence' column to character type.")
}
cluster_df_raw <- cluster_df_raw %>%
    filter(!is.na(sequence) & nchar(sequence) > 0) %>%
    mutate(seq_length = nchar(sequence))
cat("Loaded and pre-processed", nrow(cluster_df_raw), "rows with sequence data from CSV.\n")

# Original Distance Matrix
original_dist_object <- tryCatch(readRDS(args$original_dist_matrix_rds), error = function(e) stop("Error reading original distance matrix RDS: ", e$message))
# ... (rest of distance matrix loading and validation) ...
if (inherits(original_dist_object, "dist")) {
    original_dist_labels <- labels(original_dist_object)
    original_dist_matrix <- as.matrix(original_dist_object)
    cat("Loaded original dist object and converted to matrix.\n")
} else if (is.matrix(original_dist_object)) {
    original_dist_matrix <- original_dist_object
    original_dist_labels <- rownames(original_dist_matrix)
    if (!identical(rownames(original_dist_matrix), colnames(original_dist_matrix))) {
        warning("Row/col names differ in original matrix, attempting reconciliation.")
        common_names <- intersect(rownames(original_dist_matrix), colnames(original_dist_matrix))
        if(length(common_names) == 0) stop("Cannot reconcile row/col names in original matrix.")
        original_dist_matrix <- original_dist_matrix[common_names, common_names]
        original_dist_labels <- common_names
    }
     cat("Loaded original distance matrix.\n")
} else {
    stop("Original distance matrix RDS must contain 'dist' or 'matrix' object.")
}
if (is.null(original_dist_labels) || any(is.na(original_dist_labels) | original_dist_labels == "")) {
    stop("Invalid labels (missing or empty) in original distance matrix.")
}
cat("Found", length(original_dist_labels), "labels in original distance matrix.\n")

# Reference Tree
cat("Loading reference tree:", args$reference_tree_file, "\n")
ref_tree <- tryCatch(read.tree(args$reference_tree_file), error = function(e) {
    stop(paste("Error loading reference tree with read.tree:", e$message), call. = FALSE)
})
if (!inherits(ref_tree, "phylo")) {
     stop("Reference file did not load as a valid phylo object.", call. = FALSE)
}
ref_tree_tips <- ref_tree$tip.label
cat("Reference tree initially loaded with", length(ref_tree_tips), "tips.\n")

# Check for NAs introduced during loading and attempt pruning
if (any(is.na(ref_tree_tips))) {
    warning("NA values found in reference tree tip labels after initial loading. This might indicate issues with the tree file or label parsing.", call. = FALSE)
    cat("Number of NA tip labels found:", sum(is.na(ref_tree_tips)), "\n")

    valid_ref_tips <- ref_tree_tips[!is.na(ref_tree_tips)]
    num_valid_tips <- length(valid_ref_tips)
    num_na_tips <- length(ref_tree_tips) - num_valid_tips

    if (num_na_tips > 0) {
        cat("Attempting to prune reference tree to remove", num_na_tips, "tips with NA labels...\n")
        if (num_valid_tips < 3) {
             stop(paste("Fewer than 3 valid (non-NA) tips remain (", num_valid_tips, "). Cannot proceed."), call. = FALSE)
        }

        ref_tree_pruned <- tryCatch({
             keep.tip(ref_tree, valid_ref_tips)
           }, error = function(e) {
            cat("  ERROR pruning reference tree with keep.tip:", e$message, "\n")
            NULL # Indicate failure
        })

        if (!is.null(ref_tree_pruned) && inherits(ref_tree_pruned, "phylo")) {
            ref_tree <- ref_tree_pruned # Update tree object only if pruning succeeded
            ref_tree_tips <- ref_tree$tip.label # Update tip list
            cat("  Successfully pruned NA tips. Reference tree now has", length(ref_tree_tips), "tips.\n")
            # Re-check tip count just in case keep.tip did something unexpected
             if (length(ref_tree_tips) < 3) {
                 stop("Reference tree has fewer than 3 tips after successful pruning. Cannot proceed.", call. = FALSE)
            }
        } else {
            cat("  WARNING: Failed to prune NA tips from reference tree or result was invalid. Continuing with the original tree containing NAs, but downstream errors are likely.\n")
            # Keep the original ref_tree and ref_tree_tips which contain NAs
        }
    } else {
         # This case should not be reached if any(is.na(ref_tree_tips)) was TRUE, but included for completeness
         cat("  No NA tip labels found to prune.\n")
    }
} else {
     cat("  No NA tip labels found in reference tree.\n")
}

# --- Generate Hclust from Full Original Distance Matrix (for plotting reference) ---
cat("Generating Hclust tree from full original distance matrix for plotting reference...\n")
full_hclust_ref_tree <- NULL
tryCatch({
    if (!inherits(original_dist_matrix, "dist")) { # Ensure it's a dist object or convert
        if(!isSymmetric(original_dist_matrix)) {
            original_dist_matrix[lower.tri(original_dist_matrix)] <- t(original_dist_matrix)[lower.tri(original_dist_matrix)]
        }
        full_dist_obj <- as.dist(original_dist_matrix)
    } else {
        full_dist_obj <- original_dist_matrix # If it was already a dist object (though script converts to matrix)
    }
    full_hclust_result <- hclust(full_dist_obj, method = "average")
    full_hclust_ref_tree <- as.phylo(full_hclust_result)
    cat("Successfully generated Hclust tree from full distance matrix.\n")
}, error = function(e) {
    cat("Error generating Hclust tree from full distance matrix:", e$message, "\n")
    # full_hclust_ref_tree remains NULL
})


# --- Bayesian Optimization Objective Function (Core Logic) ---
evaluate_gap_for_rf_core <- function(gap_thresh_cont, 
                                   cluster_df_raw_arg, 
                                   original_dist_matrix_arg, 
                                   original_dist_labels_arg, 
                                   ref_tree_arg, 
                                   ref_tree_tips_arg,
                                   gap_min_arg,
                                   gap_max_arg) {
  # Explicitly load packages required within the function for parallel execution
  suppressPackageStartupMessages({
    library(dplyr, quietly = TRUE)
    library(stringr, quietly = TRUE)
    library(ape, quietly = TRUE)
    library(phangorn, quietly = TRUE)
    library(stats, quietly = TRUE)
  })

  current_gap_threshold <- max(0, min(100, round(gap_thresh_cont, 2)))
  cat(sprintf("\n--- Evaluating Gap Threshold = %.2f (Continuous = %.4f) ---\n", current_gap_threshold, gap_thresh_cont))

  # Default score (0 = worst score for optimizer)
  current_score <- 0.0
  normalized_rf_dist <- NA_real_
  n_final_common_for_log <- NA_integer_ # For logging
  sequence_penalty_component_for_log <- NA_real_ # For logging

  # --- Filtering and Rescue ---
  cat("Calculating gaps and filtering...\n")
  eval_df <- cluster_df_raw_arg %>% # Use argument
    mutate(
      gap_count = str_count(sequence, fixed("-")),
      gap_perc = ifelse(seq_length > 0, (gap_count / seq_length) * 100, 100)
    )
  filtered_df <- eval_df %>% filter(gap_perc <= current_gap_threshold)
  cat(nrow(filtered_df), "sequences passed initial filter.\n")

  original_clusters <- unique(eval_df$cluster)
  clusters_after_filter <- unique(filtered_df$cluster)
  lost_clusters <- setdiff(original_clusters, clusters_after_filter)
  cat(length(lost_clusters), "clusters lost after filtering.\n")

  rescued_rows_list <- list()
  if (length(lost_clusters) > 0) {
      cat("Rescuing singletons...\n")
      for (lost_cid in lost_clusters) {
          cluster_members <- eval_df %>% filter(cluster == lost_cid)
          if (nrow(cluster_members) > 0) {
              best_singleton <- cluster_members %>% arrange(gap_perc)
              best_singleton <- best_singleton[1, , drop = FALSE]
              rescued_rows_list[[as.character(lost_cid)]] <- best_singleton
          }
      }
  }
  rescued_df <- bind_rows(rescued_rows_list)
  final_kept_df <- bind_rows(filtered_df, rescued_df) %>%
                   distinct(specimen_id, .keep_all = TRUE)
  n_final_kept <- nrow(final_kept_df)
  cat("Total sequences kept after rescue:", n_final_kept, "\n")

  current_final_labels <- final_kept_df$tip_label # Use a local variable name
  current_final_labels <- intersect(current_final_labels, original_dist_labels_arg) # Use argument
  n_final_common <- length(current_final_labels)
  n_final_common_for_log <- n_final_common # Store for logging
  cat("Labels present in original distance matrix:", n_final_common, "\n")

  # Define failure return list structure for BO, including new fields
  failure_return <- list(Score = 0.0, Normalized_RF_Distance = NA_real_, Num_Sequences_Kept = n_final_common_for_log, Sequence_Penalty = NA_real_)


  if (n_final_common < 3) {
      cat("Fewer than 3 sequences remaining. Score = 0.\n")
      # Update failure_return with the actual (low) n_final_common
      failure_return$Num_Sequences_Kept <- n_final_common
      return(failure_return)
  }

  # --- Subset Data ---
  cat("Subsetting distance matrix...\n")
  subset_dist_matrix <- original_dist_matrix_arg[current_final_labels, current_final_labels] # Use argument and updated labels

  # --- Calculate RF Distance Score ---
  cat("Calculating RF distance score...\n")
  hclust_tree <- NULL
  tryCatch({
      if (!inherits(subset_dist_matrix, "dist")) {
          if(!isSymmetric(subset_dist_matrix)) {
               subset_dist_matrix[lower.tri(subset_dist_matrix)] <- t(subset_dist_matrix)[lower.tri(subset_dist_matrix)]
          }
          subset_dist_obj <- as.dist(subset_dist_matrix)
      } else {
          subset_dist_obj <- subset_dist_matrix
      }
      hclust_result <- hclust(subset_dist_obj, method = "average")
      hclust_tree <- as.phylo(hclust_result)
  }, error = function(e) {
      cat("  Hclust failed:", e$message, "\n")
  })

  if (!is.null(hclust_tree)) {
      common_tips_rf <- intersect(hclust_tree$tip.label, ref_tree_tips_arg) # Use argument
      n_common_rf <- length(common_tips_rf)
      if (n_common_rf >= 3) {
          pruned_hclust_tree <- tryCatch(keep.tip(hclust_tree, common_tips_rf), error=function(e) { cat("    Error pruning hclust tree:", e$message, "\n"); NULL })
          pruned_ref_tree <- tryCatch(keep.tip(ref_tree_arg, common_tips_rf), error=function(e) { cat("    Error pruning ref tree:", e$message, "\n"); NULL }) # Use argument

          if (!is.null(pruned_hclust_tree) && !is.null(pruned_ref_tree)) {
              cat("    Pruning successful. Forcing binary trees...\n")
              # Force binary trees to potentially avoid issues with RF.dist
              pruned_hclust_tree_binary <- tryCatch(multi2di(pruned_hclust_tree), error = function(e) { cat("    Error forcing binary hclust tree:", e$message, "\n"); NULL })
              pruned_ref_tree_binary <- tryCatch(multi2di(pruned_ref_tree), error = function(e) { cat("    Error forcing binary ref tree:", e$message, "\n"); NULL })

              if (!is.null(pruned_hclust_tree_binary) && !is.null(pruned_ref_tree_binary)) {
                  cat("    Binary trees created. Calculating Normalized RF distance...\n")
                  normalized_rf_dist <- tryCatch({
                      phangorn::RF.dist(pruned_hclust_tree_binary, pruned_ref_tree_binary, normalize=TRUE, check.labels=TRUE)
                  }, error=function(e) {
                      cat("    Error during phangorn::RF.dist (normalized):", e$message, "\n")
                      NA_real_
                  })
                  cat("    Normalized RF distance calculation finished (Result:", normalized_rf_dist, ").\n")
              } else { cat("    Failed to create binary trees for RF calculation.\n") }
          } else { cat("    Pruning for RF failed.\n") }
      } else { cat("    Fewer than 3 common tips for RF.\n") }
  } else { cat("    Hclust tree is NULL.\n") }

  # Score for maximization: 1 - normalized_rf (0 is best for norm_rf, so 1-0=1 is best for score)
  rf_score_component <- ifelse(is.na(normalized_rf_dist), 0.0, 1.0 - normalized_rf_dist)
  
  # Penalty for number of sequences
  sequence_count_penalty_weight <- 0.1 # Adjustable weight
  max_possible_sequences <- length(unique(original_dist_labels_arg))
  scaled_sequence_count <- ifelse(max_possible_sequences > 0, n_final_common / max_possible_sequences, 0)
  sequence_penalty_component <- sequence_count_penalty_weight * scaled_sequence_count
  sequence_penalty_component_for_log <- sequence_penalty_component # Store for logging

  current_score <- rf_score_component - sequence_penalty_component

  cat("  Normalized RF Distance:", normalized_rf_dist, "\n")
  cat("  RF Score Component (1 - NormRF):", rf_score_component, "\n")
  cat("  Num Sequences Kept:", n_final_common, " (Scaled:", scaled_sequence_count, ")\n")
  cat("  Sequence Count Penalty (weight=", sequence_count_penalty_weight, "):", sequence_penalty_component, "\n")
  cat("  Final Score for BO (RF_component - Seq_penalty):", current_score, "\n")

  # --- Combine Scores (only RF here) ---
  final_score <- current_score # The score to maximize

  # Add penalty for boundary values
  if (current_gap_threshold == gap_min_arg || current_gap_threshold == gap_max_arg) { # Use arguments
      final_score <- final_score - 1e-6
  }

  cat(sprintf("Gap = %.2f -> Final Score (Maximized) = %.6f\n", current_gap_threshold, final_score))

  # Return list expected by bayesOpt
  return(list(Score = final_score, 
              Normalized_RF_Distance = normalized_rf_dist,
              Num_Sequences_Kept = n_final_common_for_log,
              Sequence_Penalty = sequence_penalty_component_for_log))
}

# --- Wrapper Function for Bayesian Optimization ---
create_objective_wrapper_for_rf <- function(cluster_df_raw_arg, 
                                          original_dist_matrix_arg, 
                                          original_dist_labels_arg, 
                                          ref_tree_arg, 
                                          ref_tree_tips_arg,
                                          gap_min_arg,
                                          gap_max_arg) {
    force(cluster_df_raw_arg)
    force(original_dist_matrix_arg)
    force(original_dist_labels_arg)
    force(ref_tree_arg)
    force(ref_tree_tips_arg)
    force(gap_min_arg)
    force(gap_max_arg)

    wrapper <- function(gap_thresh_cont) {
        evaluate_gap_for_rf_core(
            gap_thresh_cont = gap_thresh_cont,
            cluster_df_raw_arg = cluster_df_raw_arg,
            original_dist_matrix_arg = original_dist_matrix_arg,
            original_dist_labels_arg = original_dist_labels_arg,
            ref_tree_arg = ref_tree_arg,
            ref_tree_tips_arg = ref_tree_tips_arg,
            gap_min_arg = gap_min_arg,
            gap_max_arg = gap_max_arg
        )
    }
    return(wrapper)
}

objective_function_for_bo_rf <- create_objective_wrapper_for_rf(
    cluster_df_raw_arg = cluster_df_raw,
    original_dist_matrix_arg = original_dist_matrix,
    original_dist_labels_arg = original_dist_labels,
    ref_tree_arg = ref_tree,
    ref_tree_tips_arg = ref_tree_tips,
    gap_min_arg = args$gap_min,
    gap_max_arg = args$gap_max
)

# --- Run Bayesian Optimization ---
cat("\n--- Starting Bayesian Optimization for Gap Threshold to Minimize RF ---\n")
gap_bounds <- list(gap_thresh_cont = c(args$gap_min, args$gap_max))

# Export to cluster if running in parallel
if (num_cores > 1 && !is.null(cl)) {
    cat("Exporting objective function and data to parallel cluster nodes...\n")
    parallel::clusterExport(cl, varlist=c("objective_function_for_bo_rf", "evaluate_gap_for_rf_core"), envir=environment())
    # Note: Data objects (cluster_df_raw, etc.) are in the closure of objective_function_for_bo_rf
    cat("Loading required libraries on cluster nodes...\n")
    parallel::clusterEvalQ(cl, { 
        suppressPackageStartupMessages({ 
            library(dplyr); library(stringr); library(ape); library(phangorn); library(stats); library(Biostrings)
        }) 
    })
    cat("Export and library loading complete for cluster.\n")
}

bo_args_list <- list(
    FUN = objective_function_for_bo_rf,
    bounds = gap_bounds,
    initPoints = args$gap_init_pts,
    iters.n = args$gap_n_iter,
    iters.k = ifelse(num_cores > 1, num_cores, 1), # Use num_cores for iters.k in parallel
    parallel = (num_cores > 1),
    verbose = 2
)

bayesOpt_results <- tryCatch({
    set.seed(222) # Different seed
    do.call(ParBayesianOptimization::bayesOpt, bo_args_list)
}, error = function(e) {
    cat("\nError during Bayesian Optimization:\n"); print(e); return(NULL)
})

cat("\n--- Bayesian Optimization Finished ---\n")
if (is.null(bayesOpt_results)) stop("Bayesian Optimization failed.", call. = FALSE)

# --- Extract and Print Initial Best Results ---
best_params_df <- getBestPars(bayesOpt_results)
best_gap_cont <- best_params_df$gap_thresh_cont
best_gap_rounded <- max(0, min(100, round(best_gap_cont, 2)))
best_bo_score <- max(bayesOpt_results$scoreSummary$Score) # Max score (1 - Normalized RF)

# Find the details corresponding to the best BO score
# Using head(1) instead of slice(1) to potentially avoid Rle list error
best_run_summary <- bayesOpt_results$scoreSummary %>% 
                    filter(Score == best_bo_score) %>% 
                    arrange(desc(Epoch)) %>% # In case of ties, take later epoch
                    head(1)

best_normalized_rf <- best_run_summary$Normalized_RF_Distance
best_num_sequences <- best_run_summary$Num_Sequences_Kept
best_sequence_penalty <- best_run_summary$Sequence_Penalty


cat("Optimal Gap Threshold (Rounded):", best_gap_rounded, "%\n")
cat("Best BO Score (1 - NormRF - SeqPenalty):", best_bo_score, "\n")
cat("  Corresponding Normalized RF Distance:", best_normalized_rf, "\n")
cat("  Corresponding Number of Sequences Kept:", best_num_sequences, "\n")
cat("  Corresponding Sequence Penalty Applied:", best_sequence_penalty, "\n")
# Note: RF per Tip % will be commented out for now.


# --- Save and Plot Optimization Results ---
cat("\nSaving optimization results...\n")
opt_rds_file <- paste0(args$output_prefix, "_bayesopt_results.rds")
saveRDS(bayesOpt_results, file = opt_rds_file)
cat("Optimization results object saved to:", opt_rds_file, "\n")

# Save summary table (includes raw RF distance)
summary_df <- bayesOpt_results$scoreSummary
summary_df$gap_threshold_rounded <- pmax(0, pmin(100, round(summary_df$gap_thresh_cont, 2)))
summary_csv_file <- paste0(args$output_prefix, "_bayesopt_summary.csv")
write.csv(summary_df, summary_csv_file, row.names = FALSE)
cat("Optimization summary saved to:", summary_csv_file, "\n")

# --- Plotting ---
if(nrow(summary_df) > 0) {
    cat("Generating optimization diagnostic plots...\n")
    p_progress <- NULL; p_explore <- NULL
    tryCatch({
        p_progress <- ggplot(summary_df, aes(x = Epoch, y = Score)) +
            geom_line(color = "red") + geom_point(color = "red") +
            labs(title = "A: BO Progress (Score = 1 - Normalized RF)", x = "Epoch", y = "Score (1 - Normalized RF)") + theme_bw()
        progress_plot_file <- paste0(args$output_prefix, "_bayesopt_progress_plot.png")
        ggsave(progress_plot_file, plot = p_progress, width = 8, height = 4)
        cat("Progress plot saved to:", progress_plot_file, "\n")
    }, error = function(e) { cat("Error generating progress plot:", e$message, "\n"); p_progress <<- NULL })
    tryCatch({
        plot_data <- summary_df %>% filter(is.finite(Score))
        best_plot_data <- data.frame(gap_thresh_cont=best_gap_cont, Score=best_bo_score) %>% filter(is.finite(Score)) # Use best_bo_score
        if (nrow(plot_data) > 0) {
             p_explore <- ggplot(plot_data, aes(x = gap_thresh_cont, y = Score)) +
                geom_point(alpha=0.6, size=2) +
                {if(nrow(best_plot_data) > 0) geom_point(data=best_plot_data, color='red', size=4, shape=17)} +
                labs(title = "B: Explored Gap Thresholds vs Score", x = "Gap Threshold (%)", y = "Score (1 - Normalized RF)") + theme_bw()
            explore_plot_file <- paste0(args$output_prefix, "_bayesopt_explore_plot.png")
            ggsave(explore_plot_file, plot = p_explore, width = 8, height = 4)
            cat("Exploration plot saved to:", explore_plot_file, "\n")
        } else { p_explore <<- NULL }
    }, error = function(e) { cat("Error generating exploration plot:", e$message, "\n"); p_explore <<- NULL })
    # No combined plot attempt
} else { cat("Optimization history is empty, cannot generate plots.\n") }


# --- Generate Final Outputs for Best Threshold ---
cat("\n--- Generating final outputs for best gap threshold:", best_gap_rounded, "% ---\n")

# Re-filter and rescue based on the best threshold found
cat("Re-filtering and rescuing singletons...\n")
final_eval_df <- cluster_df_raw %>%
    mutate(
      gap_count = str_count(sequence, fixed("-")),
      gap_perc = ifelse(seq_length > 0, (gap_count / seq_length) * 100, 100)
    )
final_filtered_df <- final_eval_df %>% filter(gap_perc <= best_gap_rounded)
final_original_clusters <- unique(final_eval_df$cluster)
final_clusters_after_filter <- unique(final_filtered_df$cluster)
final_lost_clusters <- setdiff(final_original_clusters, final_clusters_after_filter)
final_rescued_rows_list <- list()
if (length(final_lost_clusters) > 0) {
    for (lost_cid in final_lost_clusters) {
        cluster_members <- final_eval_df %>% filter(cluster == lost_cid)
        if (nrow(cluster_members) > 0) {
            best_singleton <- cluster_members %>% arrange(gap_perc)
            best_singleton <- best_singleton[1, , drop = FALSE]
            final_rescued_rows_list[[as.character(lost_cid)]] <- best_singleton
        }
    }
}
final_rescued_df <- bind_rows(final_rescued_rows_list)
final_kept_df <- bind_rows(final_filtered_df, final_rescued_df) %>%
                 distinct(specimen_id, .keep_all = TRUE)
final_labels <- intersect(final_kept_df$tip_label, original_dist_labels)
cat("Final number of sequences kept:", length(final_labels), "\n")

if (length(final_labels) >= 2) {
    # Subset final distance matrix
    final_subset_dist_matrix <- original_dist_matrix[final_labels, final_labels]
    dist_matrix_output_rds <- paste0(args$output_prefix, "_best_dist_matrix.rds")
    tryCatch({
        saveRDS(final_subset_dist_matrix, file = dist_matrix_output_rds)
        cat("Final subsetted distance matrix saved to:", dist_matrix_output_rds, "\n")
    }, error = function(e) { cat("Error saving final distance matrix:", e$message, "\n") })

    # Save final info CSV
    info_output_csv <- paste0(args$output_prefix, "_best_info.csv")
    tryCatch({
        output_info_df <- final_kept_df %>%
                          filter(tip_label %in% final_labels) %>%
                          select(tip_label, cluster, specimen_id, seq_length, gap_perc)
        write.csv(output_info_df, info_output_csv, row.names = FALSE, quote = TRUE, na = "")
        cat("Final sequence info CSV saved to:", info_output_csv, "\n")
    }, error = function(e) { cat("Error saving final info CSV:", e$message, "\n") })

    # Perform and save final Hclust tree
    cat("Performing final hierarchical clustering...\n")
    tryCatch({
        if (!inherits(final_subset_dist_matrix, "dist")) {
             if(!isSymmetric(final_subset_dist_matrix)) {
                 final_subset_dist_matrix[lower.tri(final_subset_dist_matrix)] <- t(final_subset_dist_matrix)[lower.tri(final_subset_dist_matrix)]
             }
             final_subset_dist_obj <- as.dist(final_subset_dist_matrix)
        } else { final_subset_dist_obj <- final_subset_dist_matrix }
        hclust_result <- hclust(final_subset_dist_obj, method = "average")
        hclust_tree <- as.phylo(hclust_result)

        # --- RF per Tip Percentage Calculation (Commented out as it was based on raw RF) ---
        # cat("Calculating RF distance per tip percentage for the final optimal tree...\n")
        # rf_per_tip_perc <- NA_real_ # Initialize
        # n_final_common_rf <- NA_integer_ # Initialize
        # if (!is.null(hclust_tree) && !is.na(best_normalized_rf)) { # Changed best_raw_rf to best_normalized_rf
        #     # This calculation would need to be re-thought if desired with normalized RF
        #     # For now, we'll just report the normalized RF directly.
        #     # final_common_tips_rf <- intersect(hclust_tree$tip.label, ref_tree_tips)
        #     # n_final_common_rf <- length(final_common_tips_rf)
        #     # cat("  Number of common tips between final tree and reference:", n_final_common_rf, "\n")
        #     # if (n_final_common_rf > 0) {
        #     #      # Example: To get raw RF back: raw_rf = normalized_rf * (2*n_tips - 6)
        #     #      # max_rf_possible = (2 * n_final_common_rf - 6) # for unrooted binary
        #     #      # estimated_raw_rf = best_normalized_rf * max_rf_possible
        #     #      # rf_per_tip_perc <- (estimated_raw_rf / n_final_common_rf) * 100
        #     #      # cat("  Estimated Final RF Distance per Tip (%):", rf_per_tip_perc, "%\n")
        #     #      if (n_final_common_rf < 3) {
        #     #          cat("  Warning: Fewer than 3 common tips, RF distance may be unreliable or NA.\n")
        #     #      }
        #     # } else {
        #     #     cat("  Zero common tips, cannot calculate RF per tip percentage.\n")
        #     # }
        # } else {
        #      cat("  Final Hclust tree or best Normalized RF distance is missing, cannot calculate percentage.\n")
        # }
        # --- End RF per Tip Percentage Calculation ---

        tree_output_file <- paste0(args$output_prefix, "_best_hclust.tre")
        write.tree(hclust_tree, file = tree_output_file)
        cat("Final hierarchical clustering tree saved to:", tree_output_file, "\n")

        # --- Calculate Final Patristic Distance Correlation ---
        cat("Calculating final Patristic Distance Correlation...\n")
        final_patristic_cor <- NA_real_
        if (!is.null(hclust_tree) && !is.null(ref_tree)) {
            common_tips_patristic <- intersect(hclust_tree$tip.label, ref_tree$tip.label)
            cat("  Number of common tips for patristic correlation:", length(common_tips_patristic), "\n")

            if (length(common_tips_patristic) >= 3) {
                pruned_hclust_pat <- tryCatch(keep.tip(hclust_tree, common_tips_patristic), error=function(e) NULL)
                pruned_ref_pat <- tryCatch(keep.tip(ref_tree, common_tips_patristic), error=function(e) NULL)

                if (!is.null(pruned_hclust_pat) && !is.null(pruned_ref_pat)) {
                    pruned_hclust_pat_binary <- tryCatch(multi2di(pruned_hclust_pat), error=function(e) NULL)
                    pruned_ref_pat_binary <- tryCatch(multi2di(pruned_ref_pat), error=function(e) NULL)

                    if (!is.null(pruned_hclust_pat_binary) && !is.null(pruned_ref_pat_binary)) {
                        pat_dist_hclust <- tryCatch(cophenetic.phylo(pruned_hclust_pat_binary), error=function(e) NULL)
                        pat_dist_ref <- tryCatch(cophenetic.phylo(pruned_ref_pat_binary), error=function(e) NULL)

                        if (!is.null(pat_dist_hclust) && !is.null(pat_dist_ref)) {
                            common_order_pat <- sort(rownames(pat_dist_hclust)) # Ensure consistent order
                            pat_dist_hclust_vec <- pat_dist_hclust[common_order_pat, common_order_pat][upper.tri(pat_dist_hclust[common_order_pat, common_order_pat])]
                            pat_dist_ref_vec <- pat_dist_ref[common_order_pat, common_order_pat][upper.tri(pat_dist_ref[common_order_pat, common_order_pat])]

                            if(length(pat_dist_hclust_vec) > 1 && length(pat_dist_ref_vec) > 1 &&
                               !all(is.na(pat_dist_hclust_vec)) && !all(is.na(pat_dist_ref_vec)) &&
                               (sd(pat_dist_hclust_vec, na.rm=TRUE) > 0 || sd(pat_dist_ref_vec, na.rm=TRUE) > 0) ) { # Check for non-zero variance if possible
                                final_patristic_cor <- tryCatch(cor(pat_dist_hclust_vec, pat_dist_ref_vec, method = "pearson", use = "complete.obs"), error = function(e) NA_real_)
                                cat("  Final Patristic Correlation calculated:", final_patristic_cor, "\n")
                            } else {
                                cat("  Not enough unique distance pairs or zero variance for correlation.\n")
                            }
                        } else { cat("  Failed to calculate patristic distance matrices for correlation.\n") }
                    } else { cat("  Failed to create binary trees for patristic correlation.\n") }
                } else { cat("  Pruning for patristic correlation failed.\n") }
            } else { cat("  Fewer than 3 common tips for patristic correlation.\n") }
        } else { cat("  Final Hclust tree or reference tree is missing, cannot calculate patristic correlation.\n") }
        # --- End Patristic Distance Correlation ---

    }, error = function(e) { cat("Error during final hierarchical clustering or saving tree:", e$message, "\n") })

} else {
    cat("Fewer than 2 sequences in final set, cannot generate distance matrix or tree.\n")
    rf_per_tip_perc <- NA_real_ # Ensure it's defined even if tree fails
}

cat("\n--- Script Finished ---\n")
# Final summary lines were printed above.
# Print the best normalized RF and num sequences again for clarity at the very end.
cat("\n--- Final Best Result Summary ---\n")
cat("Optimal Gap Threshold:", best_gap_rounded, "%\n")
if (!is.na(best_normalized_rf)) {
    cat("Best Minimum Normalized RF Distance achieved:", best_normalized_rf, "\n")
} else {
    cat("Best Minimum Normalized RF Distance: Not available (check log).\n")
}
if (!is.na(best_num_sequences)) {
    cat("Number of Sequences Kept at Best Threshold:", best_num_sequences, "\n")
} else {
    cat("Number of Sequences Kept at Best Threshold: Not available (check log).\n")
}

if (exists("final_patristic_cor") && !is.na(final_patristic_cor)) { # Print patristic correlation if calculated
    cat("Final Patristic Correlation:", final_patristic_cor, "\n")
} else {
    cat("Final Patristic Correlation: Not calculated (check log).\n")
}

# --- Generate Final Tree Plot ---
cat("\n--- Generating Final Tree Plot (Full Trees Side-by-Side) ---\n")

plot_left_tree_obj <- NULL
plot_right_tree_obj <- NULL
p_left_full <- NULL
p_right_full <- NULL

# Prepare left tree (Hclust of Full Dataset)
if (!is.null(full_hclust_ref_tree) && exists("cluster_df_raw")) {
    plot_left_tree_obj <- full_hclust_ref_tree
    cat("Preparing left panel: Hclust of Full Dataset (", length(plot_left_tree_obj$tip.label), " tips)\n")
    
    # Prepare annotation data, ensuring tip_label is present for joining, and specimen_id for display
    annotation_data_left <- cluster_df_raw %>%
        filter(tip_label %in% plot_left_tree_obj$tip.label) %>%
        select(tip_label, specimen_id, Cluster = cluster) %>% # Use tip_label for join, specimen_id for display
        filter(!is.na(Cluster) & Cluster != "") %>%
        mutate(Cluster = factor(Cluster)) %>%
        distinct(tip_label, .keep_all = TRUE) # Ensure one annotation row per tree tip

    p_left_full <- ggtree(plot_left_tree_obj, size = 0.2, layout = 'circular')
    if (nrow(annotation_data_left) > 0 && nlevels(annotation_data_left$Cluster) > 0) {
        p_left_full <- p_left_full %<+% annotation_data_left + aes(color = Cluster) + # Joins by matching tree$tip.label with first col of annotation_data_left (tip_label)
                       geom_tippoint(size = 0.1) + 
                       scale_color_manual(values = scales::hue_pal()(nlevels(annotation_data_left$Cluster)), guide = "none", drop = FALSE)
    }
    p_left_full <- p_left_full + 
                   geom_tiplab(aes(label = specimen_id), size = 0.3, offset = 0.001, align=TRUE, linesize=0.05, hjust = 0.5) + # Display specimen_id
                   ggtitle(paste0("Hclust of Full Dataset (", length(plot_left_tree_obj$tip.label), " tips)")) +
                   theme(legend.position = "none", plot.title = element_text(size=10))
} else {
    cat("Full Hclust reference tree or cluster data not available for left panel.\n")
}

# Prepare right tree (Optimized HCLUST Tree)
if (exists("hclust_tree") && !is.null(hclust_tree) && exists("cluster_df_raw")) {
    plot_right_tree_obj <- hclust_tree
    cat("Preparing right panel: Optimized HCLUST Tree (", length(plot_right_tree_obj$tip.label), " tips)\n")

    cat("Preparing right panel: Optimized HCLUST Tree (", length(plot_right_tree_obj$tip.label), " tips)\n")

    # Create unique cluster IDs for each tip in the optimized tree for unique coloring
    annotation_data_right <- tibble(
        tip_label = plot_right_tree_obj$tip.label, # Original full tip labels from the tree
        specimen_id = sub("-.+", "", plot_right_tree_obj$tip.label), # Extract specimen ID for display
        Cluster = factor(1:length(plot_right_tree_obj$tip.label)) # Assign a unique factor level to each tip
    )
    cat("  Generated unique cluster IDs for right panel plotting. N_levels:", nlevels(annotation_data_right$Cluster), "\n")

    p_right_full <- ggtree(plot_right_tree_obj, size = 0.5, layout = 'circular')
    # Color tip points and labels by the new unique Cluster factor
    if (nrow(annotation_data_right) > 0 && nlevels(annotation_data_right$Cluster) > 0) {
        p_right_full <- p_right_full %<+% annotation_data_right + 
                        aes(color = Cluster) + 
                        geom_tippoint(size = 1) +
                        scale_color_manual(values = scales::hue_pal()(nlevels(annotation_data_right$Cluster)), guide = "none", drop = FALSE) +
                        geom_tiplab(aes(label = specimen_id, color = Cluster), size = 1.5, offset = 0.02, align=TRUE, linesize=0.2, hjust = 0.5) # Display specimen_id, color by new Cluster
    } else {
         p_right_full <- p_right_full + 
                        geom_tiplab(aes(label = specimen_id), size = 1.5, offset = 0.02, align=TRUE, linesize=0.2, hjust = 0.5) # Display specimen_id if no cluster data
    }
    p_right_full <- p_right_full + 
                    ggtitle(paste0("Optimized HCLUST Tree (", length(plot_right_tree_obj$tip.label), " tips)")) +
                    theme(legend.position = "none", plot.title = element_text(size=10))
} else {
    cat("Optimized Hclust tree or cluster data not available for right panel.\n")
}

# Combine and save if both plots were generated
if (!is.null(p_left_full) && !is.null(p_right_full)) {
    combined_tree_plot <- p_left_full + p_right_full + plot_layout(ncol = 2)
    plot_output_file <- paste0(args$output_prefix, "_full_vs_optimized_hclust_plot.png")
    
    height_left_plot <- if(!is.null(plot_left_tree_obj)) length(plot_left_tree_obj$tip.label) * 0.03 else 8
    height_right_plot <- if(!is.null(plot_right_tree_obj)) length(plot_right_tree_obj$tip.label) * 0.06 else 8
    plot_final_height <- max(c(10, height_left_plot, height_right_plot)) 
    plot_final_height <- min(plot_final_height, 12) # Cap height at 120 inches

    tryCatch({
        ggsave(plot_output_file, plot = combined_tree_plot, width = 18, height = plot_final_height, units = "in", dpi = 300, limitsize = FALSE) # Increased width, adjusted dpi
        cat("Full vs Optimized HCLUST tree plot saved to:", plot_output_file, "\n")
    }, error = function(e) {
        cat("Error saving combined tree plot:", e$message, "\n")
    })
} else {
    cat("One or both trees could not be prepared for the full comparison plot.\n")
}
# --- End of Final Tree Plot ---

# Note: on.exit() will handle closing the sinks
