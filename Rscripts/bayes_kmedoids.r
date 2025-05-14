#!/usr/bin/env Rscript

# --- Bayesian Optimization for Optimal K-Medoids Clusters (k) ---
#
# Description:
# This script uses Bayesian Optimization to find an optimal number of clusters (k)
# for k-medoids (PAM) clustering based on multiple criteria:
#   - Median cluster size within a target range.
#   - Mean cluster size close to a target value.
#   - High average within-cluster sequence identity.
#   - High Silhouette Score based on PATRISTIC distances from an original tree.
#
# It requires:
#   1. A pre-computed distance matrix (RDS format) for the k-medoids clustering.
#      Labels (rownames/colnames) should follow 'SpecimenID-Rank1-Rank2...'.
#   2. The original phylogenetic tree file (e.g., Newick) for calculating patristic
#      distances used in the Silhouette score evaluation. Tree tip labels must match
#      the labels in the distance matrix.
#   3. A corresponding FASTA file (AA sequences) where headers match the *first part*
#      (SpecimenID, before the first hyphen) of the distance matrix/tree labels.
#      The script will internally filter this FASTA to match the distance matrix labels.
#
# Usage:
# Rscript bayes_kmedoids.r <dist_matrix_rds> <tree_file> <fasta_file> <output_prefix> \
#                          <k_min> <k_max> <init_pts> <n_iter> \
#                          <med_min> <med_max> <mean_target> <min_msa_size> \
#                          <w_median> <w_mean> <w_identity> <w_silhouette> [num_cores]
#
# Arguments:
#   dist_matrix_rds    : Path to the distance matrix object (.rds file) for PAM clustering.
#   tree_file          : Path to the phylogenetic tree file (e.g., .tre, .newick) for Silhouette evaluation.
#   fasta_file         : Path to the sequence file (.fasta, .fa) - AA sequences.
#   output_prefix      : Prefix for output files (e.g., "kmedoids_run").
#   k_min              : Minimum number of clusters (k) to search (integer >= 2).
#   k_max              : Maximum number of clusters (k) to search (integer).
#   init_pts           : Number of initial random points (k values) for BO (integer).
#   n_iter             : Number of Bayesian Optimization iterations (integer).
#   med_min            : Target minimum median cluster size (numeric).
#   med_max            : Target maximum median cluster size (numeric).
#   mean_target        : Target mean cluster size (numeric).
#   min_msa_size       : Minimum cluster size for MSA/identity calculation (integer).
#   w_median           : Weight for median size score (numeric).
#   w_mean             : Weight for mean size score (numeric).
#   w_identity         : Weight for sequence identity score (numeric).
#   w_silhouette       : Weight for patristic silhouette score (numeric).
#   num_cores          : [Optional] Number of CPU cores for parallel processing (integer).
#
# Example:
# Rscript bayes_kmedoids.r seq_dist_matrix.rds tree.tre sequences.fasta kmedoids_run \
#                          5 100 20 50 5 12 10 3 0.3 0.15 0.45 0.1 8
#

# --- Dependencies ---
cat("Loading required packages...\n")
suppressPackageStartupMessages({
    library(ParBayesianOptimization)
    library(doParallel)
    library(DECIPHER) # For AlignSeqs, DistanceMatrix, RemoveGaps
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(ggplot2)
    library(ape)      # For read.tree, cophenetic.phylo
    library(Biostrings)
    library(parallel)
    library(cluster)  # For pam, silhouette
    library(pbapply)
    library(patchwork)
    library(viridis)
})

# --- Argument Parsing ---
args <- commandArgs(trailingOnly = TRUE)
script_start_time <- Sys.time() # VERBOSE
cat(sprintf("bayes_kmedoids.r started at: %s\n", script_start_time)) # VERBOSE

# Expected number of mandatory arguments
min_args <- 16 # Removed gap_threshold
max_args <- 17 # Removed gap_threshold, added optional num_cores

if (length(args) < min_args || length(args) > max_args) {
  cat("Usage: Rscript bayes_kmedoids.r <dist_matrix_rds> <tree_file> <fasta_file> <output_prefix> \\\n")
  cat("                                <k_min> <k_max> <init_pts> <n_iter> \\\n")
  cat("                                <med_min> <med_max> <mean_target> <min_msa_size> \\\n")
  cat("                                <w_median> <w_mean> <w_identity> <w_silhouette> [num_cores]\n\n")
  cat("Error: Incorrect number of arguments provided.\n")
  cat("Received", length(args), "arguments.\n")
  quit(status = 1)
}

# Assign arguments to variables
dist_matrix_rds       <- args[1]
tree_file             <- args[2]
sequence_fasta_file   <- args[3]
output_prefix         <- args[4]
k_min                 <- as.integer(args[5])
k_max                 <- as.integer(args[6])
bayesopt_init_points  <- as.integer(args[7])
bayesopt_n_iter       <- as.integer(args[8])
target_median_min     <- as.numeric(args[9])
target_median_max     <- as.numeric(args[10])
target_mean           <- as.numeric(args[11])
min_cluster_size_for_msa <- as.integer(args[12])
weight_median         <- as.numeric(args[13])
weight_mean           <- as.numeric(args[14])
weight_identity       <- as.numeric(args[15])
weight_silhouette     <- as.numeric(args[16])
# gap_threshold removed

# Handle optional num_cores argument
num_cores_arg <- NA
if (length(args) == max_args) {
  num_cores_arg <- as.integer(args[17]) # Index updated
}

# Validate numeric/integer arguments
if (k_min < 2) stop("Error: k_min must be >= 2.")
if (k_max < k_min) stop("Error: k_max must be >= k_min.")
# gap_threshold validation removed


numeric_args <- list(k_min=k_min, k_max=k_max, bayesopt_init_points=bayesopt_init_points,
                     bayesopt_n_iter=bayesopt_n_iter, target_median_min=target_median_min,
                     target_median_max=target_median_max, target_mean=target_mean,
                     min_cluster_size_for_msa=min_cluster_size_for_msa, weight_median=weight_median,
                     weight_mean=weight_mean, weight_identity=weight_identity,
                     weight_silhouette=weight_silhouette, num_cores_arg=num_cores_arg) # Removed gap_threshold

for (name in names(numeric_args)) {
  val <- numeric_args[[name]]
  if (name == "num_cores_arg" && is.na(val)) next
  if (is.null(val) || is.na(val) || (!is.numeric(val) && !is.integer(val)) ) {
    arg_index <- which(names(numeric_args) == name)
    original_arg_pos <- c(5:16, 17)[arg_index] # Approximate mapping - indices updated
    stop("Error: Argument '", name, "' must be a valid number/integer. Received: '", args[original_arg_pos], "'")
  }
}

# Validate file existence
if (!file.exists(dist_matrix_rds)) stop("Error: Distance matrix RDS file not found: ", dist_matrix_rds)
if (!file.exists(tree_file)) stop("Error: Tree file not found: ", tree_file)
if (!file.exists(sequence_fasta_file)) stop("Error: Sequence FASTA file not found: ", sequence_fasta_file)

# Set up parallel processing cores
if (!is.na(num_cores_arg) && num_cores_arg > 0) {
    num_cores <- num_cores_arg
} else {
    num_cores <- detectCores() - 1
    if (is.na(num_cores) || num_cores < 1) num_cores <- 1
}

cat("--- Configuration ---\n")
cat("Distance Matrix RDS  :", dist_matrix_rds, "\n")
cat("Tree File            :", tree_file, "\n")
cat("Sequence FASTA File  :", sequence_fasta_file, "\n")
cat("Output Prefix        :", output_prefix, "\n")
cat("k Search Range       : [", k_min, ",", k_max, "]\n")
cat("BO Init Points       :", bayesopt_init_points, "\n")
cat("BO Iterations        :", bayesopt_n_iter, "\n")
cat("Target Median Range  : [", target_median_min, ",", target_median_max, "]\n")
cat("Target Mean Size     :", target_mean, "\n")
cat("Min Cluster for MSA  :", min_cluster_size_for_msa, "\n")
cat("Weight Median        :", weight_median, "\n")
cat("Weight Mean          :", weight_mean, "\n")
cat("Weight Identity      :", weight_identity, "\n")
cat("Weight Silhouette    :", weight_silhouette, "\n")
# Gap Threshold print removed
cat("Using", num_cores, "cores for parallel processing.\n")
cat("---------------------\n")


# --- Load Data ---
cat("VERBOSE: Starting Data Loading section.\n") # VERBOSE
cat("Loading distance matrix from RDS:", dist_matrix_rds, "\n")
clustering_dist_object <- tryCatch({
    readRDS(dist_matrix_rds)
}, error = function(e) {
    stop("Error reading distance matrix RDS file '", dist_matrix_rds, "': ", e$message)
})

# Check if it's a dist object or matrix
if (!inherits(clustering_dist_object, "dist") && !is.matrix(clustering_dist_object)) {
    stop("Input distance matrix RDS must contain an object of class 'dist' or 'matrix'. Found: ", class(clustering_dist_object)[1])
}

# --- Get Labels and Ensure Matrix Format ---
if (inherits(clustering_dist_object, "dist")) {
    # Convert dist object to a full matrix for consistent handling
    dist_labels <- labels(clustering_dist_object)
    if (is.null(dist_labels)) stop("Distance object loaded from RDS must have labels.")
    # The pam function can take a dist object directly, but converting to matrix
    # makes subsetting and label checking consistent with matrix input.
    clustering_dist_matrix <- as.matrix(clustering_dist_object)
    cat("Converted dist object to matrix.\n")
} else { # It's already a matrix
    clustering_dist_matrix <- clustering_dist_object
    dist_labels <- rownames(clustering_dist_matrix)
    if (is.null(dist_labels)) stop("Distance matrix loaded from RDS must have row names.")
    if (!identical(rownames(clustering_dist_matrix), colnames(clustering_dist_matrix))) {
         warning("Row names and column names of the loaded distance matrix are not identical. Attempting reconciliation.")
         common_names <- intersect(rownames(clustering_dist_matrix), colnames(clustering_dist_matrix))
         if (length(common_names) > 0) {
             clustering_dist_matrix <- clustering_dist_matrix[common_names, common_names]
             dist_labels <- common_names # Update labels to common set
             cat("Subsetted matrix to common row/column names:", length(common_names), "\n")
         } else {
             stop("Cannot reconcile row/column names in distance matrix.")
         }
    }
}
# Validate labels AFTER potential subsetting
if (any(is.na(dist_labels) | dist_labels == "")) stop("Distance matrix contains missing/empty labels after processing.")
if (any(duplicated(dist_labels))) stop("Distance matrix contains duplicated labels after processing.")
cat("Total labels in distance matrix:", length(dist_labels), "\n")

cat("Loading phylogenetic tree from:", tree_file, "\n")
tic_tree_load <- Sys.time() # VERBOSE
phylo_tree <- tryCatch(ape::read.tree(tree_file), error = function(e) stop("Error loading tree: ", e$message))
tree_tip_labels <- phylo_tree$tip.label
toc_tree_load <- Sys.time() # VERBOSE
cat(sprintf("VERBOSE: Tree loading took %.2f seconds.\n", as.numeric(difftime(toc_tree_load, tic_tree_load, units="secs")))) # VERBOSE
cat("Total tips in tree:", length(tree_tip_labels), "\n")
if(length(tree_tip_labels) == 0) stop("Error: No tip labels found in the tree object.")

# --- Load and Filter Sequences ---
cat("Loading sequences from:", sequence_fasta_file, "\n")
tic_fasta_load <- Sys.time() # VERBOSE
seqs_raw <- tryCatch(readAAStringSet(sequence_fasta_file),
                     error = function(e) stop("Error reading FASTA file: ", e$message))
fasta_headers <- names(seqs_raw)
seq_lengths <- width(seqs_raw)
toc_fasta_load <- Sys.time() # VERBOSE
cat(sprintf("VERBOSE: FASTA loading took %.2f seconds.\n", as.numeric(difftime(toc_fasta_load, tic_fasta_load, units="secs")))) # VERBOSE
cat("Total sequences loaded from FASTA:", length(fasta_headers), "\n")
if(length(fasta_headers) == 0) stop("Error: No sequences found in the FASTA file.")

cat("VERBOSE: Starting sequence filtering...\n") # VERBOSE
cat("Filtering sequences to match distance matrix labels...\n")
# Get specimen IDs from distance matrix labels
dist_label_ids <- sub("-.+", "", dist_labels)
dist_label_to_id_map <- setNames(dist_label_ids, dist_labels)
# Find which specimen IDs from the distance matrix are present in the FASTA headers
ids_in_dist_matrix <- unique(dist_label_ids)
ids_to_keep <- intersect(ids_in_dist_matrix, fasta_headers)
cat("Found", length(ids_to_keep), "specimen IDs common to distance matrix and FASTA.\n")
if (length(ids_to_keep) == 0) {
    stop("Error: No common specimen IDs found between distance matrix labels and FASTA headers.")
}

# Filter the raw sequences
seqs_filtered <- seqs_raw[ids_to_keep]
filtered_fasta_headers <- names(seqs_filtered)
filtered_seq_lengths <- width(seqs_filtered)
cat("Filtered sequences down to", length(filtered_fasta_headers), "entries based on distance matrix IDs.\n")

# --- Create Sequence Lookup Maps (using FILTERED data) ---
cat("Creating sequence lookup maps from filtered sequences...\n")
seq_map_df <- data.frame(
    specimen_id = filtered_fasta_headers,
    sequence = as.character(seqs_filtered),
    seq_length = filtered_seq_lengths,
    stringsAsFactors = FALSE
) %>% filter(!is.na(sequence) & nchar(sequence) > 0) # Ensure no empty sequences included

# Create a named vector map: specimen_id -> sequence
seq_map <- setNames(seq_map_df$sequence, seq_map_df$specimen_id)
cat("Sequence lookup maps created using", nrow(seq_map_df), "non-empty filtered sequences.\n")
if (nrow(seq_map_df) == 0) {
    stop("Error: No valid sequences remained after filtering based on distance matrix labels.")
}

# --- Label Consistency Checks & Common Labels (using FILTERED sequence IDs) ---
cat("VERBOSE: Starting label consistency checks.\n") # VERBOSE
cat("Checking label consistency using filtered sequences...\n")
# Find labels common to dist matrix and tree
common_labels_tree_dist <- intersect(dist_labels, tree_tip_labels)
cat("Found", length(common_labels_tree_dist), "common labels between distance matrix and tree.\n")
if (length(common_labels_tree_dist) < 3) stop("Fewer than 3 common labels between distance matrix and tree.")

# Find dist matrix labels that have a corresponding sequence ID in the FILTERED seq_map
filtered_specimen_ids <- names(seq_map) # These are the IDs we actually have sequences for
labels_with_sequences <- names(dist_label_to_id_map[dist_label_to_id_map %in% filtered_specimen_ids])
num_labels_with_seqs <- length(labels_with_sequences)
cat("Found", num_labels_with_seqs, "distance matrix labels with matching sequences in the filtered set.\n")
if (num_labels_with_seqs == 0) stop("No overlapping IDs found between distance matrix labels and the filtered FASTA headers.")
if (num_labels_with_seqs < 3) stop("Fewer than 3 labels with matching filtered sequences.")

# Define the final set of common labels across distance matrix, tree, and FILTERED sequences
common_labels_all <- intersect(common_labels_tree_dist, labels_with_sequences)
cat("Found", length(common_labels_all), "common labels across distance matrix, tree, and filtered sequences.\n")
if (length(common_labels_all) < 3) stop("Fewer than 3 common labels across all inputs after filtering.")

# --- Subset Data to Common Labels ---
cat("VERBOSE: Starting data subsetting to common labels.\n") # VERBOSE
cat("Subsetting data to common labels...\n")
# Subset clustering distance matrix
clustering_dist_matrix_common <- clustering_dist_matrix[common_labels_all, common_labels_all]

# --- IMPORTANT MEMORY OPTIMIZATION ---
# Since the input dist_matrix_rds is confirmed to be the patristic distance matrix
# (or should be treated as such for silhouette evaluation), we will use it directly
# for silhouette calculations instead of recalculating from the tree. This avoids
# holding two very large distance matrices in memory.
# The tree_file is still used for label consistency checks and to define common_labels_all.
cat("Using the input distance matrix (subsetted to common labels) for Silhouette calculations to save memory.\n")
patristic_matrix_common <- clustering_dist_matrix_common # Use the input matrix for silhouette

# Subset tree (still useful for ensuring labels match a tree structure, even if not for cophenetic)
phylo_tree_common <- keep.tip(phylo_tree, common_labels_all)
# The following check is now on clustering_dist_matrix_common which is assigned to patristic_matrix_common
if (!is.matrix(patristic_matrix_common) || !identical(sort(rownames(patristic_matrix_common)), sort(common_labels_all))) {
    stop("The provided/subsetted distance matrix for silhouette calculation is invalid or its labels do not match common labels.")
}
# Subset sequence map data using the common labels (already effectively done by filtering earlier)
# Find the specimen IDs corresponding to the common labels
common_specimen_ids <- dist_label_to_id_map[common_labels_all]
# Filter the seq_map_df and seq_map (should already be filtered, but double-check)
seq_map_df_common <- seq_map_df %>% filter(specimen_id %in% common_specimen_ids)
seq_map_common <- seq_map[names(seq_map) %in% common_specimen_ids]
# Update label-to-id map for common labels only
tip_label_to_id_map_common <- dist_label_to_id_map[common_labels_all]

# Final check: ensure the number of sequences matches the number of common labels
if (length(seq_map_common) != length(common_labels_all)) {
   warning(paste("Mismatch between number of common labels (", length(common_labels_all),
                 ") and number of sequences kept (", length(seq_map_common),
                 "). This might indicate duplicate specimen IDs mapping to the same label set.", sep=""))
   # Depending on the logic, this might be acceptable or require further investigation/error handling
}

cat("Data subsetting complete. Using", length(common_labels_all), "tips for analysis.\n")


# --- Helper Functions ---
# calculate_cluster_identity function (identical to previous versions)
calculate_cluster_identity <- function(cluster_tip_labels, label_to_id_map, sequence_lookup_map, min_size = 3) {
    # ... (function body as in bayes_treeslice.r) ...
    current_call_context <- sys.call(-1); is_final_calc <- !identical(as.character(current_call_context[[1]]), "lapply") && !identical(as.character(current_call_context[[1]]), "parLapply"); debug_prefix <- if (is_final_calc) "  [calc_id_FINAL]" else "  [calc_id_BO]"
    if (is_final_calc) cat(sprintf("%s Processing cluster with %d labels.\n", debug_prefix, length(cluster_tip_labels)))
    specimen_ids <- label_to_id_map[cluster_tip_labels]; valid_specimen_ids <- specimen_ids[specimen_ids %in% names(sequence_lookup_map)]; valid_specimen_ids <- unique(valid_specimen_ids)
    if (is_final_calc) cat(sprintf("%s Found %d valid specimen IDs.\n", debug_prefix, length(valid_specimen_ids)))
    if (length(valid_specimen_ids) < min_size) { if (is_final_calc) cat(sprintf("%s Cluster size (%d) below min_size (%d). Returning NA.\n", debug_prefix, length(valid_specimen_ids), min_size)); return(NA_real_) }
    cluster_seqs <- tryCatch({ AAStringSet(sequence_lookup_map[valid_specimen_ids]) }, error = function(e){ if (is_final_calc) cat(sprintf("%s Error creating AAStringSet: %s\n", debug_prefix, e$message)); return(NULL) })
    if(is.null(cluster_seqs) || length(cluster_seqs) < min_size){ if (is_final_calc) cat(sprintf("%s Failed to create valid AAStringSet. Returning NA.\n", debug_prefix)); return(NA_real_) }
    if (is_final_calc) cat(sprintf("%s Attempting RemoveGaps for %d sequences...\n", debug_prefix, length(cluster_seqs)))
    cluster_seqs_nogap <- tryCatch({ RemoveGaps(cluster_seqs) }, error = function(e) { if (is_final_calc) cat(sprintf("%s RemoveGaps failed: %s\n", debug_prefix, e$message)); return(NULL) })
    if(is.null(cluster_seqs_nogap) || length(cluster_seqs_nogap) < min_size){ if (is_final_calc) cat(sprintf("%s RemoveGaps returned NULL or size < min_size. Returning NA.\n", debug_prefix)); return(NA_real_) }
    if (is_final_calc) cat(sprintf("%s RemoveGaps successful. Result length: %d\n", debug_prefix, length(cluster_seqs_nogap)))
    if (is_final_calc) cat(sprintf("%s Attempting AlignSeqs for %d sequences...\n", debug_prefix, length(cluster_seqs_nogap)))
    aligned_aas <- tryCatch({ AlignSeqs(cluster_seqs_nogap, verbose = FALSE, iterations = 0, refinement = 0) }, error = function(e) { if (is_final_calc) cat(sprintf("%s AlignSeqs failed: %s\n", debug_prefix, e$message)); return(NULL) })
    if (is.null(aligned_aas)) { if (is_final_calc) cat(sprintf("%s AlignSeqs returned NULL. Returning NA.\n", debug_prefix)); return(NA_real_) }
    if (is_final_calc) cat(sprintf("%s AlignSeqs successful. Result length: %d\n", debug_prefix, length(aligned_aas)))
    dist_mat <- tryCatch({ if (is_final_calc) cat(sprintf("%s Attempting DistanceMatrix...\n", debug_prefix)); res <- suppressWarnings(DistanceMatrix(aligned_aas, type = "matrix", correction="none", verbose = FALSE)); if (is_final_calc) cat(sprintf("%s DistanceMatrix successful. Matrix dim: %s\n", debug_prefix, paste(dim(res), collapse="x"))); res }, error = function(e) { if (is_final_calc) cat(sprintf("%s DistanceMatrix failed: %s\n", debug_prefix, e$message)); return(NULL) })
    if (is.null(dist_mat) || !is.matrix(dist_mat) || nrow(dist_mat) < 2) { if (is_final_calc) cat(sprintf("%s DistanceMatrix result invalid. Returning NA.\n", debug_prefix)); return(NA_real_) }
    pairwise_identities <- 1 - dist_mat[upper.tri(dist_mat)]; if (length(pairwise_identities) == 0) { if (is_final_calc) cat(sprintf("%s No pairwise identities calculated. Returning NA.\n", debug_prefix)); return(NA_real_) }
    mean_identity <- mean(pairwise_identities, na.rm = TRUE); if (is_final_calc) cat(sprintf("%s Calculated mean identity: %.4f\n", debug_prefix, mean_identity)); return(mean_identity)
}


# --- Bayesian Optimization Objective Function ---
# Note: BO optimizes continuous variables. We optimize k_cont and round it to integer k.
evaluate_k <- function(k_cont) {
  k <- max(2, round(k_cont)) # Ensure k is integer >= 2
  # VERBOSE: Changed cat to message for BO progress to distinguish from other cats
  message(sprintf("--- Evaluating k=%d (k_cont=%.2f) ---", k, k_cont))

  tryCatch({
    # Get needed variables from global environment
    # message(sprintf("k=%d: Getting variables from global env...", k)) # VERBOSE - can be too noisy
    clustering_dist_matrix_local <- get("clustering_dist_matrix_common", envir = .GlobalEnv) # Subsetted dist matrix
    patristic_matrix_local       <- get("patristic_matrix_common", envir = .GlobalEnv) # Patristic matrix based on common_labels_all
    common_labels_local          <- get("common_labels_all", envir = .GlobalEnv) # Common labels (dist, tree, seq)
    tip_label_to_id_map_local    <- get("tip_label_to_id_map_common", envir = .GlobalEnv)
    sequence_map_local           <- get("seq_map_common", envir = .GlobalEnv)
    target_median_min_local      <- get("target_median_min", envir = .GlobalEnv)
    target_median_max_local      <- get("target_median_max", envir = .GlobalEnv)
    target_mean_local            <- get("target_mean", envir = .GlobalEnv)
    min_cluster_size_for_msa_local <- get("min_cluster_size_for_msa", envir = .GlobalEnv)
    weight_median_local          <- get("weight_median", envir = .GlobalEnv)
    weight_mean_local            <- get("weight_mean", envir = .GlobalEnv)
    weight_identity_local        <- get("weight_identity", envir = .GlobalEnv)
    weight_silhouette_local      <- get("weight_silhouette", envir = .GlobalEnv)
    cat(sprintf("k=%d: Variables retrieved.\n", k)) # DEBUG GET VARS DONE

    # Check if there are enough data points after pruning for clustering
    num_data_points <- nrow(clustering_dist_matrix_local)
    # Define default failure return list structure
    failure_return <- list(Score = -Inf, avg_silhouette_width = NA_real_, n_clusters = k)

    if (num_data_points < 3 || k > num_data_points) {
        message(sprintf("k=%d: Insufficient data points (%d) or k (%d) is too large. Returning failure list.", k, num_data_points, k)) # VERBOSE
        return(failure_return)
    }

    # 1. Perform K-Medoids (PAM) clustering
    tic_pam_local <- Sys.time() # VERBOSE
    pam_result <- tryCatch({
        # message(sprintf("k=%d: Running pam...", k)) # VERBOSE - can be too noisy
        res <- cluster::pam(clustering_dist_matrix_local, k = k, diss = TRUE) # Use the clustering matrix
        # message(sprintf("k=%d: pam finished.", k)) # VERBOSE - can be too noisy
        res
    }, error = function(e) {
        warning(sprintf("PAM failed for k=%d: %s", k, e$message))
        return(NULL)
    })
    toc_pam_local <- Sys.time() # VERBOSE
    message(sprintf("  k=%d: pam() took %.2f sec.", k, as.numeric(difftime(toc_pam_local, tic_pam_local, units="secs")))) # VERBOSE

    if (is.null(pam_result)) {
        message(sprintf("k=%d: pam_result is NULL, returning failure list.", k)) # VERBOSE
        return(failure_return)
    }

    clusters_all <- pam_result$clustering # Vector: label -> cluster_id
    # message(sprintf("k=%d: Extracted clustering vector (length %d).", k, length(clusters_all))) # VERBOSE - can be too noisy
    clusters_all_df <- data.frame(tip_label = names(clusters_all), cluster_id = clusters_all, stringsAsFactors = FALSE)

    # Ensure labels match common_labels_local (should already match due to subsetting)
    if(!identical(sort(clusters_all_df$tip_label), sort(common_labels_local))) {
         warning(sprintf("Label mismatch after PAM for k=%d. Returning failure list.", k))
         return(failure_return)
    }

    # Order for silhouette calculation (using patristic matrix)
    clusters_common_df <- clusters_all_df[match(common_labels_local, clusters_all_df$tip_label), ]
    clusters_for_sil <- clusters_common_df$cluster_id
    names(clusters_for_sil) <- clusters_common_df$tip_label

    # Filter for Size/Identity calculations (only consider tips that had sequences and passed filtering)
    # Use common_labels_all as the base for filtering, as these are the tips present in all inputs after pruning
    cluster_df_seq_filtered <- clusters_all_df %>% filter(tip_label %in% common_labels_local)


    # 2. Calculate Cluster Sizes (based on taxa present after all filtering)
    cluster_sizes <- cluster_df_seq_filtered %>% group_by(cluster_id) %>% summarise(size = n(), .groups = 'drop')
    median_size <- if(nrow(cluster_sizes) > 0) median(cluster_sizes$size) else 0
    mean_size <- if(nrow(cluster_sizes) > 0) mean(cluster_sizes$size) else 0

    # 3. Calculate Sequence Identity (based on taxa present after all filtering)
    # Only calculate identity for clusters with members that passed the gap filter and are in common_labels_all
    cluster_members_seq_filtered <- cluster_df_seq_filtered %>%
        group_by(cluster_id) %>%
        summarise(members = list(tip_label), size = n(), .groups = 'drop') %>%
        filter(size >= min_cluster_size_for_msa_local) # Filter by min size for MSA

    avg_identity <- NA_real_ # Default to NA
    if(nrow(cluster_members_seq_filtered) > 0) {
        apply_identity_func <- function(ml) calculate_cluster_identity(ml, tip_label_to_id_map_local, sequence_map_local, min_cluster_size_for_msa_local)
        cluster_identities <- lapply(cluster_members_seq_filtered$members, apply_identity_func)
        valid_identities <- na.omit(unlist(cluster_identities))
        avg_identity <- if(length(valid_identities) > 0) mean(valid_identities) else 0 # Default to 0 if no valid identities calculated
    } else {
        avg_identity <- 0 # Default to 0 if no clusters meet min size for MSA
    }

    # 4. Calculate Silhouette Score (using PATRISTIC distances)
    avg_silhouette_width <- -1 # Default to -1 (results in score_sil = 0)
    num_unique_clusters_common <- length(unique(clusters_for_sil))
    # Check if enough data points and clusters exist for silhouette calculation
    if (length(clusters_for_sil) >= 3 && num_unique_clusters_common >= 2 && num_unique_clusters_common <= length(clusters_for_sil)) {
        sil_result <- tryCatch(cluster::silhouette(clusters_for_sil, dmatrix = patristic_matrix_local), error = function(e) {
            cat(sprintf("k=%d: Silhouette calculation failed: %s\n", k, e$message)); return(NULL)
        })
        if (!is.null(sil_result)) {
            sil_summary <- summary(sil_result)
            if (!is.null(sil_summary$avg.width) && is.numeric(sil_summary$avg.width)) {
                avg_silhouette_width <- sil_summary$avg.width
            } else {
                 cat(sprintf("k=%d: Silhouette summary avg.width is invalid.\n", k))
            }
        } else {
             cat(sprintf("k=%d: Silhouette result is NULL.\n", k))
        }
    } else {
        cat(sprintf("k=%d: Not enough data points (%d) or unique clusters (%d) for silhouette calc.\n", k, length(clusters_for_sil), num_unique_clusters_common))
    } # If silhouette calculation was not possible, avg_silhouette_width remains -1

    # 5. Calculate Component Scores - Ensure scores are numeric, default to 0 if NA/NaN. Use small negative penalty if calculation not possible.
    score_med <- if (is.na(median_size) || is.nan(median_size) || median_size == 0) -0.1 else { if (median_size >= target_median_min_local && median_size <= target_median_max_local) 1.0 else { dist_to_interval <- min(abs(median_size - target_median_min_local), abs(median_size - target_median_max_local)); penalty_scale <- (target_median_max_local - target_median_min_local) / 2 + 1; max(-0.1, 1 - (dist_to_interval / penalty_scale)^2 ) } } # Added NaN check and penalty
    score_mean <- if (is.na(mean_size) || is.nan(mean_size) || mean_size == 0) -0.1 else { sigma = (target_median_max_local - target_median_min_local) / 2; if(sigma <= 0) sigma <- 5.0; max(-0.1, exp(-( (mean_size - target_mean_local)^2 / (2 * sigma^2) )) ) } # Added NaN check and penalty
    score_ident <- if (is.na(avg_identity) || is.nan(avg_identity) || avg_identity == 0) -0.1 else min(avg_identity, 1.0) # Explicit check for NaN and penalty
    # score_sil is calculated based on avg_silhouette_width, which defaults to -1 if silhouette calc failed
    score_sil <- if (is.na(avg_silhouette_width) || is.nan(avg_silhouette_width) || avg_silhouette_width == -1) -0.1 else (avg_silhouette_width + 1) / 2 # Explicit check for NaN and penalty

    # 6. Combine Scores
    final_score <- (weight_median_local * score_med) + (weight_mean_local * score_mean) + (weight_identity_local * score_ident) + (weight_silhouette_local * score_sil)

    # Add a small term based on k to ensure score variation, especially in penalized scenarios
    # This helps the optimizer explore different k values even if primary metrics are poor
    final_score <- final_score + (k * 1e-6) # Add a small, k-dependent term

    # Debugging Output
    # VERBOSE: Changed cat to message for BO progress
    message(sprintf("  k=%d -> Med:%.2f(Sc:%.2f) Mean:%.2f(Sc:%.2f) Id:%.3f(Sc:%.2f) Sil:%.3f(Sc:%.2f) NClust:%d | FinalScore: %.4f", 
                  k, 
                  ifelse(is.na(median_size), -1, median_size), score_med, 
                  ifelse(is.na(mean_size), -1, mean_size), score_mean, 
                  ifelse(is.na(avg_identity), -1, avg_identity), score_ident, 
                  ifelse(is.na(avg_silhouette_width), -999, avg_silhouette_width), score_sil, 
                  num_unique_clusters_common, final_score))

    # Prepare return list
    return_list <- list(Score = final_score,
                        avg_silhouette_width = ifelse(is.na(avg_silhouette_width), -999, avg_silhouette_width),
                        n_clusters = k)

    # Ensure the final score is not NA or NaN before returning
    if (is.na(return_list$Score) || is.nan(return_list$Score)) { # Check the Score element in the list
        message(sprintf("k=%d: Calculated return_list$Score is NA/NaN. Forcing Score to -Inf.", k)) # VERBOSE
        return_list$Score <- -Inf
    }

    # message(sprintf("k=%d: Returning Score=%.4f", k, return_list$Score)) # VERBOSE - can be too noisy

    # Return results
    return(return_list)

  }, warning = function(w) {
      # Warnings should ideally not stop the flow or change the return structure
      # If a warning implies failure, it should ideally be caught as an error or lead to NA scores
      cat(sprintf("Warning evaluating k = %d (k_cont=%.2f): %s\n", round(k_cont), k_cont, w$message));
      invokeRestart("muffleWarning")
      # Let the function continue, potentially returning NA/Inf scores which are handled later
  }, error = function(e) {
      # Ensure the error handler also returns the consistent list structure
      cat(sprintf("Error evaluating k = %d (k_cont=%.2f): %s. Returning failure list.\n", round(k_cont), k_cont, e$message));
      # Use the same k as was passed into the function evaluation attempt
      k_err <- max(2, round(k_cont))
      return(list(Score = -Inf, avg_silhouette_width = NA_real_, n_clusters = k_err))
  })
}


# --- Run Bayesian Optimization ---
cat("\n--- Starting Bayesian Optimization using ParBayesianOptimization ---\n") # VERBOSE
# Define bounds for the continuous k variable
bounds <- list(k_cont = c(as.numeric(k_min), as.numeric(k_max)))
cat(sprintf("VERBOSE: Optimization bounds for k_cont: [%.2f, %.2f]\n", bounds$k_cont[1], bounds$k_cont[2])) # VERBOSE

# --- Set up parallel cluster ---
cl <- NULL
if (num_cores > 1) {
    cat("Setting up parallel cluster with", num_cores, "cores...\n")
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    on.exit({ if(!is.null(cl)) { stopCluster(cl); registerDoSEQ(); cl <- NULL }}, add = TRUE)
    cat("Exporting necessary objects/functions to cluster nodes...\n")
    # Export variables needed by evaluate_k and calculate_cluster_identity
    clusterExport(cl,
                  varlist=c("calculate_cluster_identity", "tip_label_to_id_map_common", "seq_map_common",
                            "min_cluster_size_for_msa", "clustering_dist_matrix_common",
                            "target_median_min", "target_median_max", "target_mean",
                            "weight_median", "weight_mean", "weight_identity",
                            "weight_silhouette", "patristic_matrix_common", "common_labels_all", # Pass the subsetted patristic matrix
                            "num_cores"),
                  envir=environment())
    cat("Loading required libraries on cluster nodes...\n")
    clusterEvalQ(cl, { library(Biostrings); library(DECIPHER); library(dplyr); library(cluster); library(ape); library(stringr) }) # Added stringr
    cat("Cluster setup complete.\n")
} else {
    cat("Running in serial mode (num_cores <= 1).\n")
}

# --- Run Optimization ---
cat("VERBOSE: Calling bayesOpt()...\n") # VERBOSE
tic_bayesopt <- Sys.time() # VERBOSE
bayesOpt_results <- tryCatch({
    set.seed(123)
    bayesOpt(
        FUN = evaluate_k, # Use the k-medoids evaluation function
        bounds = bounds,  # Use bounds for k_cont
        initPoints = bayesopt_init_points,
        iters.n = bayesopt_n_iter,
        iters.k = 8, # Default iters.k for ParBayesianOptimization
        parallel = (num_cores > 1),
        verbose = 2 # ParBayesianOptimization's own verbosity
    )
}, error = function(e) {
    cat("\nError during Bayesian Optimization:\n"); print(e); return(NULL)
})
toc_bayesopt <- Sys.time() # VERBOSE
cat(sprintf("VERBOSE: bayesOpt() call took %.2f seconds.\n", as.numeric(difftime(toc_bayesopt, tic_bayesopt, units="secs")))) # VERBOSE

# --- Stop cluster ---
if (!is.null(cl)) { cat("Stopping parallel cluster...\n"); stopCluster(cl); registerDoSEQ(); cl <- NULL }

cat("\n--- Bayesian Optimization Finished ---\n") # VERBOSE
if (is.null(bayesOpt_results)) stop("Bayesian Optimization failed.", call. = FALSE)

# --- Results ---
cat("VERBOSE: Processing optimization results...\n") # VERBOSE
best_params_df <- getBestPars(bayesOpt_results)
best_k_cont <- best_params_df$k_cont
best_k <- max(2, round(best_k_cont)) # Get the best integer k
best_score <- max(bayesOpt_results$scoreSummary$Score)
cat("Best k found (rounded):", best_k, "(continuous value:", best_k_cont, ")\n") # VERBOSE
cat("Best Score found:", best_score, "\n") # VERBOSE
opt_rds_file <- paste0(output_prefix, "_bayesopt_results.rds")
cat("Saving optimization results object to:", opt_rds_file, "\n") # VERBOSE
saveRDS(bayesOpt_results, file = opt_rds_file) # VERBOSE - actually save it

# --- Plotting ---
cat("VERBOSE: Preparing data for plots...\n") # VERBOSE
history_df <- bayesOpt_results$scoreSummary
# Add rounded k to history for plotting
history_df$k <- pmax(2, round(history_df$k_cont))

has_sil_col <- "avg_silhouette_width" %in% names(history_df)
has_nclust_col <- "n_clusters" %in% names(history_df) # Should match k

if(nrow(history_df) > 0) {
    cat("Generating optimization diagnostic plots...\n")
    p_progress <- ggplot(history_df, aes(x = Epoch, y = Score)) + geom_line(color = "blue") + geom_point(color = "blue") + labs(title = "A: Bayesian Optimization Progress", x = "Epoch", y = "Objective Score") + theme_bw()
    p_explore <- NULL
    if ("k_cont" %in% names(history_df)) { # Plot continuous k vs score
        p_explore <- ggplot(history_df, aes(x = k_cont, y = Score)) + geom_point(alpha=0.6, size=2) + geom_point(data=data.frame(k_cont=best_k_cont, Score=best_score), color='red', size=4, shape=17) + labs(title = "B: Explored k values vs Score", x = "k (continuous)", y = "Objective Score") + theme_bw()
    } else { cat("Warning: Could not find 'k_cont' column for explored points plot.\n") }
    p_sil_vs_k <- NULL # Changed plot name
    if (has_sil_col && has_nclust_col) {
        plot_data_sil <- history_df %>% filter(avg_silhouette_width > -998 & n_clusters > 1) %>% group_by(k) %>% summarise(avg_sil = mean(avg_silhouette_width), .groups='drop') # Avg sil per k
        if(nrow(plot_data_sil) > 0) {
             p_sil_vs_k <- ggplot(plot_data_sil, aes(x = k, y = avg_sil)) + geom_line(color="purple") + geom_point(color="purple", size=2) + labs(title = "C: Avg Silhouette Width vs k", x = "Number of Clusters (k)", y = "Avg Silhouette Width") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6))
        } else { cat("Warning: Not enough valid data points for Silhouette vs k plot.\n") }
    } else { cat("Warning: Missing columns for Silhouette vs k plot.\n") }
    # Combine plots
    combined_plot <- NULL; plot_elements <- list(p_progress, p_explore, p_sil_vs_k); plot_elements <- plot_elements[!sapply(plot_elements, is.null)]
    if (length(plot_elements) > 0) {
        tryCatch({
            if (length(plot_elements) == 3) { combined_plot <- plot_elements[[1]] / (plot_elements[[2]] + plot_elements[[3]]) }
            else if (length(plot_elements) == 2) { combined_plot <- plot_elements[[1]] + plot_elements[[2]] }
            else { combined_plot <- plot_elements[[1]] }
            combined_plot_file <- paste0(output_prefix, "_bayesopt_combined_plots.png")
            plot_width <- ifelse(length(plot_elements) > 1, 12, 7)
            plot_height <- ifelse(length(plot_elements) == 3, 8, 5)
            ggsave(combined_plot_file, plot = combined_plot, width = plot_width, height = plot_height)
            cat("Combined diagnostic plot saved to:", combined_plot_file, "\n")
        }, error = function(e) {
             cat("Error combining or saving plots:", e$message, "\n")
             # Fallback save omitted for brevity
        }) # End tryCatch for plot saving
    } else { # Corresponds to: if (length(plot_elements) > 0)
        cat("No diagnostic plots could be generated.\n")
    } # End if/else for plot elements existence
} else { # Corresponds to: if(nrow(history_df) > 0)
     cat("Optimization history is empty.\n")
} # End if/else for history_df rows


# --- Final Output using Best k ---
if (!is.null(best_k) && !is.na(best_k)) {
    cat("\n--- Generating final output for best k =", best_k, "---\n")

    # Run PAM with the best integer k
    final_pam_result <- cluster::pam(clustering_dist_matrix_common, k = best_k, diss = TRUE)
    final_clusters_vec <- final_pam_result$clustering
    final_clusters_df <- data.frame(
      tip_label = names(final_clusters_vec),
      cluster = final_clusters_vec,
      stringsAsFactors = FALSE
    )

    # --- Calculate Cluster Identities ---
    cat("Calculating identities for final clusters...\n")
    clusters_to_process <- final_clusters_df %>% group_by(cluster) %>% summarise(members = list(tip_label), size = n(), .groups = 'drop') %>% filter(size >= min_cluster_size_for_msa)
    if(nrow(clusters_to_process) > 0) {
        cat("Processing", nrow(clusters_to_process), "clusters for identity...\n"); pboptions(type = "timer", char = "=")
        identity_results <- pblapply(clusters_to_process$members, function(ml) { calculate_cluster_identity(ml, tip_label_to_id_map_common, seq_map_common, min_cluster_size_for_msa) })
        cluster_identity_data <- clusters_to_process %>% mutate(cluster_identity = unlist(identity_results)) %>% select(cluster, cluster_identity)
        cat("\nFinished calculating final cluster identities.\n")
        final_clusters_df <- final_clusters_df %>% left_join(cluster_identity_data, by = "cluster")
    } else { cat("No clusters met min size for identity calc.\n"); if (!"cluster_identity" %in% names(final_clusters_df)) final_clusters_df$cluster_identity <- NA_real_ }

    # --- Calculate Final Silhouette Score ---
    cat("Calculating Silhouette Score for the final clustering...\n")
    final_clusters_for_sil <- final_clusters_df$cluster; names(final_clusters_for_sil) <- final_clusters_df$tip_label
    num_unique_final_clusters <- length(unique(final_clusters_for_sil))
    if (length(final_clusters_for_sil) >= 3 && num_unique_final_clusters >= 2) {
        final_sil_result <- tryCatch({ cluster::silhouette(final_clusters_for_sil, dmatrix = patristic_matrix_common) }, error = function(e) { cat("Warning: Final silhouette calc failed:", e$message, "\n"); return(NULL) })
        if (!is.null(final_sil_result)) {
             # Add per-tip silhouette width
             sil_df <- data.frame(tip_label = names(final_clusters_for_sil), silhouette_width = final_sil_result[, "sil_width"], row.names = NULL, stringsAsFactors = FALSE)
             final_clusters_df <- final_clusters_df %>% left_join(sil_df, by = "tip_label")
             # Report average
             final_sil_summary <- summary(final_sil_result)
             if (!is.null(final_sil_summary$avg.width)) cat("Final Average Silhouette Width:", round(final_sil_summary$avg.width, 4), "\n")
        }
    } else { cat("Final clustering doesn't meet criteria for silhouette calc.\n") }
    if (!"silhouette_width" %in% names(final_clusters_df)) { final_clusters_df$silhouette_width <- NA_real_ }

    # --- Add Specimen ID and Sequence ---
    cat("Adding parsed ID and sequence to final output...\n")
    final_clusters_df <- final_clusters_df %>%
        mutate(specimen_id = tip_label_to_id_map_common[tip_label]) %>%
        left_join(seq_map_df_common, by = "specimen_id")

    # --- Separate Taxonomic Ranks ---
    cat("Attempting to add taxonomic breakdown...\n")
    tax_cols_hyphen <- c("Order", "Suborder", "Infraorder", "Family", "Subfamily", "Tribe", "Genus_Species_Combined")
    first_valid_label <- final_clusters_df$tip_label[!is.na(final_clusters_df$tip_label)][1]
    if (!is.na(first_valid_label) && grepl("-", sub("^[^-]+-?", "", first_valid_label))) {
         cat("Separating taxonomic ranks based on hyphens...\n")
         suppressWarnings( final_clusters_df <- final_clusters_df %>% filter(!is.na(tip_label)) %>% separate(tip_label, into = c("specimen_id_dup", tax_cols_hyphen), sep = "-", fill = "right", extra = "merge", remove = FALSE) %>% select(-specimen_id_dup) )
         added_cols_hyphen <- intersect(tax_cols_hyphen, names(final_clusters_df)); if(length(added_cols_hyphen) > 0){ final_clusters_df <- final_clusters_df %>% mutate(across(all_of(added_cols_hyphen), trimws)) }; cat("Added columns from hyphen separation:", paste(added_cols_hyphen, collapse=", "), "\n")
         if ("Genus_Species_Combined" %in% names(final_clusters_df)) {
             cat("Separating Genus and Species based on first underscore...\n")
             suppressWarnings( final_clusters_df <- final_clusters_df %>% separate(Genus_Species_Combined, into = c("Genus", "Species"), sep = "_", fill = "right", extra = "merge", remove = TRUE) )
             final_clusters_df <- final_clusters_df %>% mutate(across(any_of(c("Genus", "Species")), trimws)); cat("Separated Genus and Species columns.\n")
         } else { cat("Warning: 'Genus_Species_Combined' column not found.\n"); if (!"Genus" %in% names(final_clusters_df)) final_clusters_df$Genus <- NA_character_; if (!"Species" %in% names(final_clusters_df)) final_clusters_df$Species <- NA_character_ }
    } else { cat("Tip labels do not seem to have further '-' separators.\n") }
    final_tax_cols_expected <- c("Order", "Suborder", "Infraorder", "Family", "Subfamily", "Tribe", "Genus", "Species"); for(col in final_tax_cols_expected) { if (!col %in% names(final_clusters_df)) final_clusters_df[[col]] <- NA_character_ }

    # --- Reorder and Save Final Output ---
    cat("Reordering columns for final CSV...\n")
    final_tax_cols <- c("Order", "Suborder", "Infraorder", "Family", "Subfamily", "Tribe", "Genus", "Species")
    base_cols <- c("specimen_id", "tip_label", "cluster")
    taxonomic_cols_present <- intersect(final_tax_cols, names(final_clusters_df))
    stat_cols <- c("cluster_identity", "silhouette_width", "seq_length"); stat_cols_present <- intersect(stat_cols, names(final_clusters_df))
    sequence_col <- intersect("sequence", names(final_clusters_df))
    ordered_known_cols <- c(intersect(base_cols, names(final_clusters_df)), taxonomic_cols_present, stat_cols_present, sequence_col)
    remaining_cols <- setdiff(names(final_clusters_df), ordered_known_cols); final_col_order <- c(ordered_known_cols, remaining_cols)
    final_clusters_df <- final_clusters_df[, final_col_order]
    output_csv_name <- paste0(output_prefix, "_best_k", best_k, "_clusters.csv") # Use k in filename
    cat("Writing final cluster data to:", output_csv_name, "\n") # VERBOSE
    tryCatch({
        final_clusters_df <- final_clusters_df %>% mutate(across(any_of(c("cluster", "seq_length")), as.integer)) %>% mutate(across(any_of(c("cluster_identity", "silhouette_width")), as.numeric))
        write.csv(final_clusters_df, output_csv_name, row.names = FALSE, quote = TRUE, na = "")
    }, error = function(e) { stop("Error writing output CSV file: ", e$message) })

    # --- Generate Final Plots ---
    cat("VERBOSE: Generating final plots...\n") # VERBOSE
    # Histogram
    final_cluster_sizes <- final_clusters_df %>% filter(!is.na(sequence) & sequence != "") %>% count(cluster, name = "size")
    if(nrow(final_cluster_sizes) > 0){
        hist_plot_file <- paste0(output_prefix, "_best_k", best_k, "_histogram.png"); hist_plot_final <- ggplot(final_cluster_sizes, aes(x = size)) + geom_histogram(binwidth = 1, fill = "steelblue", color = "black") + labs(title = paste("Cluster Sizes (Taxa w/ Seq) at Optimal k =", best_k), subtitle = paste("Total clusters:", nrow(final_cluster_sizes), "| Total taxa:", sum(final_cluster_sizes$size)), x = "Cluster Size", y = "Frequency") + theme_bw() + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)); tryCatch(ggsave(hist_plot_file, plot = hist_plot_final, width = 8, height = 6), error=function(e){cat("Warning: Could not save final histogram plot:", e$message, "\n")}); cat("Final cluster size histogram saved to:", hist_plot_file, "\n") # VERBOSE
    } else { cat("No clusters with sequences for best k.\n") }
    # Cluster Silhouette Plot
    cat("Generating plot of average silhouette width per cluster...\n") # VERBOSE
    if ("silhouette_width" %in% names(final_clusters_df) && nrow(final_clusters_df %>% filter(!is.na(silhouette_width))) > 0) {
        cluster_summary_stats <- final_clusters_df %>% filter(!is.na(silhouette_width)) %>% group_by(cluster) %>% summarise(avg_cluster_silhouette = mean(silhouette_width, na.rm = TRUE), cluster_size = n(), .groups = 'drop') %>% filter(cluster_size > 0)
        if(nrow(cluster_summary_stats) > 0) {
            cluster_summary_stats$cluster <- factor(cluster_summary_stats$cluster); cluster_summary_stats <- cluster_summary_stats %>% arrange(desc(avg_cluster_silhouette)) %>% mutate(cluster = factor(cluster, levels = unique(cluster)))
            plot_cluster_sil_file <- paste0(output_prefix, "_best_k", best_k, "_cluster_silhouette.png")
            p_cluster_sil <- ggplot(cluster_summary_stats, aes(x = cluster, y = avg_cluster_silhouette, fill = cluster_size)) + geom_col(color = "black") + scale_fill_viridis_c(option = "plasma", name = "Cluster Size") + labs(title = paste("Average Silhouette Width per Cluster at Optimal k =", best_k), subtitle = paste("Total clusters plotted:", nrow(cluster_summary_stats)), x = "Cluster ID", y = "Average Silhouette Width") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6))
            tryCatch(ggsave(plot_cluster_sil_file, plot = p_cluster_sil, width = 10, height = 6), error=function(e){cat("Warning: Could not save cluster silhouette plot:", e$message, "\n")}); cat("Cluster silhouette plot saved to:", plot_cluster_sil_file, "\n")
        } else { cat("No valid cluster summary stats for silhouette plot.\n") }
    } else { cat("Silhouette widths not available for final data.\n") }

} else {
    cat("No best 'k' value found by optimization or optimization failed.\n")
}

    cat("Script finished successfully.\n") # VERBOSE
    script_end_time <- Sys.time() # VERBOSE
    cat(sprintf("bayes_kmedoids.r finished at: %s (Total runtime: %.2f minutes)\n", script_end_time, as.numeric(difftime(script_end_time, script_start_time, units="mins")))) # VERBOSE


# Print the best score in a parseable format (useful even without wrapper)
if (exists("best_score") && !is.null(best_score) && !is.na(best_score)) {
  cat(sprintf("FINAL_BEST_SCORE: %.10f\n", best_score))
} else {
  cat("FINAL_BEST_SCORE: NA\n") # Indicate failure or no score found
}
