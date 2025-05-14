#!/usr/bin/env Rscript

# --- Bayesian Optimization for Optimal Hierarchical Clustering Cutoff (h) ---
#
# Description:
# This script uses Bayesian Optimization to find an optimal cutoff height (h)
# for hierarchical clustering. It can take a distance object (dist) and perform
# hclust internally, or take a pre-computed hclust/dendrogram object.
# Optimization is based on multiple criteria:
#   - Median cluster size within a target range.
#   - Mean cluster size close to a target value.
#   - High average within-cluster sequence identity.
#   - High Silhouette Score based on PATRISTIC distances from an original tree.
#
# It requires:
#   1. An input RDS file: EITHER a distance object (class 'dist' or a symmetric matrix)
#      OR a pre-computed hierarchical clustering object (class 'hclust' or 'dendrogram').
#      Labels should follow 'SpecimenID-Rank1-Rank2...'.
#   2. The original phylogenetic tree file (e.g., Newick) for calculating patristic
#      distances used in the Silhouette score evaluation. Tree tip labels must *exactly* match
#      the labels in the input object (or derived hclust object).
#   3. A corresponding FASTA file (AA sequences) where headers match the *first part*
#      (SpecimenID, before the first hyphen) of the input object labels.
#      The script will internally filter this FASTA.
#
# Usage:
# Rscript bayes2.r <input_rds_file> <tree_file> <fasta_file> <output_prefix> \
#                          <h_min> <h_max> <init_pts> <n_iter> \
#                          <med_min> <med_max> <mean_target> <min_msa_size> \
#                          <w_median> <w_mean> <w_identity> <w_silhouette> \
#                          [num_cores] [hclust_method]
#
# Arguments:
#   input_rds_file     : Path to the input RDS file (dist, matrix, hclust, or dendrogram).
#   tree_file          : Path to the phylogenetic tree file (e.g., .tre, .newick).
#   fasta_file         : Path to the sequence file (.fasta, .fa) - AA sequences.
#   output_prefix      : Prefix for output files (e.g., "hcut_run").
#   h_min              : Minimum cutoff height (h) to search (numeric >= 0).
#   h_max              : Maximum cutoff height (h) to search (numeric).
#   init_pts           : Number of initial random points (h values) for BO (integer).
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
#   hclust_method      : [Optional] Agglomeration method for hclust if input is dist/matrix
#                        (e.g., "average", "complete", "ward.D2"). Default: "average".
#
# Example (with dist object input):
# Rscript bayes2.r my_dist_obj.rds tree.tre sequences.fasta hcut_run \
#                          0.1 10 20 50 5 12 10 3 0.3 0.15 0.45 0.1 8 average
#

# --- Dependencies ---
cat("Loading required packages...\n")
suppressPackageStartupMessages({
    library(ParBayesianOptimization)
    library(doParallel)
    library(DECIPHER)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(ggplot2)
    library(ape)
    library(Biostrings)
    library(parallel)
    library(cluster)
    library(pbapply)
    library(patchwork)
    library(viridis)
    library(dendextend)
})

# --- Argument Parsing ---
args <- commandArgs(trailingOnly = TRUE)

min_args <- 16 # Mandatory args up to w_silhouette
max_args <- 18 # Up to num_cores and hclust_method

if (length(args) < min_args || length(args) > max_args) {
  cat("Usage: Rscript bayes2.r <input_rds_file> <tree_file> <fasta_file> <output_prefix> \\\n")
  cat("                                <h_min> <h_max> <init_pts> <n_iter> \\\n")
  cat("                                <med_min> <med_max> <mean_target> <min_msa_size> \\\n")
  cat("                                <w_median> <w_mean> <w_identity> <w_silhouette> \\\n")
  cat("                                [num_cores] [hclust_method]\n\n")
  cat("Error: Incorrect number of arguments provided.\n")
  cat("Received", length(args), "arguments.\n")
  quit(status = 1)
}

# Assign arguments to variables
input_rds_file        <- args[1]
tree_file             <- args[2]
sequence_fasta_file   <- args[3]
output_prefix         <- args[4]
h_min                 <- as.numeric(args[5])
h_max                 <- as.numeric(args[6])
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

num_cores_arg <- NA
if (length(args) >= 17) {
  num_cores_arg <- as.integer(args[17])
}

hclust_method_arg <- "average" # Default hclust method
if (length(args) == 18) {
  hclust_method_arg <- args[18]
}

if (h_min < 0) stop("Error: h_min must be >= 0.")
if (h_max < h_min) stop("Error: h_max must be >= h_min.")

numeric_args_check <- list(h_min=h_min, h_max=h_max, bayesopt_init_points=bayesopt_init_points,
                     bayesopt_n_iter=bayesopt_n_iter, target_median_min=target_median_min,
                     target_median_max=target_median_max, target_mean=target_mean,
                     min_cluster_size_for_msa=min_cluster_size_for_msa, weight_median=weight_median,
                     weight_mean=weight_mean, weight_identity=weight_identity,
                     weight_silhouette=weight_silhouette, num_cores_arg=num_cores_arg)

for (name in names(numeric_args_check)) {
  val <- numeric_args_check[[name]]
  if (name == "num_cores_arg" && is.na(val)) next # num_cores can be NA if not provided
  if (is.null(val) || is.na(val) || (!is.numeric(val) && !is.integer(val)) ) {
    # Approximate original argument position for error message
    original_arg_pos <- match(name, names(numeric_args_check)) + 4 # +4 because first 4 args are files/prefix
    if (name == "num_cores_arg") original_arg_pos <- 17
    stop("Error: Argument '", name, "' must be a valid number/integer. Received: '", args[original_arg_pos], "'")
  }
}

if (!file.exists(input_rds_file)) stop("Error: Input RDS file not found: ", input_rds_file)
if (!file.exists(tree_file)) stop("Error: Tree file not found: ", tree_file)
if (!file.exists(sequence_fasta_file)) stop("Error: Sequence FASTA file not found: ", sequence_fasta_file)

if (!is.na(num_cores_arg) && num_cores_arg > 0) {
    num_cores <- num_cores_arg
} else {
    num_cores <- detectCores() - 1
    if (is.na(num_cores) || num_cores < 1) num_cores <- 1
}

cat("--- Configuration ---\n")
cat("Input RDS File       :", input_rds_file, "\n")
cat("Tree File            :", tree_file, "\n")
cat("Sequence FASTA File  :", sequence_fasta_file, "\n")
cat("Output Prefix        :", output_prefix, "\n")
cat("h Search Range       : [", h_min, ",", h_max, "]\n")
cat("BO Init Points       :", bayesopt_init_points, "\n")
cat("BO Iterations        :", bayesopt_n_iter, "\n")
cat("Target Median Range  : [", target_median_min, ",", target_median_max, "]\n")
cat("Target Mean Size     :", target_mean, "\n")
cat("Min Cluster for MSA  :", min_cluster_size_for_msa, "\n")
cat("Weight Median        :", weight_median, "\n")
cat("Weight Mean          :", weight_mean, "\n")
cat("Weight Identity      :", weight_identity, "\n")
cat("Weight Silhouette    :", weight_silhouette, "\n")
cat("Using", num_cores, "cores for parallel processing.\n")
cat("Hclust Method (if applicable) :", hclust_method_arg, "\n")
cat("---------------------\n")

# --- Load Data ---
cat("Loading input object from RDS:", input_rds_file, "\n")
loaded_object <- tryCatch({
    readRDS(input_rds_file)
}, error = function(e) {
    stop("Error reading input RDS file '", input_rds_file, "': ", e$message)
})

# Determine if hclust needs to be performed
if (inherits(loaded_object, "dist") || (is.matrix(loaded_object) && ncol(loaded_object) == nrow(loaded_object) && isSymmetric(loaded_object))) {
    cat("Input is a distance object/matrix. Performing hierarchical clustering using method:", hclust_method_arg, "...\n")
    dist_for_hclust <- loaded_object
    if (is.matrix(dist_for_hclust) && !inherits(dist_for_hclust, "dist")) {
        cat("Converting input matrix to dist object...\n")
        dist_for_hclust <- as.dist(dist_for_hclust)
    }
    hc_object_full <- tryCatch({
        hclust(dist_for_hclust, method = hclust_method_arg)
    }, error = function(e) {
        stop("Error performing hclust on the input distance object: ", e$message)
    })
    cat("Hierarchical clustering complete.\n")
} else if (inherits(loaded_object, c("hclust", "dendrogram"))) {
    cat("Input is already a hierarchical clustering object (hclust/dendrogram).\n")
    hc_object_full <- loaded_object
} else {
    stop("Input RDS file must contain an object of class 'dist', 'matrix' (distance matrix), 'hclust', or 'dendrogram'. Found: ", class(loaded_object)[1])
}

# Now hc_object_full is guaranteed to be an hclust or dendrogram object
hc_object_labels_full <- if (inherits(hc_object_full, "hclust")) hc_object_full$labels else labels(hc_object_full)
if (is.null(hc_object_labels_full) || length(hc_object_labels_full) == 0) stop("HC object (derived or loaded) must have labels.")
cat("Total labels in HC object:", length(hc_object_labels_full), "\n")


cat("Loading phylogenetic tree from:", tree_file, "\n")
phylo_tree_full <- tryCatch(ape::read.tree(tree_file), error = function(e) stop("Error loading tree: ", e$message))
tree_tip_labels_full <- phylo_tree_full$tip.label
cat("Total tips in tree:", length(tree_tip_labels_full), "\n")
if(length(tree_tip_labels_full) == 0) stop("Error: No tip labels found in the tree object.")

cat("Loading sequences from:", sequence_fasta_file, "\n")
seqs_raw <- tryCatch(readAAStringSet(sequence_fasta_file),
                     error = function(e) stop("Error reading FASTA file: ", e$message))
fasta_headers_raw <- names(seqs_raw)
cat("Total sequences loaded from FASTA:", length(fasta_headers_raw), "\n")
if(length(fasta_headers_raw) == 0) stop("Error: No sequences found in the FASTA file.")

cat("Filtering sequences to match HC object labels...\n")
hc_label_ids_parsed <- sub("-.+", "", hc_object_labels_full)
hc_label_to_id_map_full <- setNames(hc_label_ids_parsed, hc_object_labels_full)

ids_in_hc_object <- unique(hc_label_ids_parsed)
ids_to_keep_fasta <- intersect(ids_in_hc_object, fasta_headers_raw)
cat("Found", length(ids_to_keep_fasta), "specimen IDs common to HC object and FASTA.\n")
if (length(ids_to_keep_fasta) == 0) {
    stop("Error: No common specimen IDs found between HC object labels and FASTA headers.")
}

seqs_filtered <- seqs_raw[ids_to_keep_fasta]
filtered_fasta_headers <- names(seqs_filtered)
filtered_seq_lengths <- width(seqs_filtered)
cat("Filtered sequences down to", length(filtered_fasta_headers), "entries based on HC object IDs.\n")

cat("Creating sequence lookup maps from filtered sequences...\n")
seq_map_df_filtered <- data.frame(
    specimen_id = filtered_fasta_headers,
    sequence = as.character(seqs_filtered),
    seq_length = filtered_seq_lengths,
    stringsAsFactors = FALSE
) %>% filter(!is.na(sequence) & nchar(sequence) > 0)

seq_map_filtered_lookup <- setNames(seq_map_df_filtered$sequence, seq_map_df_filtered$specimen_id)
cat("Sequence lookup maps created using", nrow(seq_map_df_filtered), "non-empty filtered sequences.\n")
if (nrow(seq_map_df_filtered) == 0) {
    stop("Error: No valid sequences remained after filtering based on HC object labels.")
}

cat("Checking label consistency using filtered sequences...\n")
common_labels_tree_hc <- intersect(hc_object_labels_full, tree_tip_labels_full)
cat("Found", length(common_labels_tree_hc), "common labels between HC object and tree.\n")
if (length(common_labels_tree_hc) < 3) stop("Fewer than 3 common labels between HC object and tree.")

labels_with_sequences <- names(hc_label_to_id_map_full[hc_label_to_id_map_full %in% names(seq_map_filtered_lookup)])
num_labels_with_seqs <- length(labels_with_sequences)
cat("Found", num_labels_with_seqs, "HC object labels with matching sequences in the filtered set.\n")
if (num_labels_with_seqs == 0) stop("No overlapping IDs found between HC object labels and the filtered FASTA headers.")
if (num_labels_with_seqs < 3) stop("Fewer than 3 HC labels with matching filtered sequences.")

common_labels_all <- intersect(common_labels_tree_hc, labels_with_sequences)
cat("Found", length(common_labels_all), "common labels across HC object, tree, and filtered sequences.\n")
if (length(common_labels_all) < 3) stop("Fewer than 3 common labels across all inputs after filtering.")

cat("Subsetting data to common labels...\n")
phylo_tree_common <- keep.tip(phylo_tree_full, common_labels_all)
patristic_matrix_common <- tryCatch(cophenetic.phylo(phylo_tree_common), error = function(e) stop("Error calculating patristic distances on subsetted tree: ", e$message))
if (!is.matrix(patristic_matrix_common) || !identical(sort(rownames(patristic_matrix_common)), sort(common_labels_all))) {
    stop("Failed to calculate valid patristic matrix for common labels.")
}
common_specimen_ids <- hc_label_to_id_map_full[common_labels_all]
seq_map_df_common <- seq_map_df_filtered %>% filter(specimen_id %in% common_specimen_ids)
seq_map_common_lookup <- seq_map_filtered_lookup[names(seq_map_filtered_lookup) %in% common_specimen_ids]
tip_label_to_id_map_common <- hc_label_to_id_map_full[common_labels_all]

if (length(seq_map_common_lookup) != length(unique(common_specimen_ids))) {
   warning(paste("Mismatch between number of common specimen IDs (", length(unique(common_specimen_ids)),
                 ") and number of sequences kept (", length(seq_map_common_lookup),
                 "). This might indicate duplicate specimen IDs mapping to the same label set or issues in filtering."))
}
cat("Data subsetting complete. Using", length(common_labels_all), "tips for analysis where applicable.\n")

# --- Helper Functions ---
calculate_cluster_identity <- function(cluster_tip_labels, label_to_id_map, sequence_lookup_map, min_size = 3) {
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
evaluate_h <- function(h_cont) {
  h <- max(0, h_cont)
  tryCatch({
    hc_obj_for_cutree            <- get("hc_object_full", envir = .GlobalEnv) # This is now guaranteed hclust/dendrogram
    patristic_matrix_local       <- get("patristic_matrix_common", envir = .GlobalEnv)
    common_labels_local          <- get("common_labels_all", envir = .GlobalEnv)
    tip_label_to_id_map_local    <- get("tip_label_to_id_map_common", envir = .GlobalEnv)
    sequence_map_local           <- get("seq_map_common_lookup", envir = .GlobalEnv)
    target_median_min_local      <- get("target_median_min", envir = .GlobalEnv)
    target_median_max_local      <- get("target_median_max", envir = .GlobalEnv)
    target_mean_local            <- get("target_mean", envir = .GlobalEnv)
    min_cluster_size_for_msa_local <- get("min_cluster_size_for_msa", envir = .GlobalEnv)
    weight_median_local          <- get("weight_median", envir = .GlobalEnv)
    weight_mean_local            <- get("weight_mean", envir = .GlobalEnv)
    weight_identity_local        <- get("weight_identity", envir = .GlobalEnv)
    weight_silhouette_local      <- get("weight_silhouette", envir = .GlobalEnv)

    failure_return <- list(Score = -Inf, avg_silhouette_width = NA_real_, n_clusters = NA_integer_)

    clusters_raw <- tryCatch({ cutree(hc_obj_for_cutree, h = h) }, error = function(e) { NULL })
    if (is.null(clusters_raw) || length(clusters_raw) == 0) return(failure_return)
    
    clusters_all_df <- data.frame(tip_label = names(clusters_raw), cluster_id = clusters_raw, stringsAsFactors = FALSE)
    clusters_for_sil_df <- clusters_all_df %>% filter(tip_label %in% common_labels_local)
    clusters_for_sil_df <- clusters_for_sil_df[match(common_labels_local, clusters_for_sil_df$tip_label), ]
    clusters_for_sil_vec <- clusters_for_sil_df$cluster_id
    names(clusters_for_sil_vec) <- clusters_for_sil_df$tip_label
    cluster_df_metrics <- clusters_all_df %>% filter(tip_label %in% common_labels_local)

    cluster_sizes_df <- cluster_df_metrics %>% group_by(cluster_id) %>% summarise(size = n(), .groups = 'drop')
    median_size <- if(nrow(cluster_sizes_df) > 0) median(cluster_sizes_df$size) else 0
    mean_size <- if(nrow(cluster_sizes_df) > 0) mean(cluster_sizes_df$size) else 0

    cluster_members_for_identity <- cluster_df_metrics %>%
        group_by(cluster_id) %>%
        summarise(members = list(tip_label), size = n(), .groups = 'drop') %>%
        filter(size >= min_cluster_size_for_msa_local)
    avg_identity <- NA_real_
    if(nrow(cluster_members_for_identity) > 0) {
        apply_identity_func <- function(ml) calculate_cluster_identity(ml, tip_label_to_id_map_local, sequence_map_local, min_cluster_size_for_msa_local)
        cluster_identities <- lapply(cluster_members_for_identity$members, apply_identity_func)
        valid_identities <- na.omit(unlist(cluster_identities))
        avg_identity <- if(length(valid_identities) > 0) mean(valid_identities) else 0
    } else { avg_identity <- 0 }

    avg_silhouette_width <- -1
    num_unique_clusters_sil <- length(unique(clusters_for_sil_vec))
    if (length(clusters_for_sil_vec) >= 3 && num_unique_clusters_sil >= 2 && num_unique_clusters_sil <= length(clusters_for_sil_vec)) {
        sil_result <- tryCatch(cluster::silhouette(clusters_for_sil_vec, dmatrix = patristic_matrix_local), error = function(e) { NULL })
        if (!is.null(sil_result)) {
            sil_summary <- summary(sil_result)
            if (!is.null(sil_summary$avg.width) && is.numeric(sil_summary$avg.width)) {
                avg_silhouette_width <- sil_summary$avg.width
            }
        }
    }

    score_med <- if (is.na(median_size) || is.nan(median_size) || median_size == 0) -0.1 else { if (median_size >= target_median_min_local && median_size <= target_median_max_local) 1.0 else { dist_to_interval <- min(abs(median_size - target_median_min_local), abs(median_size - target_median_max_local)); penalty_scale <- (target_median_max_local - target_median_min_local) / 2 + 1; max(-0.1, 1 - (dist_to_interval / penalty_scale)^2 ) } }
    score_mean <- if (is.na(mean_size) || is.nan(mean_size) || mean_size == 0) -0.1 else { sigma = (target_median_max_local - target_median_min_local) / 2; if(sigma <= 0) sigma <- 5.0; max(-0.1, exp(-( (mean_size - target_mean_local)^2 / (2 * sigma^2) )) ) }
    score_ident <- if (is.na(avg_identity) || is.nan(avg_identity) || avg_identity == 0) -0.1 else min(avg_identity, 1.0)
    score_sil <- if (is.na(avg_silhouette_width) || is.nan(avg_silhouette_width) || avg_silhouette_width == -1) -0.1 else (avg_silhouette_width + 1) / 2
    final_score <- (weight_median_local * score_med) + (weight_mean_local * score_mean) + (weight_identity_local * score_ident) + (weight_silhouette_local * score_sil)
    final_score <- final_score + (h * 1e-7)

    return_list <- list(Score = final_score,
                        avg_silhouette_width = ifelse(is.na(avg_silhouette_width) || avg_silhouette_width == -1, -999, avg_silhouette_width),
                        n_clusters = num_unique_clusters_sil)
    if (is.na(return_list$Score) || is.nan(return_list$Score)) { return_list$Score <- -Inf }
    return(return_list)
  }, warning = function(w) {
      cat(sprintf("Warning evaluating h = %.4f (h_cont=%.4f): %s\n", max(0,h_cont), h_cont, w$message));
      invokeRestart("muffleWarning")
  }, error = function(e) {
      cat(sprintf("Error evaluating h = %.4f (h_cont=%.4f): %s. Returning failure list.\n", max(0,h_cont), h_cont, e$message));
      return(list(Score = -Inf, avg_silhouette_width = NA_real_, n_clusters = NA_integer_))
  })
}

# --- Run Bayesian Optimization ---
cat("\n--- Starting Bayesian Optimization using ParBayesianOptimization ---\n")
bounds <- list(h_cont = c(as.numeric(h_min), as.numeric(h_max)))
cl <- NULL
if (num_cores > 1) {
    cat("Setting up parallel cluster with", num_cores, "cores...\n")
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    on.exit({ if(!is.null(cl)) { stopCluster(cl); registerDoSEQ(); cl <- NULL }}, add = TRUE)
    cat("Exporting necessary objects/functions to cluster nodes...\n")
    clusterExport(cl,
                  varlist=c("calculate_cluster_identity", "tip_label_to_id_map_common", "seq_map_common_lookup",
                            "min_cluster_size_for_msa", "hc_object_full",
                            "target_median_min", "target_median_max", "target_mean",
                            "weight_median", "weight_mean", "weight_identity",
                            "weight_silhouette", "patristic_matrix_common", "common_labels_all",
                            "num_cores"),
                  envir=environment())
    cat("Loading required libraries on cluster nodes...\n")
    clusterEvalQ(cl, { library(Biostrings); library(DECIPHER); library(dplyr); library(cluster); library(ape); library(stringr); library(dendextend) })
    cat("Cluster setup complete.\n")
} else { cat("Running in serial mode (num_cores <= 1).\n") }

bayesOpt_results <- tryCatch({
    set.seed(123)
    bayesOpt( FUN = evaluate_h, bounds = bounds, initPoints = bayesopt_init_points, iters.n = bayesopt_n_iter, iters.k = 8, parallel = (num_cores > 1), verbose = 2 )
}, error = function(e) { cat("\nError during Bayesian Optimization:\n"); print(e); return(NULL) })
if (!is.null(cl)) { cat("Stopping parallel cluster...\n"); stopCluster(cl); registerDoSEQ(); cl <- NULL }
cat("\n--- Bayesian Optimization Finished ---\n")
if (is.null(bayesOpt_results)) stop("Bayesian Optimization failed.", call. = FALSE)

best_params_df <- getBestPars(bayesOpt_results)
best_h_cont <- best_params_df$h_cont
best_h <- max(0, best_h_cont)
best_score <- max(bayesOpt_results$scoreSummary$Score)
cat("Best h found (raw):", best_h, "(continuous value:", best_h_cont, ")\n")
cat("Best Score found:", best_score, "\n")
opt_rds_file <- paste0(output_prefix, "_bayesopt_results.rds")
saveRDS(bayesOpt_results, file = opt_rds_file)
cat("Optimization results object saved to:", opt_rds_file, "\n")

history_df <- bayesOpt_results$scoreSummary
history_df$h <- pmax(0, history_df$h_cont)
has_sil_col <- "avg_silhouette_width" %in% names(history_df)
if(nrow(history_df) > 0) {
    cat("Generating optimization diagnostic plots...\n")
    p_progress <- ggplot(history_df, aes(x = Epoch, y = Score)) + geom_line(color = "blue") + geom_point(color = "blue") + labs(title = "A: Bayesian Optimization Progress", x = "Epoch", y = "Objective Score") + theme_bw()
    p_explore <- NULL
    if ("h_cont" %in% names(history_df)) {
        p_explore <- ggplot(history_df, aes(x = h_cont, y = Score)) + geom_point(alpha=0.6, size=2) + geom_point(data=data.frame(h_cont=best_h_cont, Score=best_score), color='red', size=4, shape=17) + labs(title = "B: Explored h values vs Score", x = "h (continuous)", y = "Objective Score") + theme_bw()
    } else { cat("Warning: Could not find 'h_cont' column for explored points plot.\n") }
    p_sil_vs_h <- NULL
    if (has_sil_col && "h" %in% names(history_df)) {
        plot_data_sil_h <- history_df %>% filter(avg_silhouette_width > -998) %>% group_by(h) %>% summarise(avg_sil = mean(avg_silhouette_width, na.rm=TRUE), .groups='drop')
        if(nrow(plot_data_sil_h) > 0) {
             p_sil_vs_h <- ggplot(plot_data_sil_h, aes(x = h, y = avg_sil)) + geom_line(color="purple") + geom_point(color="purple", size=2) + labs(title = "C: Avg Silhouette Width vs h", x = "Cutoff Height (h)", y = "Avg Silhouette Width") + theme_bw()
        } else { cat("Warning: Not enough valid data points for Silhouette vs h plot.\n") }
    } else { cat("Warning: Missing columns for Silhouette vs h plot.\n") }
    combined_plot <- NULL; plot_elements <- list(p_progress, p_explore, p_sil_vs_h); plot_elements <- plot_elements[!sapply(plot_elements, is.null)]
    if (length(plot_elements) > 0) {
        tryCatch({
            if (length(plot_elements) == 3) { combined_plot <- plot_elements[[1]] / (plot_elements[[2]] + plot_elements[[3]]) }
            else if (length(plot_elements) == 2) { combined_plot <- plot_elements[[1]] + plot_elements[[2]] }
            else { combined_plot <- plot_elements[[1]] }
            combined_plot_file <- paste0(output_prefix, "_bayesopt_combined_plots.png")
            plot_width <- ifelse(length(plot_elements) > 1, 12, 7); plot_height <- ifelse(length(plot_elements) == 3, 8, 5)
            ggsave(combined_plot_file, plot = combined_plot, width = plot_width, height = plot_height)
            cat("Combined diagnostic plot saved to:", combined_plot_file, "\n")
        }, error = function(e) { cat("Error combining or saving plots:", e$message, "\n") })
    } else { cat("No diagnostic plots could be generated.\n") }
} else { cat("Optimization history is empty.\n") }

# --- Final Output using Best h ---
if (!is.null(best_h) && !is.na(best_h)) {
    cat("\n--- Generating final output for best h =", round(best_h, 4), "---\n")
    final_clusters_raw <- cutree(hc_object_full, h = best_h)
    final_clusters_df <- data.frame( tip_label = names(final_clusters_raw), cluster = final_clusters_raw, stringsAsFactors = FALSE )
    cat("Calculating identities for final clusters...\n")
    clusters_to_process_final <- final_clusters_df %>% filter(tip_label %in% common_labels_all) %>% group_by(cluster) %>% summarise(members = list(tip_label), size = n(), .groups = 'drop') %>% filter(size >= min_cluster_size_for_msa)
    if(nrow(clusters_to_process_final) > 0) {
        cat("Processing", nrow(clusters_to_process_final), "clusters for identity...\n"); pboptions(type = "timer", char = "=")
        identity_results_final <- pblapply(clusters_to_process_final$members, function(ml) { calculate_cluster_identity(ml, tip_label_to_id_map_common, seq_map_common_lookup, min_cluster_size_for_msa) })
        cluster_identity_data_final <- clusters_to_process_final %>% mutate(cluster_identity = unlist(identity_results_final)) %>% select(cluster, cluster_identity)
        cat("\nFinished calculating final cluster identities.\n")
        final_clusters_df <- final_clusters_df %>% left_join(cluster_identity_data_final, by = "cluster")
    } else { cat("No clusters met min size for identity calc.\n"); if (!"cluster_identity" %in% names(final_clusters_df)) final_clusters_df$cluster_identity <- NA_real_ }
    cat("Calculating Silhouette Score for the final clustering...\n")
    final_clusters_for_sil_df <- final_clusters_df %>% filter(tip_label %in% common_labels_all)
    final_clusters_for_sil_df <- final_clusters_for_sil_df[match(common_labels_all, final_clusters_for_sil_df$tip_label), ]
    final_clusters_for_sil_vec <- final_clusters_for_sil_df$cluster; names(final_clusters_for_sil_vec) <- final_clusters_for_sil_df$tip_label
    num_unique_final_clusters <- length(unique(final_clusters_for_sil_vec))
    if (length(final_clusters_for_sil_vec) >= 3 && num_unique_final_clusters >= 2) {
        final_sil_result <- tryCatch({ cluster::silhouette(final_clusters_for_sil_vec, dmatrix = patristic_matrix_common) }, error = function(e) { cat("Warning: Final silhouette calc failed:", e$message, "\n"); return(NULL) })
        if (!is.null(final_sil_result)) {
             sil_df_final <- data.frame(tip_label = names(final_clusters_for_sil_vec), silhouette_width = final_sil_result[, "sil_width"], row.names = NULL, stringsAsFactors = FALSE)
             final_clusters_df <- final_clusters_df %>% left_join(sil_df_final, by = "tip_label")
             final_sil_summary <- summary(final_sil_result)
             if (!is.null(final_sil_summary$avg.width)) cat("Final Average Silhouette Width:", round(final_sil_summary$avg.width, 4), "\n")
        }
    } else { cat("Final clustering doesn't meet criteria for silhouette calc.\n") }
    if (!"silhouette_width" %in% names(final_clusters_df)) { final_clusters_df$silhouette_width <- NA_real_ }
    cat("Adding parsed ID and sequence to final output...\n")
    final_clusters_df <- final_clusters_df %>% mutate(specimen_id = hc_label_to_id_map_full[tip_label]) %>% left_join(seq_map_df_filtered, by = "specimen_id")
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
    cat("Reordering columns for final CSV...\n")
    final_tax_cols <- c("Order", "Suborder", "Infraorder", "Family", "Subfamily", "Tribe", "Genus", "Species")
    base_cols <- c("specimen_id", "tip_label", "cluster")
    taxonomic_cols_present <- intersect(final_tax_cols, names(final_clusters_df))
    stat_cols <- c("cluster_identity", "silhouette_width", "seq_length"); stat_cols_present <- intersect(stat_cols, names(final_clusters_df))
    sequence_col <- intersect("sequence", names(final_clusters_df))
    ordered_known_cols <- c(intersect(base_cols, names(final_clusters_df)), taxonomic_cols_present, stat_cols_present, sequence_col)
    remaining_cols <- setdiff(names(final_clusters_df), ordered_known_cols); final_col_order <- c(ordered_known_cols, remaining_cols)
    final_clusters_df <- final_clusters_df[, final_col_order]
    output_csv_name <- paste0(output_prefix, "_best_h", round(best_h,4), "_clusters.csv")
    cat("Writing final cluster data to:", output_csv_name, "\n")
    tryCatch({
        final_clusters_df <- final_clusters_df %>% mutate(across(any_of(c("cluster", "seq_length")), as.integer)) %>% mutate(across(any_of(c("cluster_identity", "silhouette_width")), as.numeric))
        write.csv(final_clusters_df, output_csv_name, row.names = FALSE, quote = TRUE, na = "")
    }, error = function(e) { stop("Error writing output CSV file: ", e$message) })
    final_cluster_sizes_plot <- final_clusters_df %>% filter(!is.na(sequence) & sequence != "") %>% count(cluster, name = "size")
    if(nrow(final_cluster_sizes_plot) > 0){
        hist_plot_file <- paste0(output_prefix, "_best_h", round(best_h,4), "_histogram.png"); hist_plot_final <- ggplot(final_cluster_sizes_plot, aes(x = size)) + geom_histogram(binwidth = 1, fill = "steelblue", color = "black") + labs(title = paste("Cluster Sizes (Taxa w/ Seq) at Optimal h =", round(best_h,4)), subtitle = paste("Total clusters:", nrow(final_cluster_sizes_plot), "| Total taxa:", sum(final_cluster_sizes_plot$size)), x = "Cluster Size", y = "Frequency") + theme_bw() + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)); tryCatch(ggsave(hist_plot_file, plot = hist_plot_final, width = 8, height = 6), error=function(e){cat("Warning: Could not save final histogram plot:", e$message, "\n")}); cat("Final cluster size histogram saved to:", hist_plot_file, "\n")
    } else { cat("No clusters with sequences for best h.\n") }
    cat("Generating plot of average silhouette width per cluster...\n")
    if ("silhouette_width" %in% names(final_clusters_df) && nrow(final_clusters_df %>% filter(!is.na(silhouette_width))) > 0) {
        cluster_summary_stats_final <- final_clusters_df %>% filter(!is.na(silhouette_width)) %>% group_by(cluster) %>% summarise(avg_cluster_silhouette = mean(silhouette_width, na.rm = TRUE), cluster_size = n(), .groups = 'drop') %>% filter(cluster_size > 0)
        if(nrow(cluster_summary_stats_final) > 0) {
            cluster_summary_stats_final$cluster <- factor(cluster_summary_stats_final$cluster); cluster_summary_stats_final <- cluster_summary_stats_final %>% arrange(desc(avg_cluster_silhouette)) %>% mutate(cluster = factor(cluster, levels = unique(cluster)))
            plot_cluster_sil_file <- paste0(output_prefix, "_best_h", round(best_h,4), "_cluster_silhouette.png")
            p_cluster_sil_final <- ggplot(cluster_summary_stats_final, aes(x = cluster, y = avg_cluster_silhouette, fill = cluster_size)) + geom_col(color = "black") + scale_fill_viridis_c(option = "plasma", name = "Cluster Size") + labs(title = paste("Average Silhouette Width per Cluster at Optimal h =", round(best_h,4)), subtitle = paste("Total clusters plotted:", nrow(cluster_summary_stats_final)), x = "Cluster ID", y = "Average Silhouette Width") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6))
            tryCatch(ggsave(plot_cluster_sil_file, plot = p_cluster_sil_final, width = 10, height = 6), error=function(e){cat("Warning: Could not save cluster silhouette plot:", e$message, "\n")}); cat("Cluster silhouette plot saved to:", plot_cluster_sil_file, "\n")
        } else { cat("No valid cluster summary stats for silhouette plot.\n") }
    } else { cat("Silhouette widths not available for final data.\n") }
} else { cat("No best 'h' value found by optimization or optimization failed.\n") }
cat("Script finished successfully.\n")
if (exists("best_score") && !is.null(best_score) && !is.na(best_score)) {
  cat(sprintf("FINAL_BEST_SCORE: %.10f\n", best_score))
} else { cat("FINAL_BEST_SCORE: NA\n") }
