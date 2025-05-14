#!/usr/bin/env Rscript

# --- Parallel Bayesian Optimization for K-Medoids (PAM) using Phylogenetic Subtree Partitions ---
#
# Description:
# This script partitions a large dataset by identifying disjoint clades of a specified
# size range from a main phylogenetic tree. It then runs 'bayes_kmedoids.r' (which uses PAM)
# in parallel for each partition and finally combines the clustering results.
#
# Usage:
# Rscript run_parallel_bayes_kmedoids.r \
#   <full_dist_matrix_rds> \
#   <full_tree_file> \
#   <full_fasta_file> \
#   <final_output_prefix> \
#   <min_subtree_tips> \
#   <max_subtree_tips> \
#   <path_to_bayes_kmedoids_script> \
#   <k_min> <k_max> <init_pts> <n_iter> \
#   <med_min> <med_max> <mean_target> <min_msa_size> \
#   <w_median> <w_mean> <w_identity> <w_silhouette> \
#   <num_cores_for_bayes_kmedoids> \
#   <num_parallel_partitions>
#
# Arguments:
#   1. full_dist_matrix_rds       : Path to the FULL distance matrix object (.rds).
#   2. full_tree_file             : Path to the FULL phylogenetic tree file.
#   3. full_fasta_file            : Path to the FULL (cleaned) sequence FASTA file.
#   4. final_output_prefix        : Prefix for the final combined output files.
#   5. min_subtree_tips           : Minimum number of tips for a subtree partition.
#   6. max_subtree_tips           : Maximum number of tips for a subtree partition.
#   7. path_to_bayes_kmedoids_script : Full path to the 'bayes_kmedoids.r' script.
#   --- bayes_kmedoids.r arguments (8-21) ---
#   8. k_min                      : Min k for bayes_kmedoids.r
#   9. k_max                      : Max k for bayes_kmedoids.r
#   10. init_pts                  : BO Init Points for bayes_kmedoids.r
#   11. n_iter                    : BO Iterations for bayes_kmedoids.r
#   12. med_min                   : Target Median Min Size
#   13. med_max                   : Target Median Max Size
#   14. mean_target               : Target Mean Size
#   15. min_msa_size              : Min Cluster Size for MSA
#   16. w_median                  : Weight Median
#   17. w_mean                    : Weight Mean
#   18. w_identity                : Weight Identity
#   19. w_silhouette              : Weight Silhouette
#   20. num_cores_for_bayes_kmedoids : Num cores for EACH bayes_kmedoids.r instance.
#   --- Wrapper script argument ---
#   21. num_parallel_partitions   : Number of partitions to process in parallel by this wrapper.
#
# Example:
# Rscript run_parallel_bayes_kmedoids.r \
#   full_dist.rds full_tree.tre full_seqs_cleaned.fasta \
#   final_run_kmedoids_phylo 4000 6000 path/to/bayes_kmedoids.r \
#   50 500 10 20 4 10 8 2 0.3 0.2 0.4 0.1 2 \
#   3

# --- Dependencies ---
cat("Loading required packages for wrapper...\n")
suppressPackageStartupMessages({
    library(ape)
    library(Biostrings)
    library(dplyr)
    library(stringr)
    library(parallel)
    library(tools) # For file_path_sans_ext
    library(tictoc)
    library(rlang)
    # Consider adding phytools if getDescendants is preferred and available
})

tic("Total wrapper script time")

# --- Argument Parsing for Wrapper ---
args <- commandArgs(trailingOnly = TRUE)
expected_args <- 21 # Updated number of arguments
if (length(args) != expected_args) {
  cat("Usage: Rscript run_parallel_bayes_kmedoids.r <full_dist_matrix_rds> <full_tree_file> <full_fasta_file> <final_output_prefix> <min_subtree_tips> <max_subtree_tips> <path_to_bayes_kmedoids_script> ... (13 bayes_kmedoids.r args) ... <num_parallel_partitions>\n")
  cat("Error: Incorrect number of arguments. Expected", expected_args, "but received", length(args), "\n")
  quit(status = 1)
}

full_dist_matrix_rds_file      <- args[1]
full_tree_file                 <- args[2]
full_fasta_file                <- args[3]
final_output_prefix            <- args[4]
min_subtree_tips               <- as.integer(args[5])
max_subtree_tips               <- as.integer(args[6])
path_to_bayes_kmedoids_script  <- args[7]
bayes_kmedoids_args_vector     <- args[8:20] 
num_parallel_partitions        <- as.integer(args[21])


cat("--- Wrapper Configuration ---\n")
cat("Full Distance Matrix RDS :", full_dist_matrix_rds_file, "\n")
cat("Full Tree File           :", full_tree_file, "\n")
cat("Full FASTA File          :", full_fasta_file, "\n")
cat("Final Output Prefix      :", final_output_prefix, "\n")
cat("Min Subtree Tips         :", min_subtree_tips, "\n")
cat("Max Subtree Tips         :", max_subtree_tips, "\n")
cat("Path to bayes_kmedoids.r :", path_to_bayes_kmedoids_script, "\n")
cat("Num Parallel Partitions  :", num_parallel_partitions, "\n")
cat("Arguments for bayes_kmedoids.r:\n", paste(bayes_kmedoids_args_vector, collapse=" "), "\n")
cat("---------------------------\n")

if (!file.exists(full_dist_matrix_rds_file)) stop("Full distance matrix RDS not found: ", full_dist_matrix_rds_file)
if (!file.exists(full_tree_file)) stop("Full tree file not found: ", full_tree_file)
if (!file.exists(full_fasta_file)) stop("Full FASTA file not found: ", full_fasta_file)
if (!file.exists(path_to_bayes_kmedoids_script)) stop("bayes_kmedoids.r script not found at: ", path_to_bayes_kmedoids_script)
if (min_subtree_tips < 2) stop("min_subtree_tips must be at least 2.")
if (max_subtree_tips < min_subtree_tips) stop("max_subtree_tips must be >= min_subtree_tips.")
if (num_parallel_partitions < 1) stop("Number of parallel partitions must be at least 1.")

temp_dir <- file.path(dirname(final_output_prefix), paste0(basename(final_output_prefix), "_temp_partitions_", format(Sys.time(), "%Y%m%d%H%M%S")))
if (!dir.exists(temp_dir)) {
  dir.create(temp_dir, recursive = TRUE)
  cat("Created temporary directory for partition data:", temp_dir, "\n")
}

cat("Loading full distance matrix...\n")
full_dist_obj <- readRDS(full_dist_matrix_rds_file)
if (inherits(full_dist_obj, "dist")) {
    full_dist_matrix <- as.matrix(full_dist_obj)
} else if (is.matrix(full_dist_obj)) {
    full_dist_matrix <- full_dist_obj
} else {
    stop("Loaded distance data is not a 'dist' object or a matrix.")
}
full_dist_labels <- rownames(full_dist_matrix)
if (is.null(full_dist_labels)) stop("Full distance matrix must have rownames.")

cat("Loading full FASTA sequences...\n")
full_seqs <- readAAStringSet(full_fasta_file)
full_seq_names <- names(full_seqs)

cat("Loading full tree...\n")
full_tree <- read.tree(full_tree_file)
if (!setequal(full_tree$tip.label, full_dist_labels)) {
    warning("Tip labels in the tree do not exactly match labels in the distance matrix. This may cause issues. Proceeding with tree labels for partitioning.")
}

# --- 2. Partitioning using Phylogenetic Subtrees (More Efficiently) ---
cat("DEBUG: Entered Partitioning section.\n")
cat("Partitioning data using phylogenetic clades (min:", min_subtree_tips, ", max:", max_subtree_tips, " tips)...\n")
tic("Data partitioning (phylogenetic clades)")

# Helper function to get tips for a clade defined by a node
get_clade_tips <- function(tree, node) {
  if (node <= length(tree$tip.label)) { # It's a tip
    return(tree$tip.label[node])
  }
  # Ensure phangorn is explicitly loaded or its functions are namespaced if not automatically.
  # The requireNamespace check should handle availability.
  descendants <- phangorn::Descendants(tree, node, type = "tips")[[1]] 
  return(tree$tip.label[descendants])
}
cat("DEBUG: Defined get_clade_tips function.\n")

# Ensure phangorn is loaded for Descendants
cat("DEBUG: Checking for phangorn package using requireNamespace...\n")
if (!requireNamespace("phangorn", quietly = TRUE)) {
    cat("ERROR: requireNamespace(\"phangorn\") returned FALSE. Package 'phangorn' might not be installed correctly or not found in library paths.\n")
    stop("Package 'phangorn' is needed for the Descendants function. Please install it.", call. = FALSE)
}
cat("DEBUG: requireNamespace(\"phangorn\") was successful.\n")

# Attempt to load phangorn to be certain its functions are available
cat("DEBUG: Attempting to load phangorn via library(phangorn)...\n")
if (!library(phangorn, logical.return = TRUE, quietly = TRUE)) {
    cat("ERROR: library(phangorn) failed. The package is installed but cannot be loaded. Check for DLL issues or missing dependencies like igraph.\n")
    stop("Failed to load 'phangorn'. Please check its installation and dependencies.", call. = FALSE)
}
cat("DEBUG: phangorn loaded successfully via library().\n")


cat("DEBUG: Defining internal_nodes...\n")
internal_nodes <- (length(full_tree$tip.label) + 1):(length(full_tree$tip.label) + full_tree$Nnode)
cat("DEBUG: internal_nodes defined. Count:", length(internal_nodes), "\n")

cat("DEBUG: Initializing potential_clades list...\n")
potential_clades <- list()
cat("DEBUG: potential_clades list initialized.\n")

cat("Identifying potential clades from internal nodes...\n") # Existing message
node_counter <- 0
progress_interval <- 500 # Print progress every 500 nodes

for (node in internal_nodes) {
  node_counter <- node_counter + 1
  if (node_counter %% progress_interval == 0) {
    cat("  Processed", node_counter, "of", length(internal_nodes), "internal nodes for potential clades...\n")
  }
  clade_tips <- get_clade_tips(full_tree, node)
  num_clade_tips <- length(clade_tips)
  if (num_clade_tips >= min_subtree_tips && num_clade_tips <= max_subtree_tips) {
    potential_clades[[as.character(node)]] <- list(tips = clade_tips, size = num_clade_tips, node_id = node)
  }
}
cat("Finished identifying potential clades. Found", length(potential_clades), "potential clades within the size range.\n")

# Sort potential clades by size (descending)
if (length(potential_clades) > 0) {
    potential_clades <- potential_clades[order(sapply(potential_clades, `[[`, "size"), decreasing = TRUE)]
}

tip_partitions <- list()
assigned_tips <- character(0)
partition_count <- 0

cat("Selecting disjoint clades...\n")
for (clade_info in potential_clades) {
  if (all(!clade_info$tips %in% assigned_tips)) {
    partition_count <- partition_count + 1
    partition_name <- paste0("clade_N", clade_info$node_id, "_P", partition_count)
    tip_partitions[[partition_name]] <- clade_info$tips
    assigned_tips <- c(assigned_tips, clade_info$tips)
    cat("Selected partition", partition_name, "with", clade_info$size, "tips.\n")
  }
}

remaining_tips <- setdiff(full_tree$tip.label, assigned_tips)
if (length(remaining_tips) > 0) {
  if (length(remaining_tips) >= min_subtree_tips) {
    partition_count <- partition_count + 1
    partition_name <- paste0("remainder_P", partition_count)
    tip_partitions[[partition_name]] <- remaining_tips
    cat("Created a remainder partition", partition_name, "with", length(remaining_tips), "tips.\n")
  } else {
    cat("Skipping remainder partition as it has only", length(remaining_tips), "tips (less than min_subtree_tips:", min_subtree_tips, "). These tips will not be clustered by this script.\n")
  }
}

toc() 
cat("Data partitioned into", length(tip_partitions), "actual groups based on phylogeny.\n")
if (length(tip_partitions) == 0) {
    stop("No suitable disjoint clades found with the given size constraints. Adjust min/max_subtree_tips or check your tree structure.")
}
for(i in seq_along(tip_partitions)) {
    cat("Partition", names(tip_partitions)[i], "size:", length(tip_partitions[[i]]), "tips.\n")
}

# --- 3. Function to Process a Single Partition ---
process_partition_kmedoids <- function(partition_name, tip_list, common_args) { # Changed partition_idx to partition_name
  partition_label <- partition_name # Use the actual name (e.g., clade_N123_P1)
  num_tips_in_partition <- length(tip_list)
  cat("Processing Partition:", partition_label, "with", num_tips_in_partition, "tips.\n")

  k_min_for_sub_problem <- as.integer(common_args$bayes_kmedoids_args[1]) 
  min_practical_partition_size <- max(10, k_min_for_sub_problem + 5) 
  min_size_to_process <- max(min_practical_partition_size, common_args$min_subtree_tips_arg)

  if (num_tips_in_partition < min_size_to_process) {
    cat(partition_label, ": Skipping processing as it has only", num_tips_in_partition, "tips (less than min practical/specified size of", min_size_to_process, ").\n")
    return(NULL)
  }
  
  part_dist_rds_file <- file.path(common_args$temp_dir, paste0(partition_label, "_dist.rds"))
  part_tree_file     <- file.path(common_args$temp_dir, paste0(partition_label, "_tree.tre"))
  part_fasta_file    <- file.path(common_args$temp_dir, paste0(partition_label, "_seqs.fasta"))
  part_output_prefix <- file.path(common_args$temp_dir, paste0(basename(common_args$final_output_prefix), "_", partition_label))
  
  cat(partition_label, ": Creating sub-distance matrix...\n")
  sub_dist_matrix <- common_args$full_dist_matrix[tip_list, tip_list]
  saveRDS(as.dist(sub_dist_matrix), file = part_dist_rds_file)
  
  cat(partition_label, ": Creating sub-FASTA...\n")
  sub_seqs <- common_args$full_seqs[names(common_args$full_seqs) %in% tip_list]
  writeXStringSet(sub_seqs, filepath = part_fasta_file)
  
  cat(partition_label, ": Creating sub-tree...\n")
  sub_tree <- keep.tip(common_args$full_tree, tip_list)
  write.tree(sub_tree, file = part_tree_file)
  
  cmd_args <- c(
      part_dist_rds_file,
      part_tree_file,
      part_fasta_file,
      part_output_prefix,
      common_args$bayes_kmedoids_args 
  )
  
  cat(partition_label, ": Running bayes_kmedoids.r with args:", paste(cmd_args, collapse=" "), "\n")
  log_file <- file.path(common_args$temp_dir, paste0(partition_label, "_bayes_kmedoids.log"))
  
  tic(paste(partition_label, "bayes_kmedoids.r execution"))
  tryCatch({
    system2(
      "Rscript",
      args = c(common_args$path_to_bayes_kmedoids_script, cmd_args),
      stdout = log_file,
      stderr = log_file
    )
    cat(partition_label, ": bayes_kmedoids.r finished. Log in:", log_file, "\n")
  }, error = function(e) {
    cat(partition_label, ": ERROR running bayes_kmedoids.r. Check log:", log_file, "\nError message:", e$message, "\n")
  })
  toc_val <- toc(quiet=TRUE)
  cat(sprintf("%s: bayes_kmedoids.r execution took %.3f sec\n", partition_label, toc_val$toc - toc_val$tic))
  
  bayes_opt_rds_path <- paste0(part_output_prefix, "_bayesopt_results.rds")
  if(file.exists(bayes_opt_rds_path)){
      bayes_opt_data <- readRDS(bayes_opt_rds_path)
      best_k_val <- max(2, round(getBestPars(bayes_opt_data)$k_cont)) 
      results_csv_path <- paste0(part_output_prefix, "_best_k", best_k_val, "_clusters.csv")
      if (file.exists(results_csv_path)) {
          cat(partition_label, ": Reading results from", results_csv_path, "\n")
          df <- read.csv(results_csv_path, stringsAsFactors = FALSE)
          df$original_partition_id <- partition_label 
          return(df)
      } else {
          cat(partition_label, ": Results CSV not found at", results_csv_path, "\n")
          return(NULL)
      }
  } else {
      cat(partition_label, ": Bayes Opt RDS not found at", bayes_opt_rds_path, "\n")
      return(NULL)
  }
}

# --- 4. Run Processing for Each Partition in Parallel ---
cat("\n--- Starting Parallel Processing of Partitions ---\n")
common_args_for_parallel <- list(
  temp_dir = temp_dir,
  full_dist_matrix = full_dist_matrix,
  full_seqs = full_seqs,
  full_tree = full_tree,
  final_output_prefix = final_output_prefix, 
  path_to_bayes_kmedoids_script = path_to_bayes_kmedoids_script,
  bayes_kmedoids_args = bayes_kmedoids_args_vector,
  min_subtree_tips_arg = min_subtree_tips 
)

partition_names_to_process <- names(tip_partitions) # Get the actual names

if (num_parallel_partitions > 1 && .Platform$OS.type != "windows") {
  cat("Using mclapply for parallel execution on", num_parallel_partitions, "cores.\n")
  all_partition_results <- mclapply(
    partition_names_to_process, # Iterate over names
    function(p_name) process_partition_kmedoids(p_name, tip_partitions[[p_name]], common_args_for_parallel),
    mc.cores = num_parallel_partitions
  )
} else {
  if (.Platform$OS.type == "windows" && num_parallel_partitions > 1) {
    cat("Windows detected. Using parLapply for parallel execution with", num_parallel_partitions, "cores.\n")
    cl <- makeCluster(num_parallel_partitions)
    clusterExport(cl, varlist = c("process_partition_kmedoids", "common_args_for_parallel", "tip_partitions"), envir = environment())
    clusterEvalQ(cl, { library(ape); library(Biostrings); library(dplyr); library(stringr); library(tools); library(tictoc); library(cluster); library(DECIPHER); library(ParBayesianOptimization); library(phangorn) })
    
    all_partition_results <- parLapply(
      cl,
      partition_names_to_process, # Iterate over names
      function(p_name) process_partition_kmedoids(p_name, tip_partitions[[p_name]], common_args_for_parallel)
    )
    stopCluster(cl)
  } else {
    cat("Running partitions sequentially (num_parallel_partitions=1 or non-Unix single core).\n")
    all_partition_results <- lapply(
      partition_names_to_process, # Iterate over names
      function(p_name) process_partition_kmedoids(p_name, tip_partitions[[p_name]], common_args_for_parallel)
    )
  }
}

# --- 5. Combine Results ---
cat("\n--- Combining Results from Partitions ---\n")
valid_results <- Filter(Negate(is.null), all_partition_results)

if (length(valid_results) == 0) {
  stop("No valid results obtained from any partition. Check logs in:", temp_dir)
}

global_cluster_counter <- 0
combined_df_list <- list()

for (i in seq_along(valid_results)) {
  partition_df <- valid_results[[i]]
  if (nrow(partition_df) > 0 && "cluster" %in% names(partition_df)) {
    unique_local_clusters <- unique(partition_df$cluster)
    local_to_global_map <- setNames(seq(global_cluster_counter + 1, global_cluster_counter + length(unique_local_clusters)), unique_local_clusters)
    
    partition_df$global_cluster_id <- local_to_global_map[as.character(partition_df$cluster)]
    global_cluster_counter <- global_cluster_counter + length(unique_local_clusters)
    combined_df_list[[i]] <- partition_df
  }
}

if (length(combined_df_list) > 0) {
    final_combined_df <- do.call(rbind, combined_df_list)
    if (nrow(final_combined_df) > 0) {
      final_csv_path <- paste0(final_output_prefix, "_combined_clusters.csv")
      cat("Writing final combined cluster data to:", final_csv_path, "\n")
      write.csv(final_combined_df, final_csv_path, row.names = FALSE, quote = TRUE, na = "")
    } else {
      cat("No data to combine into a final CSV after processing partitions.\n")
    }
} else {
    cat("No valid partition results to combine.\n")
}

cat("Wrapper script finished. Check temporary files in:", temp_dir, "if issues occurred.\n")
toc()
