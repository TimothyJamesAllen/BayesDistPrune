#!/usr/bin/env Rscript

# --- Calculate Patristic Distance Matrix and Save as Dist Object (RDS) & CSV ---
#
# Description:
# This script reads a phylogenetic tree, calculates the patristic distance matrix
# (cophenetic distances), converts it to a 'dist' object, and saves the
# 'dist' object to an .rds file. It also saves the full distance matrix to a .csv file.
#
# Usage:
# Rscript distance_matrix.R <input_tree_file>
#
# Arguments:
#   input_tree_file : Path to the phylogenetic tree file (e.g., Newick, Nexus).
#
# Output:
#   1. <input_tree_file>_dist_object.rds : An RDS file containing the 'dist' object.
#   2. <input_tree_file>_distance_matrix.csv : A CSV file of the full distance matrix.
#

# --- Dependencies ---
suppressPackageStartupMessages({
    library(ape)
    library(tictoc)
})

# --- Main Script ---
tic()

# Argument Parsing
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  cat("Usage: Rscript distance_matrix.R <input_tree_file>\n\n")
  cat("Error: Incorrect number of arguments. Please provide the input tree file.\n")
  quit(status = 1)
}

input_file <- args[1]

if (!file.exists(input_file)) {
  cat("Error: Input tree file not found:", input_file, "\n")
  quit(status = 1)
}

cat("Reading tree from:", input_file, "...\n")
tree <- tryCatch({
  read.tree(input_file)
}, error = function(e) {
  cat("Error reading tree file '", input_file, "': ", e$message, "\n")
  return(NULL)
})

if (is.null(tree)) {
  quit(status = 1)
}
cat("Tree successfully read.\n")

cat("Calculating patristic distance matrix (cophenetic distances)...\n")
patristic_matrix <- tryCatch({
  cophenetic.phylo(tree)
}, error = function(e) {
  cat("Error calculating cophenetic distances: ", e$message, "\n")
  return(NULL)
})

if (is.null(patristic_matrix)) {
  quit(status = 1)
}
cat("Patristic distance matrix calculated.\n")

# Convert to dist object and save as RDS
cat("Converting matrix to 'dist' object...\n")
dist_object <- tryCatch({
  as.dist(patristic_matrix)
}, error = function(e) {
  cat("Error converting matrix to dist object: ", e$message, "\n")
  return(NULL)
})

if (is.null(dist_object)) {
  cat("Skipping RDS save due to error in dist object conversion.\n")
} else {
  output_rds_file <- paste0(tools::file_path_sans_ext(basename(input_file)), "_dist_object.rds")
  cat("Saving 'dist' object to RDS:", output_rds_file, "...\n")
  tryCatch({
    saveRDS(dist_object, file = output_rds_file)
    cat("RDS file with 'dist' object saved successfully.\n")
  }, error = function(e) {
    cat("Error writing RDS file '", output_rds_file, "': ", e$message, "\n")
  })
}

toc()
cat("Script finished.\n")
