#  Author: Jonathan Samson
#  Description: Runs Sleuth and produces graphs and tables

rm(list=ls())
.libPaths("~/thesis/rlib")
setwd("/mnt/scratch/samso008/Project103470/kallisto")
design.file = "../data/design.txt" # Change me

source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
install.packages("Rcpp")  # This was required to get devtools to install
install.packages("devtools")
install.packages("scales")  # Necessary for pheatmap to install
install.packages("tidyselect")  # Necessary for dplyr
install.packages("htmltools")  # Necessary for shiny
install.packages("plyr")  # Necessary for sleuth
install.packages("pkgconfig")  # Sleuth needed a newer version
install.packages("reshape2")  # Sleuth failed without this
install.packages("ggfortify")
install.packages("stringi")  # Needed for tidyverse
install.packages("tidyverse", dependencies = TURE)  # For making the TF only heatmap
install.packages("purrr")  # Needed for tidyverse
install.packages("gplots")
devtools::install_github("pachterlab/sleuth")
devtools::install_cran("tidyverse", force=TRUE)

library(tidyverse)
library(sleuth)
library(data.table)  # Needed to filter with %like%
library(ggfortify)
library(pheatmap)
library(gplots)
library(DESeq2)  # For the clustering function


importDesign <- function(sample.symbol, base.dir = ".",
		         design.file = "design.txt"){
	# Imports the sample design to a list.
	#
	# Args:
	#   sample.symbol: (character) A symbol that differentiates the samples 
	#			       from other folders.
	#   base.dir: (character) The directory that holds the samples, usually WD.
	#   design.file: (character, /path/to/file) file with the design of the 
	#					    project, usually design.txt.
	#
	#  Returns:
	#   design: the list of samples in the design.
	sample.ids <- dir(base.dir, sample.symbol)
	design <- read.table(file.path(base.dir, design.file), header = TRUE,
			     stringsAsFactors = FALSE)
  design$path <- paste0("./", design$sample)
	return(design)
}

subsetDesign <- function(design, filter.str){
	# Filters the samples containing filter.str from design.
	#
	# Args:
	#   design: (list) the design to be filtered.
	#   filter.str: (character) the pattern contained in the unwanted files.
	#
	# Returns:
	#   short.design: (list) the design without the files containing filter.str.
	short.design <- subset(design, !(design$sample %like% filter.str))
	return(short.design)
}

runSleuth <- function(design, func = function(x) log2(X + 0.5)){
	# Runs the sleuth tools to create a sleuth object with fits and an LRT.
	#
	# Args:
	#   design: (list) The design to run sleuth on.
	#   func: (function) The function to transform the data with to prevent NA.
	#
	# Returns:
	#   sleuth: (sleuth object) The sleuth object that was generated.
	sleuth <- sleuth_prep(design, ~condition, extra_bootstrap_summary=TRUE,
			      transformation_function = func)
	sleuth <- sleuth_fit(sleuth, ~condition, 'full')
	sleuth <- sleuth_fit(sleuth, ~1, 'reduced')
	sleuth <- sleuth_lrt(sleuth, 'reduced', 'full')
	return(sleuth)
}

runSleuthGenome <- function(so, s2c, func = function(x) log2(X + 0.5)){
	# Reuns the sleuth tools to create a sleuth object with genomic alignment.
	#
	# Ags:
	#		so: (sleuth object) The sleuth object to begin with.
	#		s2c: (list) The design of the sleuth object.
	#		func: (function) The function to transform with to prevent NA.
	#
	# Returns:
	#		sleuth: (sleuth object) The sleuth object that was generated.
	sleuth_matrix <- sleuth_to_matrix(so, 'obs_norm', 'tpm') #TPM sum ~= 1million
	transcripts <- rownames(sleuth_matrix)
	t2g <- data.frame(target_id=transcripts, gene_id=substring(transcripts, 1, 9),
			  stringsAsFactors=FALSE)
	sleuth <- sleuth_prep(s2c, ~ condition, target_mapping = t2g,
			      aggregation_column = 'gene_id', gene_mode = TRUE)
	# sleuth = sleuth_fit(sleuth, ~condition, 'full')
	# sleuth = sleuth_fit(sleuth, ~1, 'reduced')
	# sleuth = sleuth_lrt(sleuth, 'reduced', "full")
	return(sleuth)
}

genesToFilter <- function(filt.path, sleuth.object){
	# Find the genes found in filt.path that are in the sleuthObject.
	#
	# Args:
	#   filt.path: (string) Path to the file containing genes to be filtered.
	#   sleuthObject: (sleuth object) the object to filter from.
	#
	# Returns:
	#   filt.short: (sleuth object) the sleuth object filtered of genes from the 
	#				file.
	#
  	# Note that the file containing the genes to be filtered needs a header, and
	# the column with the genes to be filtered must be titled "ID".
	filt <- read.table(file.path(filt.path), header = TRUE,
			   stringsAsFactors = FALSE)
	filt.short <- subset(filt, (filt$ID %in%
				    names(sleuth.object$filter_bool)))
	return(filt.short)
}

makeSleuthTable <- function(sleuth.object){
	# Creates a sleuth table from the sleuth object.
	#
	# Args:
	#   sleuth.object: (sleuth object) The so to make the table from.
	#
	# Returns:
	#   sleuth.table: (sleuth table)  The table with sleuth results.
	sleuth.table <- sleuth_results(sleuth.object, 'reduced:full', 'lrt',
				       show_all=FALSE)
	return(sleuth.table)
}

makeSleuthMatrix <- function(sleuth.object){
	# Creates a sleuth matrix from the sleuth object.
	#
	# Args:
	#   sleuth.object: (slueth object) The sleuth object to generate the matrix
	#				   from.
	#
	# Returns:
	#   sleuth.matrix: (sleuth matrix) The sleuth matrix generated from the
	#				   object.
	sleuth.matrix <- t(sleuth_to_matrix(sleuth.object, 'obs_norm', 'tpm'))
	return(sleuth.matrix)
}

makeColumns <- function(sleuth.df, reg.ex = "\\w{3,4}"){
	# Generates a column to add to a dataframe.
	#
	# Args:
	#   sleuth.df: (dataframe) A sleuth dataframe.
	#   reg.ex: (string) The regular expression from the row names with 
	#		     column names.  Default is for sample names in base data.
	# Returns:
	#   result: (list) The column for sample names.
	rgx.result <- regexpr(reg.ex, rownames(sleuth.df))
	result <- regmatches(rownames(sleuth.df), rgx.result)
	return(result)
}

writePCA <- function(sleuth.matrix, sleuth.df, PC.x = 1, PC.y = 2){
	# Creates a PCA.
	# 
	# Args:
	#   sleuth.matrix: (sleuth matrix) The matrix to use to create the PCA.
	#   sleuth.df: (dataframe) The dataframe to use for the data of the PCA.
	#   PC.x: (int) The principal component to use for the X axis, default 1.
	#   PC.y: (int) The principal component to use for the Y axis, defualt 2.
	#
	# Returns:
	#   PCA: (autoplot) The PCA generated by autoplot.
	PCA <- autoplot(prcomp(sleuth.matrix, center = TRUE, scale = FALSE), 
			data = sleuth.df, x = PC.x, y = PC.y, shape = "sample", 
			colour = "condition",
			main = paste("Unscaled PCA for", 
				     deparse(substitute(sleuth.matrix)),
				     "with PCs ", PC.x, "and", PC.y))
	return(PCA)
}

writeScaledPCA <- function(sleuth.matrix, sleuth.df, PC.x = 1, PC.y = 2){
	# Creates a scaled PCA.
	# 
	# Args:
	#   sleuth.matrix: (sleuth matrix) The matrix to use to create the PCA.
	#   sleuth.df: (dataframe) The dataframe to use for the data of the PCA.
	#   PC.x: (int) The principal component to use for the X axis, default 1.
	#   PC.y: (int) The principal component to use for the Y axis, defualt 2.
	#
	# Returns:
	#   PCA: (autoplot) The scaled PCA generated by autoplot
	PCA <- autoplot(prcomp(sleuth.matrix[, apply(sleuth.matrix, 2, var) != 0],
			       center = TRUE, scale = TRUE), data = sleuth.df, x = PC.x,
			       y = PC.y, shape = "sample", colour = "condition", 
			       main = paste("Scaled PCA for", 
			       deparse(substitute(sleuth.matrix)), "with PCs ", PC.x, "and",
			       PC.y))
	return(PCA)
}

makeColumnSortable <- function(df.col){
	# Adds "0"s to the beginning of condition names so they sort properly.
	#
	# Args:
	#   df.col: (list) A dataframe$column.
	#
	# Returns:
	#   sortable.df.col: (list) A sortable version of df.col.
	#
	# Prepends a zero for each: single digits, days, Mock.
  sortable.df.col <- ifelse(grepl("^\\d{1}(d|h)", df.col) == 
			    TRUE, paste0("0", df.col), df.col)
  sortable.df.col <- ifelse(grepl("\\d{2}h", sortable.df.col) == TRUE,
			    paste0("0", sortable.df.col),
			    sortable.df.col)
  sortable.df.col <- ifelse(grepl("\\_M", sortable.df.col) == TRUE,
			    paste0("0", sortable.df.col),
			    sortable.df.col)
  return(sortable.df.col)
}

writeHeatmap <- function(sleuth.matrix, sleuth.df){
	# Creates a heatmap sorted based on conditions.
	#
	# Args:
	#   sleuth.matrix: (sleuth matrix) The matrix to use to create the heatmap.
	#   sleuth.df: (dataframe) A dataframe with sortable columns.
	#
	# Returns:
	#   heat: (heatmap) The heatmap.
	#
	# Note: Running this function on all genes (in Arabidopsis thaliana) has
	# produced segfaults on altschul, and it is not advised to attempt this in
	# the future.  Additionally, running on all genes is not particularly
	# informative, it would be better to focus on only transcription factors.
	heat <- heatmap(sleuth.matrix[rownames(sleuth.matrix[order(sleuth.df$sample,
			sleuth.df$condition), ]), ], scale = 'column', Rowv = NA, 
			labCol = FALSE, main = paste("Heatmap for ",
			deparse(substitute(sleuth.matrix)), " sorted by condition"))
	return(heat)
}

genSleuths <- function(s2c, sig, cond){
	# Creates sleuth object, tables, matrix, and df for the given s2c.
	#
	# Args:
	#   s2c: (list) The design of the sleuths to generate.
	#   sig: (int) The significance value for the Sleuth table.
	#   cond: (character) The regex for selecting conditions in the data frame.
	#
	# Returns:
	#   sleuth.list: (list) A list of the following:
	#     so: (sleuth object) The sleuth object.
	#     st: (sleuth table) Sleuth table.
  #     st.sig: (sleuth table) Sleuth table of significant values.
  #     sm: (sleuth matrix)Sleuth matrix.
	#     sdf: (dataframe) Sleuth data frame.
	so <- runSleuth(s2c)
	st <- makeSleuthTable(so)
	st.sig <- st[st$qval <= sig, ]
	sm <- makeSleuthMatrix(so)
	sdf <- as.data.frame(sm)
	sdf$sample <- makeColumns(sdf)
	sdf$condition <- makeColumns(sdf, cond)
	sdf$condition <- makeColumnSortable(sdf$condition)
	sdf$condition <- factor(sdf$condition, levels = unique(sdf$condition))
	sleuth.list <- list(so, st, st.sig, sm, sdf)
	names(sleuth.list) <- c('so', 'st', 'st.sig', 'sm', 'df')
	return(sleuth.list)
}

genTfFilter <- function(sleuth.list, filter.list, cond){
	# Creates a list of objects that only contain the filter.list genes.
	#
	# Args:
	#   sleuth.list: (list) The list of original objects.
	#   filter.list: (list) The list containing the genes to be filtered.
	#   cond: (character) The regex for selecting conditions in the dataframe.
	#
	# Returns:
	#   filtered.list: (list) The list after filtering, containing:
	#     df.tf: (dataframe) The trasncription factor dataframe.
	#     sm.tf: (sleuth matrix) The transcription factor sleuth matrix.
	#     sm.tf.sig: (sleuth matrix)The significantly different tf sleuth matrix.
	df.tf <- sleuth.list[['df']][, filter.list$ID]
	df.tf$sample <- makeColumns(df.tf)
	df.tf$condition <- makeColumns(df.tf, cond)
	df.tf$condition <- makeColumnSortable(df.tf$condition)
	sm.tf <- sleuth.list[['sm']][, names(sleuth.list[['sm']][1, ]) %in%
				     filter.list$ID]
	sm.tf.sig <- sm.tf[, colnames(sm.tf) %in% sleuth.list[['st.sig']]$target_id]
	filtered.list <- list(sm.tf, df.tf, sm.tf.sig)
	names(filtered.list) <- c('sm.tf', 'df.tf', 'sm.tf.sig')
	return(filtered.list)
}

prepareForPython <- function(dataframe){
	# Prepares the dataframe for use in python.
	#
	# Args:
	#		dataframe: (dataframe) The dataframe to prepare.
	#
	# Returns:
	#		df.prepped: (dataframe) The prepared dataframe.
  dataframe$time = gsub(".*_(\\d+[hd]).*","\\1",rownames(dataframe))
  dataframe$treatment = gsub(".*_","",rownames(dataframe))
  dataframe$stage = paste(dataframe$treatment,dataframe$time)
  levels(dataframe$stage) = c("Mock 4h","Mock 10h","Dex 4h","Dex 10h",
                              "Dex 1d","Dex 2d","Dex 3d","Dex 7d","Dex 14d")
  return(dataframe)
}

# main()
# Assign constants
condition.rgx <- "[[:alnum:]]?\\w[hdo]\\w*"  # Regex for condition from base
sig.q <- 0.05  # Sets a significance value
tf.file <- "/mnt/scratch/samso008/Project103470/reference/Ath_TF_list"


# Design of the experiment is in design.txt loaded with the following
s2c.total <- importDesign("_")
# Creates a sleuth object of the total samples
list.total <- genSleuths(s2c.total, sig.q, condition.rgx)
tf <- genesToFilter(tf.file, list.total[['so']])
list.total <- c(list.total, genTfFilter(list.total, tf, condition.rgx))
sleuth_matrix.total <- sleuth_to_matrix(list.total$so, 'obs_norm', 'tpm') #TPM sum ~= 1million
# sleuth_matrix.total <- prepareForPython(sleuth_matrix.total)
saveRDS((sleuth_matrix.total), "total.sleuth")
# write.csv(sleuth_matrix.total, "total.sleuth", row.names=TRUE)
list.total.sleuth <- runSleuthGenome(list.total$so, s2c.total)
sleuth_matrix.total.genome <- sleuth_to_matrix(list.total.sleuth,
																							 'obs_norm', 'tpm')
# sleuth_matrix.total.genome <- prepareForPython(sleuth_matrix.total.genome)
saveRDS((sleuth_matrix.total.genome), "totalgenome.sleuth")
# write.csv(sleuth_matrix.total.genome, "totalgenome.sleuth", row.names=TRUE)
pdf("total.pdf")
# Generats PCAs for PC's 1+2 and PC's 2+3
writePCA(list.total[['sm']], list.total[['df']])
writePCA(list.total[['sm']], list.total[['df']], 2, 3)
writeScaledPCA(list.total[['sm']], list.total[['df']])
writeScaledPCA(list.total[['sm']], list.total[['df']], 2, 3)
writePCA(list.total[['sm.tf']], list.total[['df.tf']])
writePCA(list.total[['sm.tf']], list.total[['df.tf']], 2, 3)
writeScaledPCA(list.total[['sm.tf']], list.total[['df.tf']])
writeScaledPCA(list.total[['sm.tf']], list.total[['df.tf']], 2, 3)
writeHeatmap(list.total[['sm.tf']], list.total[['df.tf']])
writeHeatmap(list.total[['sm.tf.sig']], list.total[['df.tf']])
dev.off()


# Filters samples
s2c.syn <- subsetDesign(s2c.total, "Han")  # Filters Han samples
s2c.syn <- subsetDesign(s2c.syn, "PLT")  # Filters out PLT
s2c.syn <- subsetDesign(s2c.syn, "Col") # Filters out Col0
# Generates sleuths for the s2c
list.syn <- genSleuths(s2c.syn, sig.q, condition.rgx)
list.syn <- c(list.syn, genTfFilter(list.syn, tf, condition.rgx))
sleuth_matrix.syn <- sleuth_to_matrix(list.syn$so, 'obs_norm', 'tpm') #TPM sum ~= 1million
# sleuth_matrix.syn <- prepareForPython(sleuth_matrix.syn)
saveRDS((sleuth_matrix.syn), "syn.sleuth")
# write.csv(sleuth_matrix.syn, "syn.sleuth", row.names=TRUE)
list.syn.sleuth <- runSleuthGenome(list.syn$so, s2c.syn)
sleuth_matrix.syn.genome <- sleuth_to_matrix(list.syn.sleuth,
																							 'obs_norm', 'tpm')
# sleuth_matrix.syn.genome <- prepareForPython(sleuth_matrix.syn.genome)
saveRDS(sleuth_matrix.syn.genome, "syngenome.sleuth")
# write.csv(sleuth_matrix.syn.genome, "syngenome.sleuth", row.names=TRUE)
pdf("syn.pdf")
# Generats PCAs for PC's 1+2 and PC's 2+3
writePCA(list.syn[['sm']], list.syn[['df']])
writePCA(list.syn[['sm']], list.syn[['df']], 2, 3)
writeScaledPCA(list.syn[['sm']], list.syn[['df']])
writeScaledPCA(list.syn[['sm']], list.syn[['df']], 2, 3)
writePCA(list.syn[['sm.tf']], list.syn[['df.tf']])
writePCA(list.syn[['sm.tf']], list.syn[['df.tf']], 2, 3)
writeScaledPCA(list.syn[['sm.tf']], list.syn[['df.tf']])
writeScaledPCA(list.syn[['sm.tf']], list.syn[['df.tf']], 2, 3)
writeHeatmap(list.syn[['sm.tf']], list.syn[['df.tf']])
writeHeatmap(list.syn[['sm.tf.sig']], list.syn[['df.tf']])
dev.off()


# Filters samples
s2c.base <- subsetDesign(s2c.total, "PLT")  # Filters out PLT
s2c.base <- subsetDesign(s2c.base, "Han") # Filters out HAN
# Generates sleuths for the s2c
list.base <- genSleuths(s2c.base, sig.q, condition.rgx)
list.base <- c(list.base, genTfFilter(list.base, tf, condition.rgx))
sleuth_matrix.base <- sleuth_to_matrix(list.base$so, 'obs_norm', 'tpm') #TPM sum ~= 1million
# sleuth_matrix.base <- prepareForPython(sleuth_matrix.base)
saveRDS((sleuth_matrix.base), "base.sleuth")
# write.csv(sleuth_matrix.base, "base.sleuth", row.names=TRUE)
list.base.sleuth <- runSleuthGenome(list.base$so, s2c.base)
sleuth_matrix.base.genome <- sleuth_to_matrix(list.base.sleuth,
																							 'obs_norm', 'tpm')
# sleuth_matrix.base.genome <- prepareForPython(sleuth_matrix.base.genome)
saveRDS((sleuth_matrix.base.genome), "basegenome.sleuth")
# write.csv(sleuth_matrix.base.genome, "base.genome.sleuth", row.names=TRUE)
pdf("base.pdf")
# Generats PCAs for PC's 1+2 and PC's 2+3
writePCA(list.base[['sm']], list.base[['df']])
writePCA(list.base[['sm']], list.base[['df']], 2, 3)
writeScaledPCA(list.base[['sm']], list.base[['df']])
writeScaledPCA(list.base[['sm']], list.base[['df']], 2, 3)
writePCA(list.base[['sm.tf']], list.base[['df.tf']])
writePCA(list.base[['sm.tf']], list.base[['df.tf']], 2, 3)
writeScaledPCA(list.base[['sm.tf']], list.base[['df.tf']])
writeScaledPCA(list.base[['sm.tf']], list.base[['df.tf']], 2, 3)
writeHeatmap(list.base[['sm.tf']], list.base[['df.tf']])
writeHeatmap(list.base[['sm.tf.sig']], list.base[['df.tf']])
dev.off()


# Filters samples
s2c.411sh <- subsetDesign(s2c.base, "401")  # Filters out 401
# Generates sleuths for the s2c
list.411sh <- genSleuths(s2c.411sh, sig.q, condition.rgx)
list.411sh <- c(list.411sh, genTfFilter(list.411sh, tf, condition.rgx))
sleuth_matrix.411sh <- sleuth_to_matrix(list.411sh$so, 'obs_norm', 'tpm') #TPM sum ~= 1million
# sleuth_matrix.411sh <- prepareForPython(sleuth_matrix.411sh)
saveRDS((sleuth_matrix.411sh), "411sh.sleuth")
# write.csv(sleuth_matrix.411sh, "411sh.sleuth", row.names=TRUE)
list.411sh.sleuth <- runSleuthGenome(list.411sh$so, s2c.411sh)
sleuth_matrix.411sh.genome <- sleuth_to_matrix(list.411sh.sleuth,
																							 'obs_norm', 'tpm')
# sleuth_matrix.411sh.genome <- prepareForPython(sleuth_matrix.411sh.genome)
saveRDS((sleuth_matrix.411sh.genome), "411shgenome.sleuth")
# write.csv(sleuth_matrix.411sh.genome, "411shgenome.sleuth", row.names=TRUE)
pdf("411sh.pdf")
# Generats PCAs for PC's 1+2 and PC's 2+3
writePCA(list.411sh[['sm']], list.411sh[['df']])
writePCA(list.411sh[['sm']], list.411sh[['df']], 2, 3)
writeScaledPCA(list.411sh[['sm']], list.411sh[['df']])
writeScaledPCA(list.411sh[['sm']], list.411sh[['df']], 2, 3)
writePCA(list.411sh[['sm.tf']], list.411sh[['df.tf']])
writePCA(list.411sh[['sm.tf']], list.411sh[['df.tf']], 2, 3)
writeScaledPCA(list.411sh[['sm.tf']], list.411sh[['df.tf']])
writeScaledPCA(list.411sh[['sm.tf']], list.411sh[['df.tf']], 2, 3)
writeHeatmap(list.411sh[['sm.tf']], list.411sh[['df.tf']])
writeHeatmap(list.411sh[['sm.tf.sig']], list.411sh[['df.tf']])
dev.off()


# Filters samples
s2c.411 <- subsetDesign(s2c.411sh, "Col")  # Filters out Col
# Generates sleuths for the s2c
list.411 <- genSleuths(s2c.411, sig.q, condition.rgx)
list.411 <- c(list.411, genTfFilter(list.411, tf, condition.rgx))
sleuth_matrix.411 <- sleuth_to_matrix(list.411$so, 'obs_norm', 'tpm') #TPM sum ~= 1million
# sleuth_matrix.411 <- prepareForPython(sleuth_matrix.411)
saveRDS((sleuth_matrix.411), "411.sleuth")
# write.csv(sleuth_matrix.411, "411.sleuth", row.names=TRUE)
list.411.sleuth <- runSleuthGenome(list.411$so, s2c.411)
sleuth_matrix.411.genome <- sleuth_to_matrix(list.411.sleuth,
																							 'obs_norm', 'tpm')
# sleuth_matrix.411.genome <- prepareForPython(sleuth_matrix.411.genome)
saveRDS((sleuth_matrix.411.genome), "411genome.sleuth")
# write.csv(sleuth_matrix.411.genome, "411genome.sleuth", row.names=TRUE)
pdf("411.pdf")
# Generats PCAs for PC's 1+2 and PC's 2+3
writePCA(list.411[['sm']], list.411[['df']])
writePCA(list.411[['sm']], list.411[['df']], 2, 3)
writeScaledPCA(list.411[['sm']], list.411[['df']])
writeScaledPCA(list.411[['sm']], list.411[['df']], 2, 3)
writePCA(list.411[['sm.tf']], list.411[['df.tf']])
writePCA(list.411[['sm.tf']], list.411[['df.tf']], 2, 3)
writeScaledPCA(list.411[['sm.tf']], list.411[['df.tf']])
writeScaledPCA(list.411[['sm.tf']], list.411[['df.tf']], 2, 3)
writeHeatmap(list.411[['sm.tf']], list.411[['df.tf']])
writeHeatmap(list.411[['sm.tf.sig']], list.411[['df.tf']])
dev.off()


# Filters samples
s2c.401sh <- subsetDesign(s2c.base, "411")  # Filters out 411
# Generates sleuths for the s2c
list.401sh <- genSleuths(s2c.401sh, sig.q, condition.rgx)
list.401sh <- c(list.401sh, genTfFilter(list.401sh, tf, condition.rgx))
sleuth_matrix.401sh <- sleuth_to_matrix(list.401sh$so, 'obs_norm', 'tpm') #TPM sum ~= 1million
# sleuth_matrix.401sh <- prepareForPython(sleuth_matrix.401.sh)
saveRDS((sleuth_matrix.401sh), "401sh.sleuth")
# write.csv(sleuth_matrix.401sh, "401sh.sleuth", row.names=TRUE)
list.401sh.sleuth <- runSleuthGenome(list.401sh$so, s2c.401sh)
sleuth_matrix.401sh.genome <- sleuth_to_matrix(list.401sh.sleuth,
																							 'obs_norm', 'tpm')
# sleuth_matrix.401sh.genome <- prepareForPython(sleuth_matrix.401sh.genome)
saveRDS((sleuth_matrix.401sh.genome), "401shgenome.sleuth")
# write.csv(sleuth_matrix.401sh.genome, "401shgenome.sleuth", row.names=TRUE)
pdf("401sh.pdf")
# Generats PCAs for PC's 1+2 and PC's 2+3
writePCA(list.401sh[['sm']], list.401sh[['df']])
writePCA(list.401sh[['sm']], list.401sh[['df']], 2, 3)
writeScaledPCA(list.401sh[['sm']], list.401sh[['df']])
writeScaledPCA(list.401sh[['sm']], list.401sh[['df']], 2, 3)
writePCA(list.401sh[['sm.tf']], list.401sh[['df.tf']])
writePCA(list.401sh[['sm.tf']], list.401sh[['df.tf']], 2, 3)
writeScaledPCA(list.401sh[['sm.tf']], list.401sh[['df.tf']])
writeScaledPCA(list.401sh[['sm.tf']], list.401sh[['df.tf']], 2, 3)
writeHeatmap(list.401sh[['sm.tf']], list.401sh[['df.tf']])
writeHeatmap(list.401sh[['sm.tf.sig']], list.401sh[['df.tf']])
dev.off()


# Filters samples
s2c.401 <- subsetDesign(s2c.401sh, "Col")  # Filters out Col
# Generates sleuths for the s2c
list.401 <- genSleuths(s2c.401, sig.q, condition.rgx)
list.401 <- c(list.401, genTfFilter(list.401, tf, condition.rgx))
sleuth_matrix.401 <- sleuth_to_matrix(list.401$so, 'obs_norm', 'tpm') #TPM sum ~= 1million
# sleuth_matrix.401 <- prepareForPython(sleuth_matrix.401)
saveRDS((sleuth_matrix.401), "401.sleuth")
# write.csv(sleuth_matrix.401, "401.sleuth", row.names=TRUE)
list.401.sleuth <- runSleuthGenome(list.401$so, s2c.401)
sleuth_matrix.401.genome <- sleuth_to_matrix(list.401.sleuth,
																							 'obs_norm', 'tpm')
# sleuth_matrix.401.genome <- prepareForPython(sleuth_matrix.401.genome)
saveRDS((sleuth_matrix.401.genome), "401genome.sleuth")
# write.csv(sleuth_matrix.401.genome, "401genome.sleuth", row.names=TRUE)
pdf("401.pdf")
# Generats PCAs for PC's 1+2 and PC's 2+3
writePCA(list.401[['sm']], list.401[['df']])
writePCA(list.401[['sm']], list.401[['df']], 2, 3)
writeScaledPCA(list.401[['sm']], list.401[['df']])
writeScaledPCA(list.401[['sm']], list.401[['df']], 2, 3)
writePCA(list.401[['sm.tf']], list.401[['df.tf']])
writePCA(list.401[['sm.tf']], list.401[['df.tf']], 2, 3)
writeScaledPCA(list.401[['sm.tf']], list.401[['df.tf']])
writeScaledPCA(list.401[['sm.tf']], list.401[['df.tf']], 2, 3)
writeHeatmap(list.401[['sm.tf']], list.401[['df.tf']])
writeHeatmap(list.401[['sm.tf.sig']], list.401[['df.tf']])
dev.off()


# Separate out different conditions from s2c
# Separate out different conditions from s2c
condA  <-  which(s2c$condition == "411_0004h_Mock")
condB  <-  which(s2c$condition == "411_004h")
s2c.AvsB  <-  s2c[c(condA,condB),]
###############################################################################
# perform Wald test which returns a value that is analogous to the fold-change
so <- sleuth_wt(so,which_beta="conditionCol0_Shoot", which_model="full")
so <- sleuth_wt(so,which_beta="condition401_4h_Mock", which_model="full")
so <- sleuth_wt(so,which_beta="condition401_10h_Mock", which_model="full")
so <- sleuth_wt(so,which_beta="condition411_4h_Mock", which_model="full")
so <- sleuth_wt(so,which_beta="condition411_10h_Mock", which_model="full")
so <- sleuth_wt(so,which_beta="condition401_4h", which_model="full")
so <- sleuth_wt(so,which_beta="condition401_10h", which_model="full")
so <- sleuth_wt(so,which_beta="condition401_1d", which_model="full")
so <- sleuth_wt(so,which_beta="condition401_2d", which_model="full")
so <- sleuth_wt(so,which_beta="condition401_3d", which_model="full")
so <- sleuth_wt(so,which_beta="condition401_7d", which_model="full")
so <- sleuth_wt(so,which_beta="condition401_10d", which_model="full")
so <- sleuth_wt(so,which_beta="condition401_14d", which_model="full")
so <- sleuth_wt(so,which_beta="condition411_4h", which_model="full")
so <- sleuth_wt(so,which_beta="condition411_10h", which_model="full")
so <- sleuth_wt(so,which_beta="condition411_1d", which_model="full")
so <- sleuth_wt(so,which_beta="condition411_2d", which_model="full")
so <- sleuth_wt(so,which_beta="condition411_3d", which_model="full")
so <- sleuth_wt(so,which_beta="condition411_7d", which_model="full")
so <- sleuth_wt(so,which_beta="condition411_14d", which_model="full")

# Plots a sample heatmap using Jensen-Shannon divergence
# This map is not clustered
plot_sample_heatmap(so, cluster_bool = FALSE, filename =
		    "unclusteredheatmap.pdf", show_rownames = FALSE,
		    show_colnames = FALSE, annotation_row = , annotation_col = ,
		    annotation_colors = )

# Get the results of the Wald test
sleuth.table  <-  sleuth_results(so,"conditionCol0_Shoot", "wt", show_all=FALSE)

# Group density, conditions should be similar
plot_group_density(so, use_filtered = TRUE, units = "est_counts",trans="log",
		   grouping="condition", offset = 1)

# Checks loadings for PC's
pdf("loadings.pdf")
layout(matrix(1:3,nrow=1))
plot_loadings(so, pc_input = 1)
plot_loadings(so, pc_input = 2)
plot_loadings(so, pc_input = 3)
dev.off()

# Plots bootstraps involved in pc1 loadings
pdf("PC1bootstraps.pdf")
layout(matrix(1:5,nrow=1), widths = lcm(50), heights = lcm(50))
plot_bootstrap(so, "AT3G09260.1", x_axis_angle = 0)  # PYK1, long ER body
plot_bootstrap(so, "AT1G21310.1", x_axis_angle = 0)  # EXT3 Extensin 3
plot_bootstrap(so, "AT5G02500.1", x_axis_angle = 0)  # Heat Shock Cognate 70
plot_bootstrap(so, "AT1G29930.1", x_axis_angle = 0)  # CAB1 Chlorophyll AB 1
plot_bootstrap(so, "AT1G66200.1", x_axis_angle = 0)  # GLN Glutamine Synthase
dev.off()

# Plots bootstraps involved in pc2 loadings
pdf("PC2bootstraps.pdf")
layout(matrix(1:5,nrow=1))
plot_bootstrap(so, "AT1G29930.1", x_axis_angle = 0)  # CAB1 Chlorophyll AB 1
plot_bootstrap(so, "AT2G34420.1", x_axis_angle = 0)  # PSII LHC B1B2
plot_bootstrap(so, "AT1G67090.1", x_axis_angle = 0)  # Ribulose Biphosphate
plot_bootstrap(so, "AT3G09260.1", x_axis_angle = 0)  # PYK1 long ER body
plot_bootstrap(so, "AT1G29910.1", x_axis_angle = 0)  # LH Chlorophyll A/B 1.2
dev.off()

# Plots bootstraps involved in pc3 loadings
pdf("PC3bootstraps.pdf")
layout(matrix(1:5,nrow=1))
plot_bootstrap(so, "AT2G09260.1", x_axis_angle = 0) # PYK1 long ER body
plot_bootstrap(so, "AT2G21160.1", x_axis_angle = 0) # TRAP translocon 
plot_bootstrap(so, "AT5G02500.1", x_axis_angle = 0) # Heat shock cognate 70
plot_bootstrap(so, "AT3G16640.1", x_axis_angle = 0) # TCTP1
plot_bootstrap(so, "AT1G21310.1", x_axis_angle = 0) # EXT3 Extensin 3
dev.off()

# Check the expression of endogenous vs induced genes
pdf("endoinducedexp.pdf")
plot_bootstrap(list.base[['so']], "ATSG37650.1", x_axis_angle = 0) # SHR
plot_bootstrap(list.base[['so']], "AT4G37650.1", x_axis_angle = 0)
plot_bootstrap(list.base[['so']], "ATSG54220.1", x_axis_angle = 0) # SCR
plot_bootstrap(list.base[['so']], "AT3G54220.1", x_axis_angle = 0)
plot_bootstrap(list.base[['so']], "ATSG20840.1", x_axis_angle = 0) # PLT1
plot_bootstrap(list.base[['so']], "AT3G20840.1", x_axis_angle = 0)
plot_bootstrap(list.base[['so']], "ATSG17430.1", x_axis_angle = 0) # PLT4
plot_bootstrap(list.base[['so']], "AT5G17430.1", x_axis_angle = 0)
plot_bootstrap(list.base[['so']], "ATSG57390.2", x_axis_angle = 0) # PLT5
plot_bootstrap(list.base[['so']], "AT5G57390.2", x_axis_angle = 0)
plot_bootstrap(list.base[['so']], "ATSG11260.1", x_axis_angle = 0) # Wox5
plot_bootstrap(list.base[['so']], "AT3G11260.1", x_axis_angle = 0)
plot_bootstrap(list.base[['so']], "ATSG78080.1", x_axis_angle = 0) # WIND1
plot_bootstrap(list.base[['so']], "AT1G78080.1", x_axis_angle = 0)
plot_bootstrap(list.base[['so']], "ATSG99995.1", x_axis_angle = 0) # AmiGO
plot_bootstrap(list.base[['so']], "ATSG99996.1", x_axis_angle = 0)  # Bar
plot_bootstrap(list.base[['so']], "ATSG99997.1", x_axis_angle = 0)  # ALC
plot_bootstrap(list.base[['so']], "ATSG99998.1", x_axis_angle = 0)  # XVE
plot_bootstrap(list.base[['so']], "ATSG99999.1", x_axis_angle = 0)  # GVG
dev.off()

pdf("nobootstrap.pdf")
barplot(c(list.base[['sm']][,"ATSG99995.1"]), main = "Synthetic AmiGO", 
				ylab = "Counts", las = 2)
barplot(c(list.base[['sm']][,"ATSG37650.1"]), main="Synthetic SHR",
				ylab = "Counts", las = 2)
barplot(c(list.base[['sm']][,"ATSG54220.1"]), main = "Synthetic SCR", 
				ylab = "Counts", las = 2)
dev.off()

# readies slueth for comparison with DESeq2
sleuth_matrix <- sleuth_to_matrix(so, 'obs_norm', 'tpm') #TPM sum ~= 1million
transcripts <- rownames(sleuth_matrix)
t2g <- data.frame(target_id=transcripts, gene_id=substring(transcripts, 1, 9),
		  stringsAsFactors=FALSE)
so = sleuth_prep(s2c, ~ condition, target_mapping = t2g, aggregation_column =
		 'gene_id')
so = sleuth_fit(so, ~conditionm 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_w(so, which_beta="conditionST", which_model="full")
so = sleuth_wt(so, which_beta="conditionST", which_model="full")
sleuth_table <- sleuth_results(so, 'conditionST', 'wt', show_all=FALSE)

# Launches a web interface to view the results
sleuth_live(so)
