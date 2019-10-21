#! /usr/bin/R

## Author: Jonathan Samson
## 
## Runs DESeq2 on the input data, and produces relevant graphs, clusters, etc.
##
## col.data requires a phenom.csv file, with 4 columns:
##			names	condition	construct	time
## and requires those columns to be named as above.  The names must be
## the same as the names given by the gene_count_matrix.csv so the script
## can match them up properly.

rm(list=ls())
.libPaths("~/thesis/rlib")
setwd("/mnt/scratch/samso008/Project103470/hisat2/images")
count.data.file <- paste0("/mnt/scratch/samso008/Project103470/",
													"hisat2/gene_count_matrix.csv"),
col.data.file  <- "/mnt/scratch/samso008/Project103470/hisat2/phenom.csv"
tf <- "/mnt/scratch/samso008/Project103470/reference/Ath_TF_list"
goi.file <- "/mnt/scratch/samso008/Project103470/hisat2/images/geneclusters.csv"
# The GOI file was not included in the data files, and must be created.  It is
# expected to be a file with a header which has gene loci in a column.

## ----------------------------------------------------------------------------
source("http://bioconductor.org/biocLite.R")
install.packages("xml2")
install.packages("igraph")
biocLite("DESeq2")
biocLite("graph")
biocLite("RCy3")
biocLite("EBSeq")
biocLite("biomaRt")
biocLite("topGO")
biocLite("blockmodeling")
install.packages("mnormt")
biocLite("multtest")
biocLite("Hmisc")
install.packages(
	"https://cran.r-project.org/src/contrib/Archive/psych/psych_1.7.5.tar.gz",
	repos = NULL, type = "source")

## ----------------------------------------------------------------------------
library(DESeq2)
library(ggplot2)
library(dplyr)
library(data.table)
library(RCy3)
library(igraph)
library("psych")
library(EBSeq)
library("biomaRt")
library("topGO")
library("multtest")

## ----------------------------------------------------------------------------

## ----------------------------------------------------------------------------
## Define functions
x.prefix <- function(xless){
	# Appends a prefix of "X" to the string.
	#
	# Keyword Arguments:
	#   xless: (character) The string to prefix with X.
	#
	# Returns:
	#   xed: (character) The string prefixed with X, Xxless.
	#
	# This is required because rownames can not start with a number.
	xed <- paste0("X", xless)
	return(xed)
}

genDESeq <- function(count.matrix, column){
	# Generates the DESeq data set from.
	# 
	# Keyword Arguments:
	#   count.matrix: (matrix) The matrix read in from the .csv file.
	#   column: (dataframe) The column data.
	# 
	# Returns:
	#   dds: (S4) DESeq Data Set generated from the matrix.
	rownames(column)[grep("^\\d", rownames(column))] <- lapply(rownames(column[
				   grep("^\\d", rownames(column)), ]), x.prefix)
	apply(count.matrix, 2, sum)
	mx <- apply(count.matrix, 1, max)
	count.matrix <- count.matrix[mx > 10, ]
	dds <- DESeqDataSetFromMatrix(count.matrix[, order(colnames(count.matrix))],
																column[order(rownames(column)), ],
																~condition)
	dds <- dds[, colnames(dds[, order(dds$construct, -rank(dds$condition), 
																		dds$time)])]
	dds <- estimateSizeFactors(dds)
	sizeFactors(dds)
	dds <- estimateDispersions(dds)
	dds <- nbinomWaldTest(dds)
	return(dds)
}

genesToFilter <- function(filt.path, sleuth.object){
	# Find the genes found in filt.path that are in the sleuthObject.
	#
	# Args:
	#   filt.path: (character)Path to the file containing genes to be filtered.
	#   sleuth.object: (sleuth object) the object to filter from.
	#
	# Returns:
	#   filt.short: (sleuth object) the sleuth object filtered of genes from the 
	#							  file.
	#
  # Note that the file containing the genes to be filtered needs a header, and
	# the column with the genes to be filtered must be titled "ID".
	filt <- read.table(file.path(filt.path), header = TRUE,
										 stringsAsFactors = FALSE)
	filt.short <- subset(filt, (filt$ID %in%
											 names(sleuth.object$filter_bool)))
	return(filt.short)
}

createDiffList <- function(data.counts, data.columns, tfs, filename){
	# Creates a list with all of the applicable information for the plots.
	#
	# Keyword Arguments:
	#   data.counts: (integer) A matrix of count data.
	#   data.columns: (list) The column names for the count data.
	#   tfs:  (character) A file holding the transcription factor information.
	#		filename: (character) The file to output or read from.
	#
	# Returns:
	#   list.data: A list contining the following:
	#     dds: (S4) A DESeq object.
	#     dds.counts: (double) The counts of the dds.
	#			dds.dist: (double) The pearson correlation of the whole set.
	#			dds.fpm: (double) The fragments per million of the dds.
	#     rld: (S4)The rlog of the dds.
	#     rld.dists: (double) The distance data for the rld.
	#     res: (S4) The results of the dds.
	#     tf: (double) The dds filtered to only include transcription factors.
	#     tf.dist: (double) The pearson correlation of the transcription factors.
	#
	# filename is used to create a file that has the diff list.  The function
	# first checks if the file is present, and if so, reads the diff list from
	# the file rather than recalculating.  If not, it calculates and creates
	# the file.
	if(file.exists(filename)){
		print("DESeq file exists, reading...")
		list.data <- readRDS(file = filename)
		print("Read complete.")
	}
	else{
		print("DESeq file does not exist. Generating list...")
		dds <- genDESeq(data.counts, data.columns)
		dds.counts <- counts(dds, normalized = TRUE)
		dds.dist <- as.dist(0.5 - 0.5 * (cor(t(dds.counts), method = "pearson")), 
												diag = TRUE)
		dds.fpm <- fpm(dds)
		rld <- rlog(dds)
		rld.dists <- dist(t(assay(rld)))
		res <- results(dds)
		res$padj <- ifelse(is.na(res$padj), 1, res$padj)
		tf.table <- read.table(tfs, header = TRUE, row.names = 1)
		tf <- (dds.counts[tf.table$Gene.ID, ])
		tf.dist <- as.dist(0.5 - 0.5 * (cor(t(tf), method = "pearson")), diag = TRUE)
		list.data <- list(dds, dds.counts, dds.fpm, dds.dist, rld, rld.dists, res, 
											tf, tf.dist)
		names(list.data) <- c("dds", "dds.counts", "dds.fpm", "dds.dist", "rld", 
													"rld.dists", "res", "tf", "tf.dist")
		print("List generation complete.  Writing to file...")
		saveRDS(list.data, file=filename)
		print("Write complete.")
	}
	return(list.data)
}

genHC <- function(dists, cutoff){
	# Creates a hierachical cluster.
	#
	# Keyword Arguments:
	#   dists: (double) The pearson correlation.
	#   cutoff: (float) The cutoff for clusters.
	#
	# Returns:
  #   hc: (hclust) The hierarchical clustering table
	print("Creating HC tree")
	hc.tree <- hclust(dists, method = "complete")
	print("Applying cutoff to tree")
	hc <- cutree(hc.tree, h = cutoff)
	table(hc)
	print("Tree complete")
	return(hc)
}

genClusterPatterns <- function(dds.dist, dds.fpm, cutoff, prefix, goi) {
	# Generates pattern files with both average only and all samples.
	#
	# Keyword Arguments:
	#   dds.dist: (double) The distance for the wanted dds.
	#		dds.fpm: (double) The fragments per million for the wanted dds.
	#   cutoff: (float) The cutoff for clustering.
	#   prefix: (character) The prefix for the file names.
	#		goi:	(list) The genes of interest to plot.
	#
	# Returns:
	#	  clusterfilename: (character) The filename for the cluster .txt file.
	print("Calling genHC")
	hc <- genHC(dds.dist, cutoff)
	print("Cluster created")
	print(length(table(hc)))
	clusterFileName <- createClusters(dds.fpm, cutoff, prefix, hc)
	print("Cluster file created")
	genClusterPatternFile(dds.fpm, cutoff, prefix, "average", hc)
	print("Average pattern created")
	genClusterPatternFile(dds.fpm, cutoff, prefix, "all", hc)
	print("All pattern created")
	genClusterPatternFile(dds.fpm, cutoff, prefix, "genotype", hc)
	print("Genotype pattern created")
	plotGOIClusters(dds.fpm, cutoff, prefix, hc, goi)
	print("GOI pattern created")
	findGOTerms(dds.fpm, cutoff, prefix, hc)
	return(clusterFileName)
}
	
createClusters <- function(dds.fpm, cutoff, prefix, list.hc){
	# Creates a cluster .txt file for the given dds.
	#
	# Keyword Arguments:
	#   dds.fpm: (double) The CPM of the dds to be clustered.
	#   cutoff: (float) The cutoff for clustering.
	#   prefix: (character) The prefix for the file names.
	#		list.hc: (hclust) The hierarchical cluster for the list, from genHC.
	#
	# Returns:
	#	  clusterfilename: The filename for the cluster .txt file
  clusterfilename = paste0(prefix, "clusters-", cutoff, ".txt")
  sink(clusterfilename)
  for (loop in as.integer(names(sort(table(list.hc), 
			                               decr = T)))) {
  	ids <- names(list.hc)[list.hc == loop]
		df.new <- createClusterDataframe(dds.fpm, ids)
  	if(nrow(df.new) <= 1){
  		break
  	}
  	cat("cluster:", loop, paste(as.character(sort(unique(df.new$gene))), 
  			collapse = ","), "\n")
	}
	sink()
	return(clusterfilename)
}

addDataframeColumn <- function(regex, column.name, dataframe) {
	# Adds a column to the dataframe based on the sample.
	#
	# Keyword Arguments:
	#		regex: (chatracter) The regular expression to use to find the column 
	#					 within the dataframe$sample.
	#		column.name: (character) The name of the column.
	#		dataframe: (dataframe) The dataframe to add the column to.
	#
	# Returns:
	#		dataframe: (dataframe) The dataframe after the column has been added.
	regex.result <- regexpr(regex, dataframe$sample)
	matches <- regmatches(dataframe$sample, regex.result)
	dataframe[[column.name]] <- matches
	dataframe[[column.name]] <- factor(dataframe[[column.name]],
																		 levels = unique(dataframe[[column.name]]))
	return(dataframe)
}

createClusterDataframe <- function(dds.fpm, ids) {
	# Creates a dataframe for a single cluster.
	#
	# Keyword Arguments:
	#		dds.fpm: (double) The fragments per million of the appropriate dataset.
	#		ids: (list) The gene IDS (loci) of the genes in the cluster.
	#
	# Returns:
	#		cluster.df: (dataframe) The dataframe for the cluster.
  d.f <- data.frame(unique(dds.fpm[rownames(dds.fpm) %in% ids, ]))
  df.scaled <- apply(d.f, 1, scale)
  names(df.scaled) <- names(d.f)
	colnames(df.scaled) <- rownames(d.f)
  rownames(df.scaled) <- colnames(d.f)
  df.new <- as.data.frame(as.table(as.matrix(df.scaled)))
  names(df.new) <- c("sample", "gene", "z_score")
	df.new <- addDataframeColumn("^[[:alnum:]]*_", "genotype", df.new)
	df.new <- addDataframeColumn("[[:alnum:]]*_[[:alnum:]]*$", "cond", df.new)
}

genClusterPatternFile <- function(dds.fpm, cutoff, prefix, aes.group, 
																	list.hc) {
	# Creates a .pdf of the clustering patterns for the given dds.
	#
	# Keyword Arguments:
	#   dds.fpm: (double) The fragments per million of the dds to be clustered.
	#   cutoff: (float) The cutoff for clustering.
	#   prefix: (character) The prefix for the file names.
	#   aes.group: (character) The type of aesthetic grouping (genotype, average, all).
	#		list.hc: (hclust) A hierarchical cluster of the dds.fpm from genHC.
  par(mfrow = c(1, 1), mgp = c(3, 0.5, 0), mar = c(5, 3, 0.75, 0.5), tck =
  		-0.025, omi = c(0.1, 0.5, 0.5, 0.1))
  par(las = 2)
  pdf(file = paste0(prefix, aes.group, "patterns_", cutoff, "-", Sys.Date(),
			".pdf"), 12, 8)
	print("Entering loop")
  for (loop in as.integer(names(sort(table(list.hc), 
			                               decr = T)))) {
  	ids <- names(list.hc)[list.hc == loop]
		df.new <- createClusterDataframe(dds.fpm, ids)
		print("Data frame created, checking for >1 member")
  	if(ncol(df.new) <= 1){
  		break
  	}

		sample.condition <- paste0(df.new$genotype, df.new$cond)
		df.new$condition <- sample.condition
		df.new$condition <- factor(df.new$condition, 
															 levels = unique(df.new$condition))

		print("Factoring complete, setting plot.x type")
		if (aes.group == "average"){
			plot.group <- 1
			plot.x <- df.new$condition
			group.means <- df.new %>% 
				group_by(condition) %>%
				summarize(mean_z = mean(z_score))
			df.new$zmeans <- group.means$mean_z[
														group.means$condition[df.new$condition]]
		}
		else if (aes.group == "genotype") {
			plot.group <- df.new$genotype
			plot.x <- df.new$cond
			group.means <- df.new %>% 
				group_by(cond) %>%
				summarize(mean_z = mean(z_score))
			df.new$zmeans <- group.means$mean_z[group.means$cond[df.new$cond]]
		}
	  else {
			plot.group <- df.new$gene
			plot.x <- df.new$condition
			group.means <- df.new %>% 
				group_by(condition) %>%
				summarize(mean_z = mean(z_score))
			df.new$zmeans <-
				group.means$mean_z[group.means$condition[df.new$condition]]
		}

		print("Printing ggplot to file")
  	print(ggplot(df.new, aes(x = plot.x, sample, y = z_score, 
														 group = plot.group, color = genotype)) +
					geom_point() + 
					ylim(-5, 8) +
  				theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  				labs(title = paste(prefix, " cluster", loop, ", size:", length(ids))) + 
  				theme(legend.position = "none") +
					geom_line(aes(y = zmeans, group = plot.group)))
		print("Print complete, restarting loop")
  }
  dev.off()
}

findGOTerms <- function(dds.fpm, cutoff, prefix, list.hc) {
	# Finds the GO terms of the cluster and prints them to a file
	#
	# Keyword Arguments:
	#   dds.fpm: The CPM of the dds to be analyzed
	#   cutoff: The cutoff for clustering
	#   prefix: The prefix for the file names
	#		list.hc: A hierarchical cluster of the dds.fpm from genHC
	filename <- paste0(prefix, "clusters-", cutoff, "GO.txt")
	print("Creating file and entering loop")
	sink(filename)
	options("width"=10000)
	mart <- biomaRt::useMart(biomart = "plants_mart", 
														dataset = "athaliana_eg_gene",
														host = "plants.ensembl.org")
	GTOGO <- biomaRt::getBM(attributes = c("ensembl_gene_id", "go_id"), 
													mart = mart)
	GTOGO <- GTOGO[GTOGO$go_id != '', ]  # Removes blank entries
	geneID2GO <- by(GTOGO$go_id,
									GTOGO$ensembl_gene_id,
									function(x) as.character(x))
	all.genes <- sort(unique(as.character(GTOGO$ensembl_gene_id)))
	all.genes <- all.genes[all.genes %in% rownames(dds.fpm)]
  for (loop in as.integer(names(sort(table(list.hc), 
			                               decr = T)))) {
  	ids <- names(list.hc)[list.hc == loop]
		# df.new <- createClusterDataframe(dds.fpm, ids)
		int.genes <- factor(as.integer(all.genes %in% ids))
		names(int.genes) <- all.genes
		
		go.obj <- new("topGOdata", ontology = 'BP',
									allGenes = int.genes,
									annot = annFUN.gene2GO,
									gene2GO = geneID2GO)
		result <- runTest(go.obj, algorithm = "elim", statistic = "fisher")
		results.tab <- GenTable(object = go.obj, elimFisher = result, topNodes = 10)

		#Create a map of geneIDs to GO terms
		ann.genes <- genesInTerm(go.obj)
		
		# Get multiple testing correction values (In progress)
		results.tab['multTestCorr'] = p.adjust(results.tab$elimFisher, method="fdr")

		#Select a few GO terms from the Fisher analysis
		fisher.go <- names(sort(score(result)))[1:10]
		fisher.ann.genes <- genesInTerm(go.obj, whichGO=fisher.go)
		fisher.ann.genes
		fisher.ann.genes <- data.frame(GO.ID = names(fisher.ann.genes), loci =
																	 I(fisher.ann.genes))
		results.tab = merge(x = results.tab, y = fisher.ann.genes)

		### Use getPvalues() to adjust.
		print(paste("Cluster:", loop, " Number of genes:", length(ids)))
		for (i in 1:10){
			print(results.tab[i, 1:7])
			print(ids[ids %in% results.tab[[i, 8]]])
		}
	}
	sink()
}

plotGOIClusters <- function(dds.fpm, cutoff, prefix, list.hc, 
														interesting.genes) {
	# Creates a .pdf of the clustering patterns for the given dds.
	#
	# Keyword Arguments:
	#   dds.fpm: (double) The CPM of the dds to be clustered.
	#   cutoff: (float) The cutoff for clustering.
	#   prefix: (character) The prefix for the file names.
	#		list.hc: (hclust) A hierarchical cluster of the dds.list from genHC.
	#		interesting.genes: (list) The genes to plot.
  par(mfrow = c(1, 1), mgp = c(3, 0.5, 0), mar = c(5, 3, 0.75, 0.5), tck =
  		-0.025, omi = c(0.1, 0.5, 0.5, 0.1))
  par(las = 2)
  pdf(file = paste0(prefix, "goi_patterns_", cutoff, "-", Sys.Date(),
			".pdf"), 12, 8)
	print("Entering loop")
  for (loop in as.integer(names(sort(table(list.hc), 
			                               decr = T)))) {
  	ids <- names(list.hc)[list.hc == loop]
		df.new <- createClusterDataframe(dds.fpm, ids)
		for (geneofinterest in interesting.genes$Locus) {
			if (geneofinterest %in% df.new$gene) {
				group.means <- df.new %>% 
					group_by(cond) %>%
					summarize(mean_z = mean(z_score))
				df.new$zmeans <- group.means$mean_z[group.means$cond[df.new$cond]]

				average.goi <- subset(df.new, gene %in% geneofinterest) %>%
					group_by(cond) %>%
					summarize(mean_gene = mean(z_score))
				df.new$genemeans <-
					average.goi$mean_gene[average.goi$cond[df.new$cond]]

  			print(ggplot(df.new, aes(x = cond, sample, y = z_score, 
																 group = genotype, color = genotype)) +
							geom_point() + 
 							geom_point(data = subset(df.new, gene %in% geneofinterest),
													color = "blue", size = 2) +
							ylim(-5, 8) +
  						theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  						labs(title = paste(prefix, 
																 interesting.genes$Name[interesting.genes$Locus
																												%in% geneofinterest],
																 "size:", length(ids))) + 
  						theme(legend.position = "none") +
							geom_line(aes(y = zmeans), size=1) + 
							geom_line(data = subset(df.new, gene %in% geneofinterest),
												aes(y = genemeans), color = "blue", size = 1.5))
			}
		}
  }
  dev.off()
}

genPDF <- function(dds.list, list.title) {
	# Creates a .pdf with all of the relevant graphs for the dds.list.
	#
	# Keyword Arguments:
	#   dds.list: (list) The list with the dds data, as generated by
	#							createDiffList.
	#   list.title: (character) The title for the pdf.
	pdf(paste0(list.title, ".pdf"))
	plot(density(assay(dds.list$dds)[, 1]), main = "counts")
	plot(density(assay(dds.list$rld)[, 1]), main = "log counts")
	plot(hclust(dds.list$rld.dists))
	plotDispEsts(dds.list$dds)
	plotMA(dds.list$res, main = "MA Plot", ylim = c(-20, 20), alpha = 0.01)
	writeHeatmap(dds.list$tf) # Parameter may need to be t(dds.list$tf)
	dev.off()
}

writeHeatmap <- function(dds.matrix) {
	# Creates a heatmap sorted based on conditions.
	#
	# Keyword Arguments:
	#   dds.matrix: (matrix) The matrix to use to create the heatmap.
	#
	# Returns:
	#   heat: (heatmap) the heatmap.
	#
	# Note: Issues have been found when attempting to create a heatmap on the
	# whole sample size with the sleuth pipeline, so do not attempt.
	heat <- heatmap(dds.matrix[rownames(dds.matrix), ], scale = 'column', 
									Rowv = NA, labCol = FALSE, main = paste("Heatmap for ",
									deparse(substitute(dds.matrix)), " sorted by condition"))
	return(heat)
}

makeNetworkFiles <- function(cur.list, counting.data, column.data) {
	# Creates files for use in Cytoscape.
	#
	# Keyword Argumets:
	#   cur.list: (list) The list generated by createDiffList.
	#   counting.data: (integer) The raw count data for this list.
	#   column.data: (list) The col.data for this list.
	#
	# From "Step-by-Step Construction of Gene Co-expression Networks
	# from High-Throughput Arabidopsis RNA Sequencing Data",
	# Contreras-Lopez et al, 2018.
	res <- na.exclude(as.data.frame(cur.list$res))
	filtered <- res[(abs(res$log2FoldChange)>1 & res$padj<0.01), ]
	write.table(filtered, 
							paste0(deparse(substitute(cur.list)), "regulated.txt"), 
							quote = F, sep = "\t", col.names = NA)
	norm.data <- getNormalizedMat(counting.data, MedianNorm(counting.data))
	norm.data.log <- log2(norm.data+1)
	norm.interest <-norm.data.log[rownames(filtered), ]
	norm.interest.corr <- corr.test(t(norm.interest), method = "pearson", ci = F)
	norm.interest.corr$p[lower.tri(norm.interest.corr$p, diag = TRUE)] = NA
	pval.adj <- as.data.frame(as.table(norm.interest.corr$p))
	norm.interest.corr$r[lower.tri(norm.interest.corr$r, diag = TRUE)] = NA
	correlation <- as.data.frame(as.table(norm.interest.corr$r))
	cor.table <- na.exclude(cbind(correlation, pval.adj))[, c(1,2,3,6)]
	colnames(cor.table) <- c("gene1", "gene2", "cor", "p.adj")
	cor.table.filt <- cor.table[(abs(cor.table[, 3]) > 0.9 & 
															 cor.table[, 4] < 0.01), ]
	p.adj <- cor.table.filt[, 4]
	p.adj[p.adj == 0] <- as.numeric(unlist(format(.Machine)))[1]
	cor.table.filt <- cbind(cor.table.filt, log10(p.adj))
	write.table(cor.table.filt, 
							paste0(deparse(substitute(cur.list)), "cor.table.filter.txt"), 
							sep = "\t", row.names = F, quote = F)
	g <- graph.data.frame(cor.table.filt[, 1:2], directed = FALSE)
	degrees <- degree(g)
	betweennesses <- betweenness(g)
	node.nw.st <- data.frame(degrees, betweennesses)
	rank.stat <- rowMeans(cbind(rank(node.nw.st[, 1]), rank(node.nw.st[, 2])))
	node.nw.st <- cbind(node.nw.st, rank.stat)
	write.table(node.nw.st, 
							file = paste0(deparse(substitute(cur.list)), "node.nw.st.txt"), 
							sep = "\t", col.names = NA, quote = F)
}

confirmPipelines <- function(sleuthData, deseqFPM) {
	# Runs pearson correlation to check that the two pipelines are comparable.
	#
	# Keyword Arguments:
	#		sleuthData: (dataframe) The dataframe generated by Sleuth.
	#		deseqFPM: (double) The Fragments per Million generated by DESeq.
	# Returns:
	#		corr: (rcorr) The correlation and p value generated by Hmisc's rcorr.
	merged <- merge(sleuthData, deseqFPM, by="row.names", all=F)
	midway = as.integer(length(colnames(merged))/2 + 1)
	corr <- Hmisc::rcorr(unlist(merged[2:midway]),
											 unlist(merged[(midway+1):length(colnames(merged))]), 
											 type="pearson")
	return(corr)
}

## ----------------------------------------------------------------------------
count.data  <- as.matrix(read.csv(count.data.file), row.names="gene_id")
col.data  <- read.csv(col.data.file, sep="\t", row.names=1)
goi <- read.table(goi.file, header=TRUE)

# Total Data Set (401 + 411 + Col0 Shoot + PLT + Han)
# total <- createDiffList(count.data, col.data, tf.file)
# totalclusterfilename <- genClusterPatterns(total$dds.dist, total$dds.fpm, ct,
# 																					 "total")
# genPDF(total, "total")
# makeNetworkFiles(total, count.data, col.data)

# Base Data Set (401 + 411 + Col0 Shoot)
base.count.data <- (count.data[, !colnames(count.data) %like% "Han"])
base.count.data <- (base.count.data[, !colnames(base.count.data) %like% "PLT"])
base.col.data <- (col.data[!rownames(col.data) %like% "Han", ])
base.col.data <- (base.col.data[!rownames(base.col.data) %like% "PLT", ])
base <- createDiffList(base.count.data, base.col.data, tf.file, "base.deseq")
ct <- 0.645
baseclusterfilename <- genClusterPatterns(base$dds.dist, base$dds.fpm, ct, 
																					"base", goi)
ct <- 0.45
basetfclusterfilename <- genClusterPatterns(base$tf.dist, base$tf, ct,
																						"basetf", goi)
genPDF(base, "base")
makeNetworkFiles(base, base.count.data, base.col.data)

#411 Data
count.data.411 <- (count.data[, colnames(count.data) %like% "411"])
col.data.411 <- (col.data[rownames(col.data) %like% "411", ])
x411 <- createDiffList(count.data.411, col.data.411, tf.file, "x411.deseq")
ct <- 0.679
x411clusterfilename <- genClusterPatterns(x411$dds.dist, x411$dds.fpm, ct, 
																					"411", goi)
ct <- 0.447
x411tfclusterfilename <- genClusterPatterns(x411$tf.dist, x411$tf, ct, "411tf",
																						goi)
genPDF(x411, "411")
makeNetworkFiles(x411, count.data.411, col.data.411)

# 401 Data
count.data.401 <- (count.data[, colnames(count.data) %like% "401"])
col.data.401 <- (col.data[rownames(col.data) %like% "401", ])
x401 <- createDiffList(count.data.401, col.data.401, tf.file, "x401.deseq")
ct <- 0.664
x401clusterfilename <- genClusterPatterns(x401$dds.dist, x401$dds.fpm, ct, 
																					"401", goi)
ct <- .42
x401tfclusterfilename <- genClusterPatterns(x401$tf.dist, x401$tf, ct, "401tf",
																						goi)
genPDF(x401, "401")
makeNetworkFiles(x401, count.data.401, col.data.401)

# 411 + Col0 Shoot Data
count.data.411shoot <- (count.data[, !colnames(count.data) %like% "Han"])
count.data.411shoot <- (count.data.411shoot[, !colnames(count.data.411shoot) %like%
											 	"PLT"])
count.data.411shoot <- (count.data.411shoot[, !colnames(count.data.411shoot) %like%
												"401"])
col.data.411shoot <- (col.data[!rownames(col.data) %like% "Han", ])
col.data.411shoot <- (col.data.411shoot[!rownames(col.data.411shoot) %like%
											 	"PLT", ])
col.data.411shoot <- (col.data.411shoot[!rownames(col.data.411shoot) %like% 
											"401", ])
x411shoot <- createDiffList(count.data.411shoot, col.data.411shoot, tf.file,
														"x411shoot.deseq")
ct <- 0.666
x411shootclusterfilename <- genClusterPatterns(x411shoot$dds.dist,
																							 x411shoot$dds.fpm, 
																							 ct, "411shoot", goi)
ct <- 0.442
x411shoottfclusterfilename <- genClusterPatterns(x411shoot$tf.dist,
																								 x411shoot$tf, ct,
																								 "411shoottf", goi)
genPDF(x411shoot, "411shoot")
makeNetworkFiles(x411shoot, count.data.411shoot, col.data.411shoot)

# 401 + Col0 Shoot Data
count.data.401shoot <- (count.data[, !colnames(count.data) %like% "Han"])
count.data.401shoot <- (count.data.401shoot[, !colnames(count.data.401shoot) %like%
											 	"PLT"])
count.data.401shoot <- (count.data.401shoot[, !colnames(count.data.401shoot) %like%
												"411"])
col.data.401shoot <- (col.data[!rownames(col.data) %like% "Han", ])
col.data.401shoot <- (col.data.401shoot[!rownames(col.data.401shoot) %like% 
											"PLT", ])
col.data.401shoot <- (col.data.401shoot[!rownames(col.data.401shoot) %like% 
											"411", ])
x401shoot <- createDiffList(count.data.401shoot, col.data.401shoot, tf.file,
														"x401shoot.deseq")
ct <- 0.649
x401shootclusterfilename <- genClusterPatterns(x401shoot$dds.dist,
																							 x401shoot$dds.fpm, 
																							 ct, "401shoot", goi)
ct <- 0.416
x401shoottfclusterfilename <- genClusterPatterns(x401shoot$tf.dist,
																									x401shoot$tf, ct,
																									"401shoottf", goi)
genPDF(x401shoot, "401shoot")
makeNetworkFiles(x401shoot, count.data.401shoot, col.data.401shoot)

# Experimental Data Set (401 + 411)
count.data.exp <- (count.data[, !colnames(count.data) %like% "Han"])
count.data.exp <- (count.data.exp[, !colnames(count.data.exp) %like% "PLT"])
count.data.exp <- (count.data.exp[, !colnames(count.data.exp) %like% "Col"])
col.data.exp <- (col.data[!rownames(col.data) %like% "Han", ])
col.data.exp <- (col.data.exp[!rownames(col.data.exp) %like% "PLT", ])
col.data.exp <- (col.data.exp[!rownames(col.data.exp) %like% "Col", ])
experim <- createDiffList(count.data.exp, col.data.exp, tf.file,
													"experim.deseq")
ct <- 0.649
experimclusterfilename <- genClusterPatterns(experim$dds.dist, experim$dds.fpm, 
																						 ct, "exp", goi)
ct <- 0.477
experimtfclusterfilename <- genClusterPatterns(experim$tf.dist, experim$tf, ct,
																								"exptf", goi)
genPDF(experim, "exp")
makeNetworkFiles(experim, count.data.exp, col.data.exp)

# Imports sleuth data and compares pipelines using Pearson Correlation.
sleuth401 <- readRDS("~/dash/401genome.sleuth")
corr401 = confirmPipelines(sleuth401, x401$dds.fpm)
sleuth411 <- readRDS("~/dash/411genome.sleuth")
corr411 = confirmPipelines(sleuth411, x411$dds.fpm)
sleuth401shoot <- readRDS("~/dash/401shgenome.sleuth")
corr401shoot = confirmPipelines(sleuth401shoot, x401shoot$dds.fpm)
sleuth411shoot <- readRDS("~/dash/411shgenome.sleuth")
corr411shoot = confirmPipelines(sleuth411shoot, x411shoot$dds.fpm)
sleuthbase <- readRDS("~/dash/basegenome.sleuth")
corrbase = confirmPipelines(sleuthbase, base$dds.fpm)
sleuthexperim <- readRDS("~/dash/syngenome.sleuth")
correxperim = confirmPipelines(sleuthexperim, experim$dds.fpm)
## ----------------------------------------------------------------------------
