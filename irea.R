library(Seurat)
library(plyr)
library(dplyr)
library(reshape2)
library(qusage)
library(Hmisc)
library(pheatmap)
library(XLConnect)
library(ggplot2)
library(knitr)
library(xlsx)
library(RColorBrewer)
library(pdist)
library(philentropy)
set.seed(0)






#### Immune response enrichment score functions #####

#' Get enrichment score
#'
#' \code{GetEnrichmentScore} Compute enrichment score using the pathway overrepresentation test
#' which is commonly involves the hypergeometric test
#'
#' @param    degs              A list of differentially expresssed genes that user would like to investigate
#' @param    input_celltype          Choose from one of the listed cell types that most resemble the input
#'
#'

GeneSetEnrichmentScore = function(degs, input_celltype) {
  `%notin%` = Negate(`%in%`)
  celltypes = c("B_cell", "cDC1", "cDC2", "eTAC", "ILC", "Macrophage", "MigDC",
                "Monocyte", "Neutrophil", "NK_cell", "pDC", "T_cell_CD4", "T_cell_CD8",
                "T_cell_gd", "Treg")
  if (input_celltype %notin% celltypes) stop(paste0("cell type must be one of the following: ", paste0(celltypes, collapse = ", ")))

  # Load reference data
  ## TODO: to speed up, can save each cell type into a different Rds file
  ## TODO: to speed up, subset the gene list to only include the genes that are >0.01 (or other threshold) between cytokine and PBS
  load("~/Dropbox (Personal)/Hacohen/ligands/rdata/201209-ligands-alldata-seurat-p3-subsample-100cells.RData")
  lig_seurat = subset(lig_seurat, celltype == input_celltype)
  cytokines = setdiff(sort(unique(lig_seurat@meta.data$sample)), "PBS")
  profiles_cc = as.matrix(lig_seurat@assays[['RNA']][,])

  # only include the DEGs that overlap between user input and reference data
  degs_to_include = degs[degs %in% rownames(profiles_cc)]
  profiles_cc_genes = profiles_cc[degs_to_include, ]

  # Calculate enrichment scores
  scores_tmp = aggregate(t(profiles_cc_genes), by = list(lig_seurat@meta.data$sample), sum)
  rownames(scores_tmp) = scores_tmp$Group.1
  scores_tmp$Group.1 = NULL
  scores = apply(scores_tmp, 1, sum)
  scores = scores - scores["PBS"]
  scores_df = data.frame(scores)
  scores_df = scores_df[cytokines,, drop = FALSE]
  scores_pvals = sapply(cytokines, function(x){
    test_res = wilcox.test(apply(profiles_cc_genes[, rownames(subset(lig_seurat@meta.data, sample == x))], 2, sum),
                           apply(profiles_cc_genes[, rownames(subset(lig_seurat@meta.data, sample == "PBS"))], 2, sum));
    return(test_res$p.value)
  })
  scores_df$pval = as.vector(scores_pvals)
  scores_df$padj = p.adjust(scores_df$pval, method = "fdr")

  # Add pseudocount for log transform
  scores_df$nlog10_padj = -log10(scores_df$padj+0.000001)

  # Sort results by p-value
  scores_df$cytokine = rownames(scores_df)
  scores_df$cytokine = factor(scores_df$cytokine, levels = scores_df$cytokine[order(scores_df$nlog10_padj)])

  ggplot(scores_df, aes(x = cytokine, y = nlog10_padj)) +
    geom_hline(yintercept = 2, color = "blue", linetype = 'dotted') +
    geom_bar(stat = "identity", fill = "orange") +
    coord_flip() +
    theme_classic() +
    xlab("Cytokine response") +
    ylab("IREA-GeneSet -log10 (FDR)")

}



#### Immune response enrichment score functions #####

#' Get enrichment score - Fisher's test method
#'
#' \code{GeneSetEnrichmentHyperTest} Compute enrichment score using the pathway overrepresentation test
#' which is commonly involves the hypergeometric test
#'
#' @param    degs              A list of differentially expresssed genes that user would like to investigate
#' @param    celltype          Choose from one of the listed cell types that most resemble the input
#'

GeneSetEnrichmentHyperTest = function(degs,
                                      celltype) {
  `%notin%` = Negate(`%in%`)
  celltypes = c("B_cell", "cDC1", "cDC2", "eTAC", "ILC", "Macrophage", "MigDC",
                "Monocyte", "Neutrophil", "NK_cell", "pDC", "T_cell_CD4", "T_cell_CD8",
                "T_cell_gd", "Treg")
  #if (celltype %notin% celltypes) stop("cell type must be one of the following: ")

  #ref_deg = readRDS("../data/ref_deg.Rda")
  #ref_deg_sig = subset(ref_deg, p_val < 0.01 & avg_logFC > 0)
  #saveRDS(ref_deg_sig, file = "../data/ref_deg_sig.Rda")

  ref_deg_sig = readRDS("../data/ref_deg_sig.Rda")
  ref_expressed_genes = readRDS("../data/ref_expressed_genes_per_celltype.Rda")

  # subset into the celltype of interest
  ref_deg_sig_celltype = subset(ref_deg_sig, celltype == celltype)
  ref_expressed_genes_celltype = subset(ref_expressed_genes, celltype == celltype)

  # examples
  # degs = c("Isg15", "Irf7", "Il1b", "Isg20")

  samples = sort(unique(ref_deg_sig$cytokine))

  test_res = c()
  for (ss in samples) {
    markers_ss = subset(ref_deg_sig_celltype, cytokine == ss)

    num_overlap = length(intersect(degs, markers_ss$gene))
    num_only_user = length(degs) - num_overlap
    num_only_ref = nrow(markers_ss) - num_overlap
    num_neither = length(setdiff(ref_expressed_genes_celltype$genes_expressed,
                                 c(markers_ss$gene, degs)))

    mat_test = matrix(c(num_overlap, num_only_user, num_only_ref, num_neither), nrow = 2)

    res = fisher.test(mat_test)$p.value
    test_res = c(test_res, res)
  }

  # construct a result matrix and order by p-value
  df_fisher = data.frame(sample = samples, test_pval = test_res)
  df_fisher = df_fisher[order(df_fisher$test_pval), ]

  return(df_fisher)
}




#### Immune response enrichment score functions #####

#' IREA analysis for transcriptome matrix input
#'
#' \code{GetEnrichmentScoreProjection} Compute enrichment score using the pathway overrepresentation test
#' which is commonly involves the hypergeometric test
#'
#' @param    input_profile    Gene expression matrix.
#' @param    input_celltype    Choose from one of the listed cell types that most resemble the input
#'


GetEnrichmentScoreProjection = function(input_profile,
                                        input_celltype) {
  `%notin%` = Negate(`%in%`)
  celltypes = c("B_cell", "cDC1", "cDC2", "eTAC", "ILC", "Macrophage", "MigDC",
                "Monocyte", "Neutrophil", "NK_cell", "pDC", "T_cell_CD4", "T_cell_CD8",
                "T_cell_gd", "Treg")
  if (input_celltype %notin% celltypes) stop(paste0("cell type must be one of the following: ", paste0(celltypes, collapse = ", ")))


  # Load reference data
  # TODO: make better reference data file after re-identifying cell types
  load("../../../Hacohen/ligands/rdata/200417-ligands-alldata-seurat-p3-subsample.RData")
  lig_seurat_sub = subset(lig_seurat, celltype == input_celltype)
  profiles_cc = as.matrix(lig_seurat_sub@assays[['RNA']][,])

  # Choose intersection genes
  genes_common = intersect(rownames(input_profile), rownames(profiles_cc))

  # Make the final input and final reference by selecting the common genes
  profiles_cc = profiles_cc[genes_common, ]
  input_profile = input_profile[genes_common, , drop = FALSE]

  # Select the genes that are much different between cytokine treatment and PBS
  input_profile_agg = aggregate(t(profiles_cc), by = list(lig_seurat_sub@meta.data$sample), mean)
  rownames(input_profile_agg) = input_profile_agg$Group.1
  input_profile_agg$Group.1 = NULL

  gene_diff = apply(input_profile_agg, 2, function(x){max(x-x['PBS'])})
  genes_large_diff = names(gene_diff)[gene_diff>0.25]

  profiles_cc = profiles_cc[genes_large_diff, ]
  input_profile = input_profile[genes_large_diff, , drop = FALSE]

  dist_input_mat = cbind(input_profile, profiles_cc)

  # Project all input vectors onto reference panel
  library(philentropy)
  library(reshape2)
  projection_scores_raw = distance(t(dist_input_mat), method = "cosine")

  # Select the part of the distance matrix that correspond to the results
  projection_scores = projection_scores_raw[1:ncol(input_profile), (ncol(input_profile)+1):ncol(dist_input_mat)]
  rownames(projection_scores) = colnames(input_profile)
  colnames(projection_scores) = colnames(profiles_cc)

  input_samples = colnames(input_profile)
  reference_samples = setdiff(sort(unique(lig_seurat_sub@meta.data$sample)), "PBS")

  mat_pval = matrix(NA, length(input_samples), length(reference_samples),
                    dimnames = list(input_samples, reference_samples))
  df_irea = melt(mat_pval)
  df_irea$ES = NA
  names(df_irea) = c("Sample", "Cytokine", "pval", "ES")

  for (ii in input_samples) {
    for (x in reference_samples) {
      cols_cytokine = which(lig_seurat_sub@meta.data$sample == x);
      cols_pbs = which(lig_seurat_sub@meta.data$sample == "PBS");

      if (length(cols_cytokine) > 10) {
        test_res = wilcox.test(projection_scores[ii, cols_cytokine], projection_scores[ii, cols_pbs])
        pval = test_res$p.value
        meandiff = mean(projection_scores[ii, cols_cytokine]) - mean(projection_scores[ii, cols_pbs])
      } else {
        pval = NULL;
        meandiff = NULL}

      df_irea[df_irea$Sample == ii & df_irea$Cytokine == x, "pval"] = pval
      df_irea[df_irea$Sample == ii & df_irea$Cytokine == x, "ES"] = meandiff

    }
  }

  df_irea$padj = p.adjust(df_irea$pval, method = "fdr")

  df_irea_sig = subset(df_irea, padj < 0.01)

  return(df_irea)
}





