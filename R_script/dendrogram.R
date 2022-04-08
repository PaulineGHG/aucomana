###############################################################################
# install packages
###############################################################################

install.packages("pvclust")
install.packages("dendextend")
install.packages("dplyr")
install.packages("ape")
install.packages("ade4")

###############################################################################
# load packages
###############################################################################

library(pvclust)
library(dendextend)
library(dplyr)
library(ape)
library(ade4)

###############################################################################
# load data
###############################################################################

setwd("~/Documents/Analysis_runs")
name = "Run 04"
reactions_dataframe = read.table("data/runs/run04/analysis/all/reactions.tsv", sep = "\t", header = T, row.names = "reaction")
reactions_dataframe = select(reactions_dataframe, -colnames(select(reactions_dataframe, ends_with("..sep...") | ends_with("_formula"))))

###############################################################################
# vector of algae categories
###############################################################################

tsv_org_file = read.table("data/species_group.tsv", sep = "\t", header = F)

brown = tsv_org_file[,1][tsv_org_file[,2] == "brown"]
diatoms = tsv_org_file[,1][tsv_org_file[,2] == "diatoms"]
nano = tsv_org_file[,1][tsv_org_file[,2] == "nano"]
long_read = tsv_org_file[,1][tsv_org_file[,3] == "LR"]

###############################################################################
# select pvalue
###############################################################################

# pval = result$edges$au
# pval
# signif = which(pval > 0.95)
# signif

###############################################################################
# functions
###############################################################################

get_lab_colors = function(dendrogram){
  lab = labels(dendrogram)
  col_lab = c()
  for (i in 1:length(lab)){
    if (lab[i] %in% brown){
      col_lab[i] = "#a3724b"
    }
    else if(lab[i] %in% diatoms){
      col_lab[i] = "#c9ad4b"
    }
    else if(lab[i] %in% nano){
      col_lab[i] = "#e7baaa"
    }
    else{
      col_lab[i] = "black"
    }
  }
  return(col_lab)
}

get_lr_leaves = function(dendrogram){
  lab = labels(dendrogram)
  col_pch = rep(0, length(lab))
  pch = rep(NA, length(lab))
  for (i in 1:length(lab)){
    if (lab[i] %in% long_read){
      col_pch[i] = "blue"
      pch[i] = 19
    }
  }
  return(list(pch, col_pch))
}

get_nodes = function(dendrogram){
  n_nodes = nnodes(dendrogram)
  leaf = get_nodes_attr(dendrogram, "leaf")
  nodes_col = rep(1, n_nodes)
  nodes_pch = rep(NA, n_nodes)
  iter = c()
  for (i in 2:n_nodes){
    if (is.na(leaf[i])){
      iter = c(iter, i)
    }
  }
  print(length(iter))
  for (i in 1:length(iter)){
    if (i %in% signif){
      nodes_col[iter[i]] = 2
      nodes_pch[iter[i]] = 19
    }
  }
  return(list(nodes_col, nodes_pch))
}

get_line_colors = function(phylo, dendro){
  lab_dend = labels(dendro)
  lab = labels(phylo)
  col_lab = c()
  for (i in 1:length(lab)){
    if (lab[i] %in% brown & lab[i] %in% lab_dend){
      col_lab = c(col_lab, "#a3724b")
    }
    else if (lab[i] %in% diatoms & lab[i] %in% lab_dend){
      col_lab = c(col_lab, "#c9ad4b")
    }
    else if (lab[i] %in% nano & lab[i] %in% lab_dend){
      col_lab = c(col_lab, "#e7baaa")
    }
    else if (lab[i] %in% lab_dend){
      col_lab = c(col_lab, "black")
    }
  }
  return(col_lab)
}

###############################################################################
# dendrogram
###############################################################################

result = pvclust(reactions_dataframe, method.dist="binary", method.hclust="complete", nboot=10000, parallel=TRUE)
dend = as.dendrogram(result)

d_col_lab = get_lab_colors(dend)
# d_nodes_list = get_nodes(dend)
# d_nodes_col = d_nodes_list[[1]]
# d_nodes_pch = d_nodes_list[[2]]
d_leaf = get_lr_leaves(dend) ; d_leaf_pch = d_leaf[[1]] ; d_leaf_pchcol = d_leaf[[2]] 

dend = dend %>% 
  set("branches_lwd", 2) %>% 
  set("by_labels_branches_col", value = brown, TF_values = c("#a3724b",Inf)) %>%
  set("by_labels_branches_col", value = diatoms, TF_values = c("#c9ad4b", Inf)) %>%
  set("by_labels_branches_col", value = nano, TF_values = c("#e7baaa",Inf)) %>%
  set("labels_col", d_col_lab) %>%
  set("leaves_pch", d_leaf_pch) %>%
  set("leaves_col", d_leaf_pchcol) %>%
  # set("nodes_pch", d_nodes_pch) %>%
  # set("nodes_col", d_nodes_col) %>%
  sort(type = "nodes") 


par(mar=c(3,3,1,15))
plot_horiz.dendrogram(dend, side = F, edge.root = T, main = paste(name, " Metabolic Dendrogram"))
plot(result)

###############################################################################
# phylo tree
###############################################################################

phylo = ape::read.nexus("data/Phaeoexplorer_MLtree_rooted.nex")
# phylo = ape::read.nexus("data/Ref_tree_run01.nex.txt")
phylo = chronos(phylo)
phylo = as.dendrogram(phylo)

p_col_lab = get_lab_colors(phylo)
p_leaf = get_lr_leaves(phylo)
p_leaf_pch = p_leaf[[1]] ; p_leaf_pchcol = p_leaf[[2]] 

phylo = phylo %>% 
  set("branches_lwd", 2) %>%
  set("by_labels_branches_col", value = brown, TF_values = c("#a3724b",Inf)) %>%
  set("by_labels_branches_col", value = diatoms, TF_values = c("#c9ad4b", Inf)) %>%
  set("by_labels_branches_col", value = nano, TF_values = c("#e7baaa",Inf)) %>%
  set("labels_col", p_col_lab) %>%
  set("leaves_pch", p_leaf_pch) %>%
  set("leaves_col", p_leaf_pchcol) %>%
  sort(type = "nodes") 

###############################################################################
# tanglegram
###############################################################################

col_lines = get_line_colors(phylo, dend)
plot_horiz.dendrogram(phylo, side = F)
d1 = intersect_trees(phylo, dend)
tanglegram(d1, margin_inner = 13, color_lines = col_lines, lwd = 2, 
           main_left = "Original phylogeny", main_right = paste(name, "\nMetabolic Dendrogram"),
           margin_outer = 2, margin_top = 5)

###############################################################################
# corelation
###############################################################################

# From cophenetic: The cophenetic distance between two observations that have 
# been clustered is defined to be the intergroup dissimilarity at which the two 
# observations are first combined into a single cluster. Note that this distance 
# has many ties and restrictions.
# 
# cor_cophenetic calculates the correlation between two cophenetic distance 
# matrices of the two trees.
# 
# The value can range between -1 to 1. With near 0 values meaning that the two 
# trees are not statistically similar. For exact p-value one should result to a 
# permutation test. One such option will be to permute over the labels of one 
# tree many times, and calculating the distriubtion under the null hypothesis 
# (keeping the trees topologies constant).
# 
# Notice that this measure IS affected by the height of a branch.

cor_cophenetic(d1, method_coef = "pearson")

# Baker's Gamma (see reference) is a measure of accosiation (similarity) between
# two trees of heirarchical clustering (dendrograms).
# 
# It is calculated by taking two items, and see what is the heighst possible 
# level of k (number of cluster groups created when cutting the tree) for which 
# the two item still belongs to the same tree. That k is returned, and the same 
# is done for these two items for the second tree. There are n over 2 
# combinations of such pairs of items from the items in the tree, and all of 
# these numbers are calculated for each of the two trees. Then, these two sets 
# of numbers (a set for the items in each tree) are paired according to the 
# pairs of items compared, and a spearman correlation is calculated.
# 
# The value can range between -1 to 1. With near 0 values meaning that the two
# trees are not statistically similar. For exact p-value one should result to a
# permutation test. One such option will be to permute over the labels of one 
# tree many times, and calculating the distriubtion under the null hypothesis 
# (keeping the trees topologies constant).
# 
# Notice that this measure is not affected by the height of a branch but only of
# its relative position compared with other branches.

cor_bakers_gamma(d1, to_plot = T)

# Fowlkes-Mallows index (see references) is an external evaluation method that 
# is used to determine the similarity between two clusterings (clusters obtained 
# after a clustering algorithm). This measure of similarity could be either 
# between two hierarchical clusterings or a clustering and a benchmark 
# classification. A higher the value for the Fowlkes-Mallows index indicates a 
# greater similarity between the clusters and the benchmark classifications.
# 
# The default Bk plot comes with a line with dots (type "b") of the Bk values. 
# Also with a fragmented (lty=2) line (of the same color) of the expected Bk 
# line under H0, And a solid red line of the upper critical Bk values for 
# rejection.

Bk_plot(d1[[1]], d1[[2]], rejection_line_permutation = T)
