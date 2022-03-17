# install packages

install.packages("pvclust")
install.packages("dendextend")
install.packages("dplyr")
install.packages("ape")

# load packages

library(pvclust)
library(dendextend)
library(dplyr)
library(ape)


setwd("~/Documents/Reactions_loss")

# load data

name = "run03"
reactions_dataframe = read.table("data/reactions_data/run03_reactions.tsv", sep = "\t", header = T, row.names = "reaction")
reactions_dataframe = select(reactions_dataframe, -colnames(select(reactions_dataframe, ends_with("..sep...") | ends_with("_formula"))))

result = pvclust(reactions_dataframe, method.dist="binary", method.hclust="complete", nboot=10000, parallel=TRUE)

# vector of algae categories

brown = c("Laminarionema_elsbetiae", "Undaria_pinnatifida_Kr2015", "Cladosiphon_okamuranus", "Nemacystus_decipiens", "Saccharina_japonica", 
          "Pleurocladia_lacustris", "Ectocarpus_siliculosus", "Ectocarpus_subulatus", "Saccharina_latissima_FEMALE", "Dictyota_dichotoma_m",
          "Fucus_serratus_MALE", "Desmarestia_herbacea_m", "Ectocarpus_siliculosus_m", "Porterinema_fluviatile", "Ectocarpus_fasciculatus_m",
          "Chordaria_linearis", "Scytosiphon_promiscuus_MALE", "Ectocarpus_crouaniorum_m", "Ectocarpus_species7")

diatoms = c("Thalassiosira_pseudonana", "Fragilariopsis_cylindrus", "Fistulifera_solaris", "Phaeodactylum_tricornutum")

nano = c("Nannochloropsis_gaditana")

long_read = c("Pleurocladia_lacustris", "Saccharina_latissima_FEMALE", "Dictyota_dichotoma_m", "Fucus_serratus_MALE", "Desmarestia_herbacea_m",
              "Ectocarpus_siliculosus_m", "Porterinema_fluviatile", "Ectocarpus_fasciculatus_m", "Chordaria_linearis", "Scytosiphon_promiscuus_MALE",
              "Ectocarpus_crouaniorum_m", "Schizocladia_ischiensis")

# pval = result$edges$au
# pval
# signif = which(pval > 0.95)
# signif

# functions

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

# dendrogram

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
plot_horiz.dendrogram(dend, side = F, edge.root = T, main = paste(name, " Predicted phylogeny"))
plot(result)

# phylo tree

phylo = ape::read.nexus("data/Ref_tree_run01.nex.txt")
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

col_lines = get_lab_colors(phylo)
plot_horiz.dendrogram(phylo, side = F)
d1 = dendlist(phylo, dend)
tanglegram(d1, margin_inner = 13, color_lines = col_lines, lwd = 2, 
           main_left = "Original phylogeny", main_right = paste(name, "\nPredicted phylogeny"),
           margin_outer = 2, margin_top = 5)

