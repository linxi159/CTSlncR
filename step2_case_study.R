##########################################################################################
########## Running scripts for exploring cell-type-specific lncRNA regulation ############
##########################################################################################

## Load required R packages, please firstly install the following R packages before running scripts
library(pracma)
library(WGCNA)
library(igraph)
library(miRspongeR)
library(biclique)
library(corrplot)
library(dendextend)
library(ggplot2)
library(ggsignif)
library(cowplot)
library(clusterProfiler)
library(msigdbr)
library(vroom)
library(doParallel)

## Load utility functions
#setwd("~/r_workplace/9.gene_regulation_lncRNA-mRNA/methods/CTSlncR-master")
setwd("~/CTSlncR")
source("CTSlncR.R")

## Load prepared datasets
ASD <- read.csv("./data/ASD_gene_lncRNAs_mRNAs.csv")
ASD <- as.matrix(ASD)

## Load the preprocessing single-cell sequencing data including filtering genes with expression <= 10 cells, and log e(x+1).
lncRNA_scRNA_raw  <- read.csv("./data/GSE71315/GSE71315_scell_ncounts_genes_thresh_lncRNA.csv")
rownames(lncRNA_scRNA_raw ) <- lncRNA_scRNA_raw[,1]
lncRNA_scRNA_raw  <- lncRNA_scRNA_raw[,-1]
lncRNA_scRNA_raw   <- t(lncRNA_scRNA_raw)
lncRNA_scRNA_raw   <- as.data.frame(lncRNA_scRNA_raw)

mRNA_scRNA_raw <- read.csv("./data/GSE71315/GSE71315_scell_ncounts_genes_thresh_mRNA.csv")
rownames(mRNA_scRNA_raw) <- mRNA_scRNA_raw[,1]
mRNA_scRNA_raw <- mRNA_scRNA_raw[,-1]
mRNA_scRNA_raw  <- t(mRNA_scRNA_raw)
mRNA_scRNA_raw <- as.data.frame(mRNA_scRNA_raw)

lncRNA_gene_pre <- read.csv("./data/lncRNA_gene-pre.csv")
lncRNA_mRNA_exp <- read.csv("./data/lncRNA_mRNA-exp.csv")
lncRNA_mRNA_exp <- as.matrix(lncRNA_mRNA_exp)

## compute the average expression values of duplicate genes
## and remove genes with constant expression values in all cells

    # Transformation using log2(x+1)
    lncRNA_scRNA_norm <- log2(lncRNA_scRNA_raw+1)
    mRNA_scRNA_norm <- log2(mRNA_scRNA_raw+1) 

    # Compute the average expression values of duplicate genes
    lncRNA_scRNA_norm_average <- Averg_Duplicate(lncRNA_scRNA_norm)
    mRNA_scRNA_norm_average <- Averg_Duplicate(mRNA_scRNA_norm)

    # Remove genes with constant expression values in all cells
    lncRNA_scRNA_norm_sd <- unlist(lapply(seq(dim(lncRNA_scRNA_norm_average)[2]), function(i) sd(lncRNA_scRNA_norm_average[, i])))
    lncRNA_scRNA_norm_filter <- lncRNA_scRNA_norm_average[, which(lncRNA_scRNA_norm_sd > 0)]
    mRNA_scRNA_norm_sd <- unlist(lapply(seq(dim(mRNA_scRNA_norm_average)[2]), function(i) sd(mRNA_scRNA_norm_average[, i])))
    mRNA_scRNA_norm_filter <- mRNA_scRNA_norm_average[, which(mRNA_scRNA_norm_sd > 0)]

# save    
# save(lncRNA_scRNA_norm, mRNA_scRNA_norm, lncRNA_scRNA_norm_filter,mRNA_scRNA_norm_filter, lncRNA_gene_pre, lncRNA_mRNA_exp,file = "Exp_247_lncRNAs_10208_mRNAs_276_single_cells_GSE71315.RData")
load("Exp_247_lncRNAs_10208_mRNAs_276_single_cells_GSE71315.RData")

# # 多核并行计算
# library(parallel)
# #CSlncR_network_bootstrap_null <- CSlncR_net_bootstrap(lncRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, boxsize = 0.1, bootstrap_betw_point = 5, bootstrap_num = 1, p.value.cutoff = 0.01)
# detectCores(logical = F)  # 8
# mc <- getOption("mc.cores", 8)
# system.time({ CSlncR_network_bootstrap_null <- mclapply(list(1), CSlncR_net_bootstrap,lncR=lncRNA_scRNA_norm_filter,mR=mRNA_scRNA_norm_filter, mc.cores = mc);});
# stopCluster(mc);
# 
# #load("CSlncR_network_bootstrap_null_GSE71315.RData")
# CSlncR_network_bootstrap_null = CSlncR_network_bootstrap_null[[1]]

## Discovering cell-specific lncRNA-mRNA regulatory network   
    CSlncR_network_bootstrap_null <- CSlncR_net_bootstrap(lncRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, 
                                                          boxsize = 0.1, bootstrap_betw_point = 5, 
                                                          bootstrap_num = 1, p.value.cutoff = 0.01)
    prior_graph <- make_graph(c(t(lncRNA_gene_pre[, 1:2])), directed = TRUE)
    CSlncR_network_bootstrap_null_graph <- lapply(seq(CSlncR_network_bootstrap_null), function(i) make_graph(c(t(CSlncR_network_bootstrap_null[[i]][, 1:2])), directed = TRUE))
    CSlncR_network_bootstrap <- lapply(seq(CSlncR_network_bootstrap_null), function(i) as_data_frame(CSlncR_network_bootstrap_null_graph[[i]] %s% prior_graph))
    # without prior knowledge
    CSlncR_network_bootstrap_ <- lapply(seq(CSlncR_network_bootstrap_null), function(i) as_data_frame(CSlncR_network_bootstrap_null_graph[[i]]))

# save    
save(CSlncR_network_bootstrap_null_graph, CSlncR_network_bootstrap,file = "CSlncR_network_bootstrap_null_graph_CSlncR_network_bootstrap_GSE71315.RData")
save(CSlncR_network_bootstrap_null_graph,file = "CSlncR_network_bootstrap_null_graph_GSE71315.RData")
save(CSlncR_network_bootstrap,file = "CSlncR_network_bootstrap_GSE71315.RData")
save(CSlncR_network_bootstrap_,file = "CSlncR_network_bootstrap_GSE71315_.RData")
load("CSlncR_network_bootstrap_null_graph_CSlncR_network_bootstrap_GSE71315.RData")       
load("CSlncR_network_bootstrap_GSE71315.RData")  # ***
load("CSlncR_network_bootstrap_GSE71315_.RData") 

# divide different developmental stages for GSE71315
CSlncR_network_bootstrap_src <- CSlncR_network_bootstrap
# GW16  GW21     A1      A2       S44      S46
# GW16  GW21     GW20.5  GW20.5   GW19.5   GW23.5
# 1:26  27:50    51:115  116:173  174:199  200:276
#  26    24            123           26       77
# without prior knowledge    
CSlncR_network_bootstrap_src_  <- CSlncR_network_bootstrap_


#############################################################################################################################
# [1]GW16 to GW23.5
# Merging single-cell networks from 5 different developmental stages, separately, [CSlncR_network_bootstrap]
CSlncR_network_bootstrap_GW16 <- CSlncR_network_bootstrap[c(1:26)]
CSlncR_network_bootstrap_GW19.5 <- CSlncR_network_bootstrap[c(174:199)]
CSlncR_network_bootstrap_GW20.5 <- CSlncR_network_bootstrap[c(51:173)]
CSlncR_network_bootstrap_GW21 <- CSlncR_network_bootstrap[c(27:50)]
CSlncR_network_bootstrap_GW23.5 <- CSlncR_network_bootstrap[c(200:276)]

# After merging, from GW16 to GW23.5
tmp_1 <- append(CSlncR_network_bootstrap_GW16,CSlncR_network_bootstrap_GW19.5)
tmp_1 <- append(tmp_1,CSlncR_network_bootstrap_GW20.5)
tmp_1 <- append(tmp_1,CSlncR_network_bootstrap_GW21)
CSlncR_network_bootstrap_GW16_to_GW23.5 <-  append(tmp_1,CSlncR_network_bootstrap_GW23.5)
                                                
CSlncR_merge_GW <- function(tmp) {
  df_GW <- tmp[[1]]
  for (i in seq(length(tmp))) {
    if(i != 1 ){
      df_GW <- rbind(df_GW,tmp[[i]])
    }
  }
  return(df_GW)
}

df_GW16 <- CSlncR_merge_GW(CSlncR_network_bootstrap_GW16)
df_GW19.5 <- CSlncR_merge_GW(CSlncR_network_bootstrap_GW19.5)
df_GW20.5 <- CSlncR_merge_GW(CSlncR_network_bootstrap_GW20.5)
df_GW21 <- CSlncR_merge_GW(CSlncR_network_bootstrap_GW21)
df_GW23.5 <- CSlncR_merge_GW(CSlncR_network_bootstrap_GW23.5)
CSlncR_network_bootstrap_GW16to23.5 <- list(df_GW16, df_GW19.5, df_GW20.5, df_GW21, df_GW23.5)
CSlncR_network_bootstrap <- CSlncR_network_bootstrap_GW16to23.5

cell_number <- 5 # 5 development stages 
#without prior knowledge
# CSlncR_network_bootstrap_GW16_ <- CSlncR_network_bootstrap_[c(1:26)]
# CSlncR_network_bootstrap_GW19.5_ <- CSlncR_network_bootstrap_[c(174:199)]
# CSlncR_network_bootstrap_GW20.5_ <- CSlncR_network_bootstrap_[c(51:173)]
# CSlncR_network_bootstrap_GW21_ <- CSlncR_network_bootstrap_[c(27:50)]
# CSlncR_network_bootstrap_GW23.5_ <- CSlncR_network_bootstrap_[c(200:276)]
#
# df_GW16_ <- CSlncR_merge_GW(CSlncR_network_bootstrap_GW16_)
# df_GW19.5_ <- CSlncR_merge_GW(CSlncR_network_bootstrap_GW19.5_)
# df_GW20.5_ <- CSlncR_merge_GW(CSlncR_network_bootstrap_GW20.5_)
# df_GW21_ <- CSlncR_merge_GW(CSlncR_network_bootstrap_GW21_)
# df_GW23.5_ <- CSlncR_merge_GW(CSlncR_network_bootstrap_GW23.5_)
# CSlncR_network_bootstrap_GW16to23.5_ <- list(df_GW16_, df_GW19.5_, df_GW20.5_, df_GW21_, df_GW23.5_)
# CSlncR_network_bootstrap_ <- CSlncR_network_bootstrap_GW16to23.5_
# save(CSlncR_network_bootstrap_,file = "CSlncR_network_bootstrap_GW20.5toGW23.5_.RData")
#load("CSlncR_network_bootstrap_GW20.5toGW23.5_.RData")
#rm(list= "CSlncR_network_bootstrap_src_","CSlncR_network_bootstrap_GW20.5_","CSlncR_network_bootstrap_GW23.5_","CSlncR_network_bootstrap_GW16to23.5_","df_GW16_","df_GW19.5_","df_GW20.5_","df_GW21_","df_GW23.5_")

# Integrate expression values in the order GW16 to GW23.5
lncRNA_scRNA_norm_filter <- rbind(lncRNA_scRNA_norm_filter[1:26,],lncRNA_scRNA_norm_filter[174:199,],lncRNA_scRNA_norm_filter[51:173,],lncRNA_scRNA_norm_filter[27:50,],lncRNA_scRNA_norm_filter[200:276,])
mRNA_scRNA_norm_filter <- rbind(mRNA_scRNA_norm_filter[1:26,],mRNA_scRNA_norm_filter[174:199,],mRNA_scRNA_norm_filter[51:173,],mRNA_scRNA_norm_filter[27:50,],mRNA_scRNA_norm_filter[200:276,])
lncRNA_mRNA_scRNA_norm_filter <- cbind(lncRNA_scRNA_norm_filter,mRNA_scRNA_norm_filter)
##############################################################################################################################


##############################################################################################################################
# CSlncR_network_bootstrap_random
# random method # 无放回抽取
lncR_names <- colnames(lncRNA_scRNA_norm_filter)[sample(1:247, 247/8, replace=F)] # 247
mR_names <- colnames(mRNA_scRNA_norm_filter)[sample(1:10208, 10208/8, replace=F)] # 10208

# CSlncR_network_bootstrap_random_ <- data.frame(from="0",to="0")
# cnt <- 0
# for (i in seq(length(lncR_names))) {
#   for (j in seq(length(mR_names))) {
#     cnt  <- cnt + 1
#     CSlncR_network_bootstrap_random_[cnt,1] <- lncR_names[i]
#     CSlncR_network_bootstrap_random_[cnt,2] <- mR_names[j]
#   }
# }
# 
# save(CSlncR_network_bootstrap_random_,file = "CSlncR_network_bootstrap_random_.RData")
load("CSlncR_network_bootstrap_random_.RData")
# CSlncR_network_bootstrap_random <-list()
# for (i in 1:cell_number){
#   CSlncR_network_bootstrap_random[[i]] <- CSlncR_network_bootstrap_random_
# }
CSlncR_merge_5 <- function(tmp,cell_number) {
  CSlncR_network_bootstrap_random <-list()
  for (i in 1:cell_number){
    if(i==1){
      cells_num <- 26
    }
    if(i==2){
      cells_num <- 26
    } 
    if(i==3){
      cells_num <- 123
    } 
    if(i==4){
      cells_num <- 24
    } 
    if(i==5){
      cells_num <- 77
    }
  
    df_random <- tmp
    for (j in seq(cells_num)) {
      if(j != 1 ){
        df_random <- rbind(df_random,tmp)
      }
    }
  
    CSlncR_network_bootstrap_random[[i]] <- df_random
  }

  return(CSlncR_network_bootstrap_random)
}

CSlncR_network_bootstrap_random  <- CSlncR_merge_5(CSlncR_network_bootstrap_random_,cell_number)
# save(CSlncR_network_bootstrap_random,file = "CSlncR_network_bootstrap_random.RData")
# load("CSlncR_network_bootstrap_random.RData")

# CSlncR_network_bootstrap_LncRNA2Target
CSlncR_network_bootstrap_LncRNA2Target_ <- read.csv("./data/LncRNA2Target.csv")
# CSlncR_network_bootstrap_LncRNA2Target <-list()
# for (i in 1:cell_number){
#   CSlncR_network_bootstrap_LncRNA2Target[[i]] <- CSlncR_network_bootstrap_LncRNA2Target_
# }
CSlncR_network_bootstrap_LncRNA2Target  <- CSlncR_merge_5(CSlncR_network_bootstrap_LncRNA2Target_,cell_number)


###########################################################################################################################
################################### 11 downstream computation - 3 result analysis ########################################
###########################################################################################################################
# downstream [1] validated
## Experimentally validated cell-specific lncRNA-mRNA interactions, the ground-truth (lncRTarget variable) is from the literature of NPInter database.
    lncRTarget_graph <- make_graph(c(t(lncRNA_mRNA_exp[, 1:2])), directed = TRUE)
    CSlncR_network_bootstrap_graph <- lapply(seq(CSlncR_network_bootstrap), function(i) make_graph(c(t(CSlncR_network_bootstrap[[i]][, 1:2])), directed = TRUE))
    CSlncR_network_bootstrap_validated <- lapply(seq(CSlncR_network_bootstrap), function(i) as_data_frame(CSlncR_network_bootstrap_graph[[i]] %s% lncRTarget_graph))

    # GW16_to_GW23.5
    CSlncR_network_bootstrap_GW16_to_GW23.5_graph <- lapply(seq(CSlncR_network_bootstrap_GW16_to_GW23.5), function(i) make_graph(c(t(CSlncR_network_bootstrap_GW16_to_GW23.5[[i]][, 1:2])), directed = TRUE))
    
    # without prior knowledge 
    #CSlncR_network_bootstrap_graph_ <- lapply(seq(CSlncR_network_bootstrap_), function(i) make_graph(c(t(CSlncR_network_bootstrap_[[i]][, 1:2])), directed = TRUE))
    #CSlncR_network_bootstrap_validated_ <- lapply(seq(CSlncR_network_bootstrap_), function(i) as_data_frame(CSlncR_network_bootstrap_graph_[[i]] %s% lncRTarget_graph))
    #save(CSlncR_network_bootstrap_validated_,file = "CSlncR_network_bootstrap_validated_.RData")
    load("CSlncR_network_bootstrap_validated_.RData")  
    
    # random
    CSlncR_network_bootstrap_graph_random <- lapply(seq(CSlncR_network_bootstrap_random), function(i) make_graph(c(t(CSlncR_network_bootstrap_random[[i]][, 1:2])), directed = TRUE))
    CSlncR_network_bootstrap_validated_random <- lapply(seq(CSlncR_network_bootstrap_random), function(i) as_data_frame(CSlncR_network_bootstrap_graph_random[[i]] %s% lncRTarget_graph))
    # LncRNA2Target 
    CSlncR_network_bootstrap_graph_LncRNA2Target <- lapply(seq(CSlncR_network_bootstrap_LncRNA2Target), function(i) make_graph(c(t(CSlncR_network_bootstrap_LncRNA2Target[[i]][, 1:2])), directed = TRUE))
    CSlncR_network_bootstrap_validated_LncRNA2Target <- lapply(seq(CSlncR_network_bootstrap_LncRNA2Target), function(i) as_data_frame(CSlncR_network_bootstrap_graph_LncRNA2Target[[i]] %s% lncRTarget_graph))
    
# downstream [2] ASD
## ASD-related cell-specific lncRNA-mRNA interactions. the list of ASD-related lncRNAs and mRNAs (ASD variable) is from SFARI tools, respectively    
    CSlncR_network_bootstrap_ASD <- lapply(seq(CSlncR_network_bootstrap), function(i) CSlncR_network_bootstrap[[i]][intersect(which(CSlncR_network_bootstrap[[i]][, 1] %in% as.matrix(ASD)), which(CSlncR_network_bootstrap[[i]][, 2] %in% as.matrix(ASD))), ])
    
# downstream [3] interactions  ####[1] Result analysis: The lncRNA regulation in each cell type is unique####
    ## Overlap of cell-type-specific lncRNA-mRNA interactions across cell types  
    Overlap_network_bootstrap <- Overlap_net_interaction(CSlncR_network_bootstrap, Intersect_num = round(length(CSlncR_network_bootstrap)*1))# default:1
    Overlap_network_bootstrap_union <- Overlap_net_interaction(CSlncR_network_bootstrap, Intersect_num = 1)
    Overlap_network_bootstrap_two <- Overlap_net_interaction(CSlncR_network_bootstrap, Intersect_num = 2) 
    
    Overlap_network_bootstrap_union_graph <- make_graph(c(t(Overlap_network_bootstrap_union[, 1:2])), directed = TRUE)
    Overlap_network_bootstrap_two_graph <- make_graph(c(t(Overlap_network_bootstrap_two[, 1:2])), directed = TRUE)
    Overlap_network_bootstrap_rewired <- as_data_frame(Overlap_network_bootstrap_union_graph %m% Overlap_network_bootstrap_two_graph)

    Overlap_network_bootstrap_graph <- make_graph(c(t(Overlap_network_bootstrap[, 1:2])), directed = TRUE)
    Overlap_network_bootstrap_conserved_validated <- as_data_frame(Overlap_network_bootstrap_graph %s% lncRTarget_graph)
    Overlap_network_bootstrap_rewired_graph <- make_graph(c(t(Overlap_network_bootstrap_rewired[, 1:2])), directed = TRUE)
    Overlap_network_bootstrap_rewired_validated <- as_data_frame(Overlap_network_bootstrap_rewired_graph %s% lncRTarget_graph)

    #write.csv(Overlap_network_bootstrap_rewired,"~/data/tmp.csv")

# downstream [4] ASD-related regulatory network 
    ## ASD-related conserved and rewired lncRNA-mRNA regulatory network    
    Overlap_network_bootstrap_ASD <- Overlap_network_bootstrap[intersect(which(Overlap_network_bootstrap[, 1] %in% as.matrix(ASD)), which(Overlap_network_bootstrap[, 2] %in% as.matrix(ASD))), ]
    Overlap_network_bootstrap_rewired_ASD <- Overlap_network_bootstrap_rewired[intersect(which(Overlap_network_bootstrap_rewired[, 1] %in% as.matrix(ASD)), which(Overlap_network_bootstrap_rewired[, 2] %in% as.matrix(ASD))), ]
    
# downstream [5] hub  ####[1] Result analysis: The lncRNA regulation in each cell type is unique####
    ## Identifying cell-specific hub lncRNAs    
    CSlncR_network_bootstrap_outdegree <- lapply(seq(CSlncR_network_bootstrap), function(i) degree(CSlncR_network_bootstrap_graph[[i]], mode="out"))
    hub_lncRNAs_bootstrap <- lapply(seq(CSlncR_network_bootstrap), function(i) names(sort(CSlncR_network_bootstrap_outdegree[[i]][which(CSlncR_network_bootstrap_outdegree[[i]]!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSlncR_network_bootstrap_outdegree[[i]]!=0)))])

    # GW16_to_GW23.5
    CSlncR_network_bootstrap_GW16_to_GW23.5_outdegree <- lapply(seq(CSlncR_network_bootstrap_GW16_to_GW23.5), function(i) degree(CSlncR_network_bootstrap_GW16_to_GW23.5_graph[[i]], mode="out"))
    hub_lncRNAs_bootstrap_GW16_to_GW23.5 <- lapply(seq(CSlncR_network_bootstrap_GW16_to_GW23.5), function(i) names(sort(CSlncR_network_bootstrap_GW16_to_GW23.5_outdegree[[i]][which(CSlncR_network_bootstrap_GW16_to_GW23.5_outdegree[[i]]!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSlncR_network_bootstrap_GW16_to_GW23.5_outdegree[[i]]!=0)))])
    
    ## Overlap of cell-type-specific hub lncRNAs across cell types    
    Overlap_hub_lncRNAs_bootstrap_conserved <- Overlap_hub(hub_lncRNAs_bootstrap, Intersect_num = round(length(hub_lncRNAs_bootstrap)*1.0))# default:1
    Overlap_hub_lncRNAs_bootstrap_union <- Overlap_hub(hub_lncRNAs_bootstrap, Intersect_num = 1)
    Overlap_hub_lncRNAs_bootstrap_two <- Overlap_hub(hub_lncRNAs_bootstrap, Intersect_num = 2)
    Overlap_hub_lncRNAs_bootstrap_rewired <- setdiff(Overlap_hub_lncRNAs_bootstrap_union, Overlap_hub_lncRNAs_bootstrap_two)
      
    # write.csv(Overlap_hub_lncRNAs_bootstrap_rewired,"~/data/tmp.csv")
    
# downstream [6] ASD-related hub       
    ## ASD-related cell-type-specific hub lncRNAs    
    hub_lncRNAs_bootstrap_ASD <- lapply(seq(hub_lncRNAs_bootstrap), function(i) hub_lncRNAs_bootstrap[[i]][which(hub_lncRNAs_bootstrap[[i]] %in% as.matrix(ASD))])
    
    ## ASD-related conserved and rewired hub lncRNAs    
    Overlap_hub_lncRNAs_bootstrap_conserved_ASD <- Overlap_hub_lncRNAs_bootstrap_conserved[which(Overlap_hub_lncRNAs_bootstrap_conserved %in% as.matrix(ASD))]    
    Overlap_hub_lncRNAs_bootstrap_rewired_ASD <- Overlap_hub_lncRNAs_bootstrap_rewired[which(Overlap_hub_lncRNAs_bootstrap_rewired %in% as.matrix(ASD))]
    
# downstream [7]  similarity matrix    
    ## Calculating similarity matrix of cell-type-specific lncRNA-mRNA regulatory network across cell types    
    CSlncR_network_bootstrap_Sim <- Sim.network(CSlncR_network_bootstrap, CSlncR_network_bootstrap, directed = TRUE)
    CSlncR_network_bootstrap_GW16_to_GW23.5_Sim <- Sim.network(CSlncR_network_bootstrap_GW16_to_GW23.5, CSlncR_network_bootstrap_GW16_to_GW23.5, directed = TRUE)
    
    ## Calculating similarity matrix of cell-type-specific hub lncRNAs across cell types    
    CSlncR_hub_bootstrap_Sim <- Sim.hub(hub_lncRNAs_bootstrap, hub_lncRNAs_bootstrap)
    CSlncR_hub_bootstrap_GW16_to_GW23.5_Sim <- Sim.hub(hub_lncRNAs_bootstrap_GW16_to_GW23.5, hub_lncRNAs_bootstrap_GW16_to_GW23.5)
    
    ## Calculating similarity matrix of expression across cell types
    CSlncR_expression_Sim <- cor(t(lncRNA_mRNA_scRNA_norm_filter))
    
    #load("sim.RData") 
    #save(CSlncR_network_bootstrap_GW16_to_GW23.5_Sim,CSlncR_hub_bootstrap_GW16_to_GW23.5_Sim,CSlncR_expression_Sim,file = "GW16_to_GW23.5_sim.RData")
    load("GW16_to_GW23.5_sim.RData") 
#######################################################################################################################

    
# downstream [8] cell type-specific lncRNAs analysis 细胞类型特异的lncRNAs分析
###### [2] Result analysis: The cell type-specific regulation across single‑cells of human brain region NCX  ######
## The cell type-specific lncRNAs analysis
    
    # 关于细胞类型及marker的来源：原始数据的文章
    # endothelia "LINC00339", "TRIM52-AS1"
    lncR_endothelia  <- c("LINC-MILR1-3","SLC38A3","LINC00152","RP11-401P9.4","MIR4435-1HG",
                          "LINC00339","RP11-483C6.1","AP000459.4","AC127904.2","RP11-161M6.2",
                          "RP11-417F21.1", "TRIM52-AS1","CTD-2081C10.7","RP11-296I10.3","RP11-532M24.1")
    
    # radial glia "LINC00943","MAGI2-AS3", "RUSC1-AS1"
    lncR_radial_glia <- c("Z83001.1", "RP11-731J8.2","LINC00943", "RP3-418C23.2", "RP11-1002K11.1",
                          "MAGI2-AS3","RP11-421L21.3", "LINC-FZD8-3", "LINC-FZD8-1", "LINC00263", 
                          "EIF3J-AS1", "LOC646329","LINC-KREMEN1-1","RUSC1-AS1","DGKK")
    
    # dividing radial glia "THAP9-AS1"
    lncR_dividing_radial_glia <- c("UHRF1", "CTRR-175P5.4","RP11-138A9.1", "RP11-849F2.9", "RP11-143K11.1",
                          "AC004447.2","SNORA59B", "CTC-503J8.6", "RP11-138A9.2", "RP11-95D17.1", 
                          "THAP9-AS1", "SNHG1","CTD-2017D11.1","RP11-58B17.2","DYNLL1-AS1")
    #放射状胶质细胞特异的lncRNAs分析
    lncR_radial_glia_family <- c("LINC00943","MAGI2-AS3", "RUSC1-AS1", "THAP9-AS1")
    
    # intermediate progenitors "DGCR11"
    lncR_intermediate_progenitors <- c("LINC-TMEM200C-1","RP11-798G7.8","RP11-351J23.1-AS1","RP3-326L13.3","CTD-2245E15.3",
                                       "C1orf132","AC084018.1","RP11-73O6.3","RP11-594N15.3","RP11-436D23.1",
                                       "AC083884.8","DGCR11","RP11-456K23.1","RP6-24A23.3","RP1-20C7.6")
    
    # newborn neurons "INHBA-AS1","MYT1L-AS1","KIF9-AS1",
    lncR_newborn_neurons <- c("RP5-1024G6.8","LINC-PTCHD2-3","RP11-513M16.8","RP11-661O13.1","RP11-524C21.2",
                              "RP11-356K23.1","LINC01105","INHBA-AS1","MYT1L-AS1","KIF9-AS1",
                              "RP11-1006G14.4","RP11-296O14.3","RP11-452L6.5","CTD-3099C6.9","RP11-452H21.4")
    
    # maturing excitatory "MIR137HG","PWAR6","SIK3-IT1","NAV2-AS3","DAPK1-IT1",
    lncR_maturing_excitatory_neurons <- c("MIR137HG","LINC00599","PWAR6","SIK3-IT1","RP11-53O19.3",
                                  "RP11-402L6.1","RP11-18I14.10","RP11-486F17.1","NAV2-AS3","DAPK1-IT1",
                                  "RP11-397O4.1","RP11-64K12.10","LINC00643","RP3-462E2.5","LINC-TMEM182-5")
    
    # inhibitory interneurons  "DLX6-AS1","SOX2-OT","MEG3", 
    lncR_inhibitory_interneurons <- c("DLX6-AS1","RP11-588P7.1","SOX2-OT","GS1-18A18.1","MEG3",
                                      "LINC-DKFZP761K2322-2","GRIP2","AC087393.1","LINC00966","RP11-450H6.3",
                                      "RP13-514E23.1","RP11-379H18.1","RP11-69E11.4","AC012358.8","LINC-TBCC-1")
    
    
    lncR_cell_type_ <- c("LINC00339", "TRIM52-AS1","LINC00943","MAGI2-AS3", "RUSC1-AS1",
                        "THAP9-AS1","DGCR11","INHBA-AS1","MYT1L-AS1","KIF9-AS1",
                        "MIR137HG","PWAR6","SIK3-IT1","NAV2-AS3","DAPK1-IT1","DLX6-AS1","SOX2-OT","MEG3")
    
    lncR_cell_type <- c(lncR_endothelia,lncR_radial_glia,lncR_dividing_radial_glia,lncR_intermediate_progenitors,lncR_newborn_neurons,lncR_maturing_excitatory_neurons,lncR_inhibitory_interneurons)

    # Extracting lncR_cell_type related cell-type-specific lncRNA-mRNA regulatory networks
    CSlncR_network_bootstrap_lncRcelltype <- lapply(seq(CSlncR_network_bootstrap), function(i) CSlncR_network_bootstrap[[i]][which(CSlncR_network_bootstrap[[i]][, 1] %in% as.matrix(lncR_cell_type)), ])
    
    # Extracting conserved and rewired lncRNA-mRNA regulatory networks associated with the lncR_cell_type
    Overlap_network_bootstrap_lncRcelltype <- Overlap_network_bootstrap[which(Overlap_network_bootstrap[, 1] %in% as.matrix(lncR_cell_type)), ]
    Overlap_network_bootstrap_rewired_lncRcelltype <- Overlap_network_bootstrap_rewired[which(Overlap_network_bootstrap_rewired[, 1] %in% as.matrix(lncR_cell_type)), ]
    
    # write.csv(Overlap_network_bootstrap_lncRcelltype,"~/data/tmp.csv")
    
    # Extracting ASD-related conserved and rewired lncRNA-mRNA regulatory networks associated with the lncR_cell_type
    CSlncR_network_bootstrap_lncRcelltype_ASD <- lapply(seq(CSlncR_network_bootstrap_lncRcelltype), function(i) CSlncR_network_bootstrap_lncRcelltype[[i]][intersect(which(CSlncR_network_bootstrap_lncRcelltype[[i]][, 1] %in% as.matrix(ASD)), which(CSlncR_network_bootstrap_lncRcelltype[[i]][, 2] %in% as.matrix(ASD))), ])
    
    # Experimentally validated lncRNA-mRNA interactions associated with the lncRcelltype
    CSlncR_network_bootstrap_lncRcelltype_graph <- lapply(seq(CSlncR_network_bootstrap_lncRcelltype), function(i) make_graph(c(t(CSlncR_network_bootstrap_lncRcelltype[[i]][, 1:2])), directed = TRUE))
    CSlncR_network_bootstrap_lncRcelltype_validated <- lapply(seq(CSlncR_network_bootstrap_lncRcelltype), function(i) as_data_frame(CSlncR_network_bootstrap_lncRcelltype_graph[[i]] %s% lncRTarget_graph))
    
    ## Difference of pvalue from the lncRcelltype regulation    
  CSlncR_pvalue <- function(tmp,cell_number) {
    res_pvalue <- matrix(NA, nrow = cell_number, ncol = cell_number)
    for (i in seq(cell_number)) {
      for (j in seq(cell_number)) {
      
        cell_1 <- tmp[[i]]
        cell_2 <- tmp[[j]]

        cell_1_ <- paste(cell_1[,1], cell_1[,2], sep = "&")
        cell_2_ <- paste(cell_2[,1], cell_2[,2], sep = "&")
        cell_1_2_union <- union(cell_1_,cell_2_)
        cell_1_2_intersect <- intersect(cell_1_,cell_2_)
        cell_1_remaining_num <-  length(cell_1_) -  length(cell_1_2_intersect)
        cell_2_remaining_num <-  length(cell_2_) -  length(cell_1_2_intersect)
        
        intersect_ <- rep(c(1),length(cell_1_2_intersect))
       
        cell_1_remaining_1_ <- rep(c(1),cell_1_remaining_num)
        cell_1_remaining_0_ <- rep(c(0),cell_1_remaining_num)
        
        cell_2_remaining_1_ <- rep(c(1),cell_2_remaining_num)
        cell_2_remaining_0_ <- rep(c(0),cell_2_remaining_num)
        
        cell_1_regulation <- c(intersect_,cell_1_remaining_1_,cell_2_remaining_0_)
        cell_2_regulation <- c(intersect_,cell_1_remaining_0_,cell_2_remaining_1_)
        
        res_cell1 <- cell_1_regulation
        res_cell2 <- cell_2_regulation

        res_pvalue[i, j] <- ks.test(res_cell1, res_cell2)$p.value
      }
    }
    
    return(res_pvalue) 
  } 
    
    res_predict_pvalue <- CSlncR_pvalue(tmp= CSlncR_network_bootstrap_lncRcelltype,cell_number=cell_number) 
    res_validated_pvalue <- CSlncR_pvalue(tmp= CSlncR_network_bootstrap_lncRcelltype_validated,cell_number=cell_number) 
    res_ASD_pvalue <- CSlncR_pvalue(tmp= CSlncR_network_bootstrap_lncRcelltype_ASD,cell_number=cell_number) 

    rownames(res_predict_pvalue) <- colnames(res_predict_pvalue) <- paste("Cell",c(1:cell_number),sep=" ")
    rownames(res_validated_pvalue) <- colnames(res_validated_pvalue) <- paste("Cell",c(1:cell_number),sep=" ")
    rownames(res_ASD_pvalue) <- colnames(res_ASD_pvalue) <- paste("Cell",c(1:cell_number),sep=" ")
    test.value <- matrix(0, nrow = cell_number, ncol = cell_number)
    rownames(test.value) <- colnames(test.value) <- c("GW16","GW19.5","GW20.5", "GW21", "GW23.5")
    
    library(corrplot)
    par(mfrow=c(1,2)) # 形成1行、2列的图形矩阵
    corrplot(test.value, p.mat = res_predict_pvalue, method = "square", diag = FALSE, type = "upper", title="A  Difference in predicted targets",
             cl.pos="n", sig.level = c(.001, .01, .05), pch.cex = 1, mar=c(0,0,1.5,0), insig = "label_sig", pch.col = "black")
    
    corrplot(test.value, p.mat = res_ASD_pvalue, method = "square", diag = FALSE, type = "upper", title="B  Difference in ASD-related targets",
             cl.pos="n", sig.level = c(.001, .01, .05), pch.cex = 1, mar=c(0,0,1.5,0), insig = "label_sig", pch.col = "black")    
    
    corrplot(test.value, p.mat = res_validated_pvalue, method = "square", diag = FALSE, type = "upper", title="C  Difference in validated targets of cell type-specific",
             cl.pos="n", sig.level = c(.001, .01, .05), pch.cex = 1, mar=c(0,0,1.5,0), insig = "label_sig", pch.col = "black")
    
    pred_ <- data.frame(id=c("conserved","rewired"),value=c(dim(Overlap_network_bootstrap_lncRcelltype)[1],dim(Overlap_network_bootstrap_rewired_lncRcelltype)[1]))
    ggplot(pred_, aes(x = id, y = value)) +
      geom_bar(fill = ifelse(pred_$id=="conserved","blue","red"), stat = "identity", width = 0.6) +
      xlab("lncRNA regulation type") +
      ylab("#Predicted targerts of cell type-specific") +
      geom_text(aes(label=value,y=value+8),size=3) #+
      #labs(title='C') # GW16: A, GW19.5: B, GW20.5: C, GW21: D, GW23.5: E
    
    # each plot is  5.5 5.5           550 550
    
    ## Enrichment analysis of conserved and rewired the lncRcelltype regulation   
    CSlncR_network_lncRcelltype_bootstrap_conserved = list(Overlap_network_bootstrap_lncRcelltype)
    CSlncR_network_lncRcelltype_bootstrap_rewired = list(Overlap_network_bootstrap_rewired_lncRcelltype)
    
    lncRcelltype_bootstrap_conserved_gene_list=c()
    for (i in seq(nrow(CSlncR_network_lncRcelltype_bootstrap_conserved[[1]]))) {
      for (j in seq(ncol(CSlncR_network_lncRcelltype_bootstrap_conserved[[1]]))) {
        lncRcelltype_bootstrap_conserved_gene_list=append(x=lncRcelltype_bootstrap_conserved_gene_list, CSlncR_network_lncRcelltype_bootstrap_conserved[[1]][i,j])
      }
      
    }
    lncRcelltype_bootstrap_conserved_gene_list <- list(unique(lncRcelltype_bootstrap_conserved_gene_list))
    
    lncRcelltype_bootstrap_rewired_gene_list=c()
    for (i in seq(nrow(CSlncR_network_lncRcelltype_bootstrap_rewired[[1]]))) {
      for (j in seq(ncol(CSlncR_network_lncRcelltype_bootstrap_rewired[[1]]))) {
        lncRcelltype_bootstrap_rewired_gene_list=append(x=lncRcelltype_bootstrap_rewired_gene_list, CSlncR_network_lncRcelltype_bootstrap_rewired[[1]][i,j])
      }
      
    }
    lncRcelltype_bootstrap_rewired_gene_list <- list(unique(lncRcelltype_bootstrap_rewired_gene_list))
    
    # GO, KEGG and Reactome enrichment analysis
    lncRcelltype_bootstrap_conserved_FEA <- moduleFEA(lncRcelltype_bootstrap_conserved_gene_list, padjustvaluecutoff = 0.05)
    lncRcelltype_bootstrap_rewired_FEA <- moduleFEA(lncRcelltype_bootstrap_rewired_gene_list, padjustvaluecutoff = 0.05)
    
    # ASD enrichment analysis   
    lncRcelltype_bootstrap_conserved_ASD_EA <- module_ASD_EA(lncRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, ASD, lncRcelltype_bootstrap_conserved_gene_list)
    lncRcelltype_bootstrap_rewired_ASD_EA <- module_ASD_EA(lncRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, ASD, lncRcelltype_bootstrap_rewired_gene_list)
    
    # Hallmark enrichment analysis
    m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, human_gene_symbol)
    
    lncRcelltype_bootstrap_conserved_Hallmark <- lapply(seq(lncRcelltype_bootstrap_conserved_gene_list), function(i) enricher(lncRcelltype_bootstrap_conserved_gene_list[[i]], TERM2GENE=m_t2g, minGSSize=1) %>% as.data.frame)
    lncRcelltype_bootstrap_rewired_Hallmark <- lapply(seq(lncRcelltype_bootstrap_rewired_gene_list), function(i) enricher(lncRcelltype_bootstrap_rewired_gene_list[[i]], TERM2GENE=m_t2g, minGSSize=1) %>% as.data.frame)
    
    # Cell marker enrichment analsyis
    cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
      tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
      dplyr::select(cellMarker, geneSymbol)    
    
    lncRcelltype_bootstrap_conserved_Cellmarker <- lapply(seq(lncRcelltype_bootstrap_conserved_gene_list), function(i) enricher(lncRcelltype_bootstrap_conserved_gene_list[[i]], TERM2GENE=cell_markers, minGSSize=1) %>% as.data.frame)
    lncRcelltype_bootstrap_rewired_Cellmarker <- lapply(seq(lncRcelltype_bootstrap_rewired_gene_list), function(i) enricher(lncRcelltype_bootstrap_rewired_gene_list[[i]], TERM2GENE=cell_markers, minGSSize=1) %>% as.data.frame)
    
    #write.csv(lncRcelltype_bootstrap_conserved_FEA[[3]][[1]]@result,"~/data/tmp.csv")
    #write.csv(lncRcelltype_bootstrap_conserved_Hallmark[[1]],"~/data/tmp.csv")
    #write.csv(lncRcelltype_bootstrap_conserved_Cellmarker[[1]],"~/data/tmp.csv")
    
    # write.csv(lncRcelltype_bootstrap_rewired_FEA[[3]][[1]]@result,"~/data/tmp.csv")
    # write.csv(lncRcelltype_bootstrap_rewired_Hallmark[[1]],"~/data/tmp.csv")
################################################################################################################################################ 
    
    
# downstream [9] Similarity network plot   #### [1] Result analysis: The lncRNA regulation in each cell is unique ####   
## Similarity network plot in terms of cell-type-specific lncRNA-mRNA regulatory netowork    
    par(mfrow=c(1,2)) # 形成1行、2列的图形矩阵
    rownames(CSlncR_network_bootstrap_Sim) <- colnames(CSlncR_network_bootstrap_Sim) <- c("GW16","GW19.5","GW20.5","GW21","GW23.5")#rownames(lncRNA_scRNA_norm_filter) # c("GW16","GW19.5","GW20.5","GW21","GW23.5")#
    corrplot(CSlncR_network_bootstrap_Sim, method = "pie", type = "upper",title = "A  Cell-type-specific lncRNA-mRNA interactions", mar=c(0,0,2,0),diag = FALSE, cl.lim = c(0, 1), tl.cex = 1)
    
## Similarity network plot in terms of cell-type-specific hub lncRNAs
    rownames(CSlncR_hub_bootstrap_Sim) <- colnames(CSlncR_hub_bootstrap_Sim) <- c("GW16","GW19.5","GW20.5","GW21","GW23.5")#rownames(lncRNA_scRNA_norm_filter)# c("GW16","GW19.5","GW20.5","GW21","GW23.5")
    corrplot(CSlncR_hub_bootstrap_Sim, method = "pie", type = "upper", title =    "B  Cell-type-specific hub lncRNAs", diag = FALSE, mar=c(0,0,2,0),cl.lim = c(0, 1), tl.cex = 1)
    
    # save 7 7 inches-pdf, 700 700 tiff (single images)
    #     14 5            950 450 (GW16toGW23.5)
################################################################################################################################################    

    
# downstream [10] ggplot plot patchwork   #### [1] Result analysis: The lncRNA regulation in each cell type (development) is unique ####   
## Stem plots
#id_1=c(16,19.5,20.5,21,23.5)
id_=c("GW16","GW19.5","GW20.5","GW21","GW23.5")
# id_=seq(cell_number)
index_net <- data.frame(value = unlist(lapply(seq(CSlncR_network_bootstrap), function(i) nrow(CSlncR_network_bootstrap[[i]]))), id = id_)

col1 <- rep("#FF9999", cell_number)
p1 <- ggplot(index_net, aes(x = id, y = value)) +
  geom_point(aes(color = col1), size = 5) +
  geom_bar(aes(fill = col1), stat = "identity", width = 0.03) +
  xlab("Development Stages") +
  ylab("Predicted numbers") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position="none",
        panel.border = element_blank(),
        axis.text.x = element_text(face = "bold",size=10),
        axis.text.y = element_text(face = "bold",size=10),
        axis.title.x = element_text(face = "bold",size=12),
        axis.title.y = element_text(face = "bold",size=12)) 

index_validated_net <- data.frame(value = unlist(lapply(seq(CSlncR_network_bootstrap_validated), function(i) nrow(CSlncR_network_bootstrap_validated[[i]])/nrow(CSlncR_network_bootstrap[[i]])*100)), id = id_)

col2 <- rep("plum4", cell_number)
p2 <- ggplot(index_validated_net, aes(x = id, y = value)) +
  geom_point(aes(color = col2), size = 5) +
  geom_bar(aes(fill = col2), stat = "identity", width = 0.03) +
  scale_fill_manual(values=c("plum4"), aesthetics = "fill") +
  scale_colour_manual(values=c("plum4"), aesthetics = "colour") +
  xlab("Development Stages") +
  ylab("%Validated numbers") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position="none",
        panel.border = element_blank(),
        axis.text.x = element_text(face = "bold",size=10),
        axis.text.y = element_text(face = "bold",size=10),
        axis.title.x = element_text(face = "bold",size=12),
        axis.title.y = element_text(face = "bold",size=12))

index_ASD_net <- data.frame(value = unlist(lapply(seq(CSlncR_network_bootstrap_ASD), function(i) nrow(CSlncR_network_bootstrap_ASD[[i]])/nrow(CSlncR_network_bootstrap[[i]])*100)), id = id_)

col3 <- rep("blue", cell_number)
p3 <- ggplot(index_ASD_net, aes(x = id, y = value)) +
  geom_point(aes(color = col3), size = 5) +
  geom_bar(aes(fill = col3), stat = "identity", width = 0.03) +
  scale_fill_manual(values=c("blue"), aesthetics = "fill") +
  scale_colour_manual(values=c("blue"), aesthetics = "colour") +
  xlab("Development Stages") +
  ylab("%ASD-related numbers") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position="none",
        panel.border = element_blank(),
        axis.text.x = element_text(face = "bold",size=10),
        axis.text.y = element_text(face = "bold",size=10),
        axis.title.x = element_text(face = "bold",size=12),
        axis.title.y = element_text(face = "bold",size=12))

index_ASD_hub <- data.frame(value = unlist(lapply(seq(hub_lncRNAs_bootstrap_ASD), function(i) length(hub_lncRNAs_bootstrap_ASD[[i]])/length(hub_lncRNAs_bootstrap[[i]])*100)), id = id_)

col4 <- rep("green", cell_number)
p4 <- ggplot(index_ASD_hub, aes(x = id, y = value)) +
  geom_point(aes(color = col4), size = 5) +
  geom_bar(aes(fill = col4), stat = "identity", width = 0.03) +
  scale_fill_manual(values=c("green"), aesthetics = "fill") +
  scale_colour_manual(values=c("green"), aesthetics = "colour") +
  xlab("Development Stages") +
  ylab("%ASD-related numbers") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position="none",
        panel.border = element_blank(),
        axis.text.x = element_text(face = "bold",size=10),
        axis.text.y = element_text(face = "bold",size=10),
        axis.title.x = element_text(face = "bold",size=12),
        axis.title.y = element_text(face = "bold",size=12)) 

library(patchwork)
(p1+p2)/(p3+p4) + plot_annotation(tag_levels = 'A')
# save 8.5 5 inches-pdf, 850 500 tiff  # GW16, GW19.5, GW20.5, GW21, GW23.5
################################################################################################################################################ 


# downstream [11] comparison of CSlncR-Random-TargetScan   #### [3] Result analysis: CSlncR is effective in predicting cell-specific lncRNA targets  ####   
## Stem plots
index_validated_net <- data.frame(value = unlist(lapply(seq(CSlncR_network_bootstrap_validated), function(i) nrow(CSlncR_network_bootstrap_validated[[i]])/nrow(CSlncR_network_bootstrap[[i]])*100)), id = c("GW16","GW19.5","GW20.5", "GW21", "GW23.5"))
index_validated_net_ <- data.frame(value = unlist(lapply(seq(CSlncR_network_bootstrap_validated_), function(i) nrow(CSlncR_network_bootstrap_validated_[[i]])/nrow(CSlncR_network_bootstrap_[[i]])*100)), id = c("GW16","GW19.5","GW20.5", "GW21", "GW23.5"))
index_validated_net_random <-  data.frame(value = unlist(lapply(seq(CSlncR_network_bootstrap_validated_random), function(i) nrow(CSlncR_network_bootstrap_validated_random[[i]])/nrow(CSlncR_network_bootstrap_random[[i]])*100)), id = c("GW16","GW19.5","GW20.5", "GW21", "GW23.5"))
index_validated_net_LncRNA2Target <- data.frame(value = unlist(lapply(seq(CSlncR_network_bootstrap_validated_LncRNA2Target), function(i) nrow(CSlncR_network_bootstrap_validated_LncRNA2Target[[i]])/nrow(CSlncR_network_bootstrap_LncRNA2Target[[i]])*100)), id = c("GW16","GW19.5","GW20.5", "GW21", "GW23.5"))

data1 <- rbind(index_validated_net,index_validated_net_)
data1$Methods <- c(rep("CSlncR",cell_number), rep("CSlncR without prior knowledge",cell_number))
data2 <- rbind(index_validated_net,index_validated_net_random)
data2$Methods <- c(rep("CSlncR",cell_number), rep("Random",cell_number))
data3 <- rbind(index_validated_net,index_validated_net_LncRNA2Target)
data3$Methods <- c(rep("CSlncR",cell_number), rep("LncRNA2Target",cell_number))
  
p1 <- ggplot(data1,aes(x=id,y=value,group=Methods,color=Methods))+
  geom_line(size=0.5)+
  geom_point(aes(shape=Methods,color=Methods),size=2)+
  xlab("Development Stages") +
  ylab("%Validated numbers") +
  theme(legend.position="top",
        axis.text.x = element_text(face = "bold",size=10),
        axis.text.y = element_text(face = "bold",size=10),
        axis.title.x = element_text(face = "bold",size=12),
        axis.title.y = element_text(face = "bold",size=12)) 

p2 <- ggplot(data2,aes(x=id,y=value,group=Methods,color=Methods))+
  geom_line(size=0.5)+
  geom_point(aes(shape=Methods,color=Methods),size=2)+
  xlab("Development Stages") +
  ylab("%Validated numbers") +
  theme(legend.position="top",
        axis.text.x = element_text(face = "bold",size=10),
        axis.text.y = element_text(face = "bold",size=10),
        axis.title.x = element_text(face = "bold",size=12),
        axis.title.y = element_text(face = "bold",size=12)) 

p3 <- ggplot(data3,aes(x=id,y=value,group=Methods,color=Methods))+
  geom_line(size=0.5)+
  geom_point(aes(shape=Methods,color=Methods),size=2)+
  xlab("Development Stages") +
  ylab("%Validated numbers") +
  theme(legend.position="top",
      axis.text.x = element_text(face = "bold",size=10),
      axis.text.y = element_text(face = "bold",size=10),
      axis.title.x = element_text(face = "bold",size=12),
      axis.title.y = element_text(face = "bold",size=12)) 

library(patchwork)
(p1 + p2 + plot_layout(ncol=2)) / (p3) + plot_annotation(tag_levels = 'A') 
# save 10 7.2 inches-pdf, 1000 720 tiff
################################################################################################################################################ 


