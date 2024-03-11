### This script load most of the functions and libraries that I made or source for my DEG analysis
library(ggplot2)
library(magrittr)
library(crayon)
### Because Im lazy....
len <- function(x) (length(x)) # len instead of length
#not in function
'%not_in%' <- function(x,y)!('%in%'(x,y))

####
#detach stuborn packages
####
detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}




# Get the row names or a data frame that matches, a value useful for edgeR pipeline as the results from glmLRT (and other tests) are expressed in 
#the 1, 0, and -1, which are upregulated, not significant, or down regulated respectively
get.rownames <- function(x,y) rownames(x[x == y,]) 
# Flat violin plot function taken from https://datavizpyr.com/rain-cloud-plots-using-half-violin-plot-with-jittered-data-points-in-r/
# They take this plot from somewhere else I think.
# This function uses a matrix and plot its values according to a provided gene list
#
plot_z_raincloud_of_genes <- function(gene_list = STOP1_targets_in_low_Pi_up, matrix = mean_z_score_of_CPM){
  matrix_of_gene_list <- melt(matrix[which(rownames(matrix) %in% c(gene_list)),])
  matrix_of_gene_list <- 
    transform(matrix_of_gene_list,
              variable = factor(variable, levels = rev(colnames(matrix))))
  model <- lm(value~variable, data = matrix_of_gene_list)
  test <- HSD.test(model,"variable")
  test$groups <- test$groups[colnames(matrix),]
  
  ggplot(matrix_of_gene_list, aes(x = variable, y = value, fill = variable)) +
    geom_hline(yintercept = 0, size = 1.5) +
    geom_jitter(width=0.1, alpha = 0.7, aes(color = variable)) +
    geom_boxplot(width=0.4, alpha = 0.15,
                 aes(fill = variable, color = variable, color = after_scale(darken(color, .4, space = "HLS")))) +
    stat_summary(geom = 'text', label = rev(test$groups$groups), fun.y = max, vjust = 0.25, hjust = -0.1, size = 10) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = 0.7, color = "#00000000") +
    theme_light() +
    theme(axis.text=element_text(size=24),
                axis.title=element_text(size=24,face="bold"),
          legend.position = "top",
          legend.title = element_blank()) +
    ggtitle(deparse(substitute(gene_list))) +
    ylab("Z score") +
    xlab("") +
    coord_flip() +
    scale_fill_manual(values = c(viridis(2, option = "G", end = 0.6, begin = 0.4, direction = -1), viridis(2, option = "A", end = 0.6, begin = 0.4))) +
    scale_color_manual(values = c(viridis(2, option = "G", end = 0.6, begin = 0.4, direction = -1), viridis(2, option = "A", end = 0.6, begin = 0.4))) +
    scale_y_continuous(limits = c(-1.9, 1.9))
}
#### only for meta analisys 
edgeR_pipeline_per_study <- function(DGE_matrix_with_all_batches = complete_DGE_matrix,
                                      library_groups_with_all_batches = complete_library_groups,
                                     study = 51) {
  DGE_matrix_study_X <- NULL
  DGE_matrix_study_X$samples <- DGE_matrix_with_all_batches$samples[which(batch == study),]
  if (mean(DGE_matrix_study_X$samples$lib.size) == 0){
    warning(paste0("lib.size == 0 in study ", study, " consider remapping or removing this samples"))
    return(NULL)
  } else{
  DGE_matrix_study_X$counts <- DGE_matrix_with_all_batches$counts[,
                                                                  which(colnames(DGE_matrix_with_all_batches$counts) %in%
                                                                          rownames(DGE_matrix_study_X$samples))
                                                                  ]
  library_groups_study_X <- library_groups_with_all_batches[which(batch == study)]
  DGE_matrix_study_X <- DGEList(DGE_matrix_study_X$counts, group = library_groups_study_X)
  ###remove 0 libs
  DGE_matrix_study_X$counts <- DGE_matrix_study_X$counts[,
                                                         !colnames(DGE_matrix_study_X$counts) %in%
                                                           names(which(colSums(DGE_matrix_study_X$counts) == 0))
                                                         ]
  DGE_matrix_study_X$samples <- DGE_matrix_study_X$samples[which(DGE_matrix_study_X$samples$lib.size != 0),]
  ##normlib size for edgeR pipeline
  groups <- as.character(DGE_matrix_study_X$samples$group)
  filtered_DGE <- filterByExpr(DGE_matrix_study_X) #$counts, group = groups
  
  DGE_matrix_study_X <- DGE_matrix_study_X[filtered_DGE, , keep.lib.sizes=FALSE]
  DGE_matrix_study_X %<>% normLibSizes(method = "TMM")
  groups <- DGE_matrix_study_X$samples$group
  library_membership_matrix <- model.matrix(~0+groups)
  colnames(library_membership_matrix) <- levels(groups)
  DGE_matrix_study_X %<>% estimateDisp(library_membership_matrix)
  return(DGE_matrix_study_X)
  }
}

####Get All pairwise comparisons
library(gtools)
all_pairwise_vector <- function(count_matrix_groups = study_matrix$samples$group){
  levels(count_matrix_groups) <- unique(count_matrix_groups)
  all_combinations_matrix <- combinations(len(unique(count_matrix_groups)),
                                          2,
                                          levels(count_matrix_groups),
                                          repeats.allowed = FALSE)
  pairwise_comparisons_direction_1 <- paste0(all_combinations_matrix[,1],
                                             " - ",
                                             all_combinations_matrix[,2])
  pairwise_comparisons_direction_2 <- paste0(all_combinations_matrix[,2],
                                             " - ",
                                             all_combinations_matrix[,1])
  all_pairwise_comparisons <- append(pairwise_comparisons_direction_1, pairwise_comparisons_direction_2)
  return(all_pairwise_comparisons)
}
####
get_pairwise_DEGs <- function(count_matrix = DGE_matrix,
                              contrast_vector = contrast_to_make_vector,
                              fit = fit_model,
                              lfc = 0.2
                              ) {
  levels(count_matrix$samples$group) <- unique(count_matrix$samples$group)
  all_combinations_matrix <- combinations(len(unique(count_matrix$samples$group)),
                                          2,
                                          levels(count_matrix$samples$group),
                                          repeats.allowed = FALSE)
  contrast_result <- list()
  for(counter in 1:(len(contrast_vector))){
    contrast <- do.call(makeContrasts, list(contrast_vector[counter], levels=unique(count_matrix$samples$group)))
    glmTreat <- glmTreat(fit, contrast = contrast, lfc = lfc)
    if(counter <= len(contrast_vector)/2){
      comparison_name <- paste0("DEGs_in_", all_combinations_matrix[counter,1], "_relative_to_", all_combinations_matrix[counter,2])
      contrast_result[[comparison_name]] <- list(Comparison = comparison_name, glmTreat = glmTreat)
      text <- paste0("Done: ", crayon::bold(crayon::red(all_combinations_matrix[counter,1])), " VS ", crayon::bold(crayon::cyan(all_combinations_matrix[counter,2])), "\n")
      cat(text)
      }  else {
        comparison_name <- paste0("DEGs_in_", all_combinations_matrix[counter-len(contrast_vector)/2,2], "_relative_to_", all_combinations_matrix[counter-len(contrast_vector)/2,1])
        contrast_result[[comparison_name]] <- list(Comparison = comparison_name, glmTreat = glmTreat)
        text <- paste0("Done: ",crayon::bold(crayon::cyan(all_combinations_matrix[counter-len(contrast_vector)/2,2])), " VS ", crayon::bold(crayon::red(all_combinations_matrix[counter-len(contrast_vector)/2,1])), "\n")
        cat(text)
      }
    }
  return(contrast_result)
}





### This part load some data for ensemble so some functions works

mart <- biomaRt::useMart(biomart = "plants_mart",
                         dataset = "athaliana_eg_gene",
                         host = 'https://plants.ensembl.org')

# Get ensembl gene ids and GO terms
GTOGO <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                        "go_id", "name_1006", "tair_symbol", "namespace_1003"), mart = mart)
#Remove blank entries
GTOGO <- GTOGO[GTOGO$go_id != '',]
# convert from table format to list format
geneID2GO <- by(GTOGO$go_id,
                GTOGO$ensembl_gene_id,
                function(x) as.character(x))
all.genes <- sort(unique(as.character(GTOGO$ensembl_gene_id)))


#Gohistograms

mart <- biomaRt::useMart(biomart = "plants_mart",
                         dataset = "athaliana_eg_gene",
                         host = 'https://plants.ensembl.org')

# Get ensembl gene ids and GO terms
histogram.of.anotated.genes <-
function(mart = mart, ontology = c("biological_process", "molecular_function", "cellular_component")){
GTOGO <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                       "go_id", "name_1006", "tair_symbol", "namespace_1003"), mart = mart)
#Remove blank entries
GTOGO <- GTOGO[GTOGO$go_id != '',]
#select ontology
GTOGO <- filter(GTOGO, namespace_1003 == ontology)

# convert from table format to list format
geneID2GO <- by(GTOGO$go_id,
                GTOGO$ensembl_gene_id,
                function(x) as.character(x))
hist.data <- NULL
hist.data$ID <- rownames(geneID2GO)
no.of.GO.counter <- NULL
for(i in 1:dim(geneID2GO)) {
  no.of.GO.counter <- c(no.of.GO.counter,
  str_count(geneID2GO[i], "GO"))
}
hist.data$no.of.Terms <- no.of.GO.counter
hist.data <- as.data.frame(hist.data)
return(hist.data)

ggplot(hist.data, aes(x=hist.data$no.of.Terms)) + 
  geom_histogram(fill="lightgreen", color="grey50", bins=max(hist.data$no.of.Terms)+1) +
  ggtitle(ontology)

summarise(hist.data)
}
#### This function get all the genes annotated for a GO, this is essential of many other function I wrote
get_genes_of_a_go <- function(x = "GO:0006099") {
  cbind(GTOGO$ensembl_gene_id[GTOGO$go_id == x], GTOGO$tair_symbol[GTOGO$go_id == x], GTOGO$name_1006[GTOGO$go_id == x]) %>%
    as.data.frame() -> result_table
  colnames(result_table) <- c("ensembl_gene_id", "tair_symbol", "GO_name")
  if(len(c(which(duplicated(result_table$ensembl_gene_id)))) != 0){
    result_table[-c(which(duplicated(result_table$ensembl_gene_id))),] -> result_table
  } else{}
  return(result_table)
}
get_symbol_from_ID <- function(list = other_TFs$V2) {
  symbols <- c()
  cbind(GTOGO$ensembl_gene_id, GTOGO$tair_symbol) %>%
    as.data.frame() -> result_table
  for (i in 1:(list%>%len())){
    multi_symbol_handler <- c()
    result_table[which(list[i] == result_table$V1),] ->
    multi_symbol_handler
    if (is.na(multi_symbol_handler$V2[1] != TRUE) ) {
     c(symbols, list[i]) -> symbols
    } else if (multi_symbol_handler$V2[1] != "" ){
      c(symbols, multi_symbol_handler$V2[1]) -> symbols
    } else {
      c(symbols, multi_symbol_handler$V1[1]) -> symbols
    }
    }
  cbind(list, symbols)
}
##### This function make a heatmap of a GOs, it defalut to mean_z_score_of_CPM but you can add any matrix that have locus ID as rownames
GO.heatmap <- function(DE_matrix = mean_z_score_of_CPM, GO = "GO:0016036", option1 = "G", option2 ="A"){
genes_of_go <- get_genes_of_a_go(GO)
genes_of_go <- genes_of_go[order(genes_of_go$ensembl_gene_id), ]
DE_matrix <- DE_matrix[order(rownames(DE_matrix)), ]
genes_of_go_z_score <- DE_matrix[which(rownames(DE_matrix) %in% genes_of_go$ensembl_gene_id), ]
genes_of_go %<>% mutate(ID = paste0(ensembl_gene_id, "_" ,tair_symbol))
pheatmap(genes_of_go_z_score %>% t(),
         cluster_rows = FALSE,
         color = c(viridis(12, option = option1, end = 0.8, begin = 0.2, direction = -1), viridis(12, option = option2, end = 0.8, begin = 0.2)),
         main = genes_of_go$GO_name[1],
         border_color = "#00000000",
         labels_col = genes_of_go$ID[(which(genes_of_go$ensembl_gene_id %in%
                                              row.names(genes_of_go_z_score)))])
}

GO_PCA <- function( GO = "GO:0009737", DE_matrix = z_score_of_CPM){
  genes_of_go <- get_genes_of_a_go(GO)
  genes_of_go <- genes_of_go[order(genes_of_go$ensembl_gene_id), ]
  DE_matrix <- DE_matrix[order(rownames(DE_matrix)), ]
  genes_of_go_z_score <- DE_matrix[which(rownames(DE_matrix) %in% genes_of_go$ensembl_gene_id), ]
  genes_of_go %<>% mutate(ID = paste0(ensembl_gene_id, "_" ,tair_symbol))
  plotMDS(DE_matrix[which(rownames(DE_matrix) %in% genes_of_go$ensembl_gene_id), ], top = 20000)
}
#the cutoff is the percentage that is discarded ex: 0.75 means that the lower 75% get discarded

venn_of_treatments <- function(z_score_matrix = mean_z_score_of_CPM, list_of_genes = get_genes_of_a_go("GO:0009737")[[1]], cutoff = 0.75, upordown = "up") {
    gene_list_zscore <- z_score_matrix[which(rownames(z_score_matrix) %in% c(list_of_genes)),]
  for(counter in 1:dim(z_score_matrix)[[2]]){
    if (upordown == "up"){
    assign(paste0("treatment_", colnames(gene_list_zscore)[counter]), 
           gene_list_zscore[which(gene_list_zscore[[counter]] >=
                                    (max(mean_z_score_of_CPM)+sqrt(min(mean_z_score_of_CPM)^2)*(2*cutoff))-(max(mean_z_score_of_CPM)+sqrt(min(mean_z_score_of_CPM)^2))),])
    } else {
      assign(paste0("treatment_", colnames(gene_list_zscore)[counter]), 
             gene_list_zscore[which(gene_list_zscore[[counter]] <=
                                      (max(mean_z_score_of_CPM)+sqrt(min(mean_z_score_of_CPM)^2)*(2*cutoff))-(max(mean_z_score_of_CPM)+sqrt(min(mean_z_score_of_CPM)^2))),])
    }
    
    }
  list_handler <- list()
  names_handler <- c()
  for(counter in 1:dim(mean_z_score_of_CPM)[[2]]){
    list_handler <- append(list_handler,list(rownames(get(paste0("treatment_", colnames(gene_list_zscore)[counter])))))
    names_handler <- append(names_handler, paste0("treatment_", colnames(gene_list_zscore)[counter]))
  }
  names(list_handler) <- names_handler
  ggvenn_plot <- ggvenn(list_handler,
                        fill_alpha = 0.5,
                        stroke_color = "#00000033",
                        stroke_size = 1,
                        fill_color = viridis(dim(mean_z_score_of_CPM)[[2]], end = 0.95, begin = 0.075, option = "B"),
                        text_color = "#000000CC"
  )
  plot(ggvenn_plot)
  
}

# This functions return the enriched GOs in a gene list and apply some calculations to the table to get the recall,
# -log of p val, and some color thresholds for ploting and better visualization. default ontology is BP but you can select CC or MF
get_enriched_go_table <- function(genelist, aspect = "BP") {
  int.genes <- factor(as.integer(all.genes %in% genelist))
  names(int.genes) <- all.genes
  go.obj <- new("topGOdata", ontology = aspect
                , allGenes = int.genes
                , annot =annFUN.gene2GO
                , gene2GO = geneID2GO)
  results <- runTest(go.obj, algorithm = "elim", statistic = "fisher")
  results.tab <- GenTable(object = go.obj, elimFisher = results, topNodes = 100)
  results.tab$Term <- str_extract(results.tab$Term, ".{1,45}")
  results.tab$Term[which(duplicated(results.tab$Term))] <- paste0(results.tab$Term[which(duplicated(results.tab$Term))],".")
  results.tab <- mutate(results.tab,
                        "logpval" = -log10(as.numeric(results.tab$elimFisher)),
                        "recall" = results.tab$Significant/results.tab$Annotated)
  results.tab <- results.tab[rev(order(results.tab$logpval)),]
  results.tab$logpval[which(is.na(results.tab$logpval))] <- max(results.tab$logpval[-which(is.na(results.tab$logpval))])+2 
  results.tab$Term <- factor(results.tab$Term, levels=rev(unique(results.tab$Term)))
  return(results.tab)
}

boxplot_of_gene <- function(data = DGE_matrix$counts, gene = "AT2G35690"){
  matrix <- c()
  matrix <- melt(data)
  #matrix$X2 %<>% str_extract("[^\\d]{1,6}")
  matrix <- matrix[which(matrix[1] == gene),]
  if (identical(which(matrix[1] == gene), integer(0)) == TRUE){
    stop("Gene not found, maybe it was filted out by you analisys?")
  }
  matrix$Var2 <- str_sub(matrix$Var2, end = -2)
  model <- lm(value~matrix$Var2, data = matrix)
  test <- HSD.test(model,"matrix$Var2")
  test$groups <- test$groups[c("hp", "hp_phi","lp", "lp_phi"),]
  ggplot(matrix, aes(x = Var2, y = value, fill = Var2)) +
    geom_boxplot(alpha = 0.5, lwd=1.5) +
    stat_summary(geom = 'text', label = (test$groups$groups), fun.y = max, vjust = -0.4, hjust = 0.5, size = 10) +
    #scale_color_manual(values=c("#050505", "#550000")) +
    scale_fill_manual(values = viridis(len(levels(factor(matrix$Var2))), option = "B")) +
    xlab(paste(matrix$Tags[1], "_", GTOGO$tair_symbol[which(GTOGO$ensembl_gene_id == gene)][1])) + ylab("cpm") +
    theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),
          axis.line=element_line(color="#333333", size = 1),
          axis.ticks.y = element_line(color = "#222222", size = 1),
          axis.ticks.x = element_blank(),
          axis.ticks.length = unit(.2, "cm"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.ontop = FALSE,
          panel.background = element_rect(fill = "transparent"),
          axis.text = element_text(face="bold", color = "#444444"),
          axis.text.x = element_text(face="bold", color = "#444444",
                                     size = 28, vjust = 0.5),
          axis.text.y = element_text(face="bold", color = "#444444", size = 28),
          axis.title.x = element_text(size = 32, vjust = 0.5),
          axis.title.y = element_text(size = 32, vjust = 0.5),
          legend.key.size = unit(1, 'cm')
    )
}


### this  function plot the enriched GOs ordered by  p val. The up or down only change the color scheme
bubbleplot <- function(data = gos_lpVShp_up, color = "D", Terms = 30) {
      GOs_plot <- 
      ggplot(
      data[1:Terms,],
      aes(x=Term, y=logpval)
    ) + 
      scale_color_viridis_c(option = color, begin = 0.4,end =  0.8) +
      geom_point(aes(color = recall, size = log10(Significant)), alpha = 0.9) +
        scale_size(range = c(2.5, 10)) +
      stat_summary(geom = 'text', label = c(paste0(data$Annotated[Terms:1],"/", data$Significant[Terms:1])),
                   fun.y = max, vjust = 0.25, hjust = -0.3, size = 5) +
      theme_light() +
        theme(plot.margin = margin(0.5,0.5,0.5,0, "mm"),
              axis.line=element_line(color="#333333", size = 1),
              axis.ticks.y = element_line(color = "#222222", size = 1),
              axis.ticks.x = element_blank(),
              axis.ticks.length = unit(.2, "cm"),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.ontop = FALSE,
              panel.background = element_rect(fill = "transparent"),
              axis.text.x = element_text(face="bold", color = "#444444",
                                         size = 15, vjust = 0.5),
              axis.text.y = element_text(face="plain", color = "#444444", size = 18),
              axis.title.x = element_text(size = 16, vjust = 0.5),
              axis.title.y = element_blank(),
              legend.position = "top",
              legend.text=element_text(size=12.5),
              legend.title=element_text(size=20)
              #legend.justification = "left",
        ) +
      ylim((min(data$logpval[1:Terms])-1), max(data$logpval[1:Terms])*1.6) +
        coord_flip() +
        guides(size = "none")
    return(GOs_plot)
}

volcano_plot <- function(lrt_table = lpVSlp_phi_lrt, option1 = "A", option2 = "D", n_of_scales = 2,
                         scales_breaks = c(0.8,0.2), color_breaks = c(0, 0.225, 0.45, 0.675, 0.9), reverse_scales = 1, n_tags = 10){
  
  data <- as.data.frame(lrt_table$table)
  data$upordown <- decideTests(lrt_table)
  data$PValue[which(data$PValue == "Inf")] <- min(data$PValue[which(data$PValue == 0)])
  data %<>% mutate("-log10_pval" = -log10(PValue))
  data$`-log10_pval`[which(data$`-log10_pval` == "Inf")] <- max(data$`-log10_pval`[-which(data$`-log10_pval` == "Inf")]) + 20
  #data %<>% mutate("euclidean_distance" = sqrt((`-log10_pval`^2)+(logFC^2))) #c^2=a^2+b^2
  data %<>% mutate("color_index" = (logFC*(`-log10_pval`*0.1)))
  
  
  data <- data[order(data$color_index, decreasing = TRUE), ]
  
  if (n_of_scales == 2){
    data$color <- 
      data$color <- 
      c(
        viridis(len(data$color_index[which(data$color_index >= 0)]) - round(len(data$color_index[which(data$color_index >= 0)])*scales_breaks[1]),
                begin = color_breaks[4], end = color_breaks[5], direction = -1 * reverse_scales, option = option1),
        viridis(len(data$color_index[which(data$color_index >= 0)]) - round(len(data$color_index[which(data$color_index >= 0)])*scales_breaks[2]),
                end = color_breaks[4], direction = -1 * reverse_scales, option = option1),
        viridis(len(data$color_index[which(data$color_index < 0)]) - round(len(data$color_index[which(data$color_index < 0)])*scales_breaks[2]),
                end = color_breaks[4], option = option2, direction = 1 * reverse_scales),
        viridis(len(data$color_index[which(data$color_index < 0)]) - round(len(data$color_index[which(data$color_index < 0)])*scales_breaks[1]),
                begin = color_breaks[4], end = color_breaks[5], option = option2, direction = 1 * reverse_scales)
      ) 
  } else if (n_of_scales == 1){
    data$color <-
      c(
        viridis(len(data$color_index[which(data$color_index >= 0)]) - round(len(data$color_index[which(data$color_index >= 0)])*scales_breaks[1]),
                begin = color_breaks[4], end = color_breaks[5], direction = -1 * reverse_scales, option = option1),
        viridis(len(data$color_index[which(data$color_index >= 0)]) - round(len(data$color_index[which(data$color_index >= 0)])*scales_breaks[2]),
                begin = color_breaks[3], end = color_breaks[4], direction = -1 * reverse_scales, option = option1),
        viridis(len(data$color_index[which(data$color_index < 0)]) - round(len(data$color_index[which(data$color_index < 0)])*scales_breaks[2]),
                begin = color_breaks[2], end = color_breaks[3], direction = -1 * reverse_scales, option = option1),
        viridis(len(data$color_index[which(data$color_index < 0)]) - round(len(data$color_index[which(data$color_index < 0)])*scales_breaks[1]),
                begin = color_breaks[1], end = color_breaks[2], direction = -1 * reverse_scales, option = option1)
      )
  }else {
    stop("n_of_scales must be 1 or 2")
  }
  data$diff <- "No. sig"
  data$diff[data$logFC > min(data$logFC[data$upordown == 1]) & 
              data$PValue < max(data$PValue[data$upordown == 1])] <- "Upregulated"
  data$diff[data$logFC < max(data$logFC[data$upordown == -1]) &
              data$PValue < max(data$PValue[data$upordown == -1])] <- "Downregulated"
  
  data$color[which(data$diff == "No. sig")] <- "#666666"
  
  
  upregulated <- data[which(data$diff == "Upregulated"),]
  downregulated <- data[which(data$diff == "Downregulated"),] 
  upregulated <- upregulated[order(upregulated$logFC, decreasing = TRUE), ]
  downregulated <- downregulated[order(downregulated$logFC, decreasing = FALSE), ]
  upregulated <- upregulated[1:n_tags,]
  downregulated <- downregulated[1:n_tags,]
  upregulated$ID <- upregulated %>% rownames()
  downregulated$ID <- downregulated %>% rownames()
  ggplot(data=data, aes(x=logFC, y=`-log10_pval`)) +
    #scale_color_manual(values = c(data$color))+
    labs(color = "") +
    geom_point(color = data$color) + 
    geom_text(aes(
      y = max(data$`-log10_pval`)*1.1, x = -max(data$logFC) , label = paste(summary(decideTests(lrt_table))[1], "downregulated")),
      size = 5) +
    geom_text(aes(
      y = max(data$`-log10_pval`)*1.1, x = max(data$logFC) , label = paste(summary(decideTests(lrt_table))[3], "upregulated")),
      size = 5) +
    geom_vline(xintercept=(min(data$logFC[data$upordown == 1])), col="#00000088", size = 1.2, linetype="dotted") +
    geom_vline(xintercept=(max(data$logFC[data$upordown == -1])), col="#00000088", size = 1.2, linetype="dotted") +
    geom_hline(yintercept=-log10(max(data$PValue[data$upordown == -1])), col="#00000088", size = 1.25, linetype="dotted") +
    geom_text_repel(aes(label = upregulated$ID, x=upregulated$logFC, y=upregulated$`-log10_pval`),
                    data = upregulated, max.overlaps = 1000,
                    ylim = c(0, NA),
                    xlim = c(max(data$logFC)*1,  max(data$logFC)*2),
                    #nudge_y = min(data$`-log10_pval`),
                    #force = 2,
                    #force_pull = 0,
                    min.segment.length = 0,
                    #point.padding = 0.2,
                    #box.padding = 0.3,
                    segment.linetype = 1,
                    #segment.curvature = 0.025,
                    #segment.ncp = 3,
                    #segment.angle = 70,
                    segment.size = 0.5,
                    segment.color = "#00000055",
                    #arrow = arrow(length = unit(0.010, "npc")),
                    max.iter = 10000000,
                    max.time = 3,
                    #hjust = 0,
                    direction = "y",
                    color = "#000000aa",
                    size = 4.5,
                    seed = 42,
                    
      ) +
    geom_text_repel(aes(label = downregulated$ID, x=downregulated$logFC, y=downregulated$`-log10_pval`),
                    data = downregulated, max.overlaps = 1000,
                    ylim = c(0, NA),
                    xlim = c(-max(data$logFC)*2,  -max(data$logFC)*1),
                    #nudge_y = min(data$`-log10_pval`),
                    #force = 2,
                    #force_pull = 0,
                    min.segment.length = 0,
                    #point.padding = 0.2, 
                    #box.padding = 0.3,
                    segment.linetype = 1,
                    #segment.curvature = -0.025,
                    #segment.ncp = 3,
                    #segment.angle = 70,
                    segment.size = 0.5,
                    segment.color = "#00000055",
                    #arrow = arrow(length = unit(0.010, "npc")),
                    max.iter = 10000000,
                    max.time = 3,
                    #hjust = 1,
                    direction = "y",
                    color = "#000000aa",
                    size = 4.5,
                    seed = 42,
                    
    ) +
    theme_light() +
    theme(plot.margin = margin(0.5,0.5,0.5,0, "mm"),
          axis.line=element_line(color="#333333", size = 1),
          axis.ticks.y = element_line(color = "#222222", size = 1),
          axis.ticks.x = element_blank(),
          axis.ticks.length = unit(.2, "cm"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.ontop = FALSE,
          panel.background = element_rect(fill = "transparent"),
          axis.text.x = element_text(face="bold", color = "#444444",
                                     size = 12, vjust = 0.5),
          axis.text.y = element_text(face="bold", color = "#444444", size = 14),
          axis.title.x = element_text(size = 16, vjust = 0.5),
          axis.title.y = element_text(size = 16, vjust = 0.5)
    ) +
    xlab("Log2FC") + ylab("-log10 p-value") +
    xlim(-max(data$logFC)*2, max(data$logFC)*2) 
    #ylim(-20, max(data$`-log10_pval`)*1.3)
}

# copied from
# https://gist.github.com/dgrtwo/eb7750e74997891d7c20

# somewhat hackish solution to:
# https://twitter.com/EamonCaddigan/status/646759751242620928
# based mostly on copy/pasting from ggplot2 geom_violin source:
# https://github.com/hadley/ggplot2/blob/master/R/geom-violin.r



"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x,
                     xmax = x + width / 2)
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                             plyr::arrange(transform(data, x = xmaxv), -y))
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1,])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),
          
          required_aes = c("x", "y")
  )

#### this one plot an histogram of the No of interactions
histogram.of.tf.interactions <- function(tf.data = tf.data, top.n = 30, option = "A", begin = 0.8,end =  0.8){
read.csv("C:/Users/Admin/OneDrive - Texas Tech University/Phi arabidopsis/Phi DEG/lp_network/Ath_TF_list.txt",
         sep = "\t") -> tf.family
tf.data <- tf.data[1:top.n,]
tf.data$TF <- factor(tf.data$TF, levels = tf.data$TF)
#factor(x$name, levels = x$name[order(x$val)])
IDs <- GTOGO[which(GTOGO$ensembl_gene_id %in% tf.data$TF[1:top.n]),] %>% mutate(ID = paste0(ensembl_gene_id, "_", tair_symbol))
IDs <- IDs[,c(1,5)]
IDs <- IDs[match(unique(IDs$ensembl_gene_id), IDs$ensembl_gene_id),]
IDs %<>% arrange(factor(ensembl_gene_id, levels = tf.data$TF))
non.id.genes <- as.character(tf.data[which(!(tf.data$TF %in% IDs$ensembl_gene_id)),2])
if(len(non.id.genes) != 0){
  for (i in 1:len(non.id.genes)){
    dplyr::add_row(
      IDs,
      ensembl_gene_id = non.id.genes[i],
      ID = non.id.genes[i],
      .before = which(!(tf.data$TF %in% IDs$ensembl_gene_id))
    ) -> IDs
  }
}
cbind(tf.data, IDs$ID) -> tf.data
tf.data$`IDs$ID` <- factor(tf.data$`IDs$ID`, levels = tf.data$`IDs$ID`)
colnames(tf.data)[3] <- "ID"
tf.family[which(tf.family$Gene_ID %in% tf.data$TF),] -> tf.family.in.data
tf.family.in.data %<>% arrange(factor(Gene_ID, levels = tf.data$TF))
tf.family.in.data[-which(duplicated(tf.family.in.data$Gene_ID)),] -> tf.family.in.data
cbind(tf.data, tf.family.in.data$Family) -> tf.data
colnames(tf.data)[4] <- "Family"
ggplot(tf.data, aes(x = ID, y=Counts, fill = Counts)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_c(option = option, begin = begin, end =  end) +
  stat_summary(geom = 'text', label = tf.data$Family, fun.y = max, vjust = 0.5, hjust = 1.1, size = 5, angle = 90) +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),
        axis.line=element_line(color="#222222", size = 1),
        axis.ticks.y = element_line(color = "#222222", size = 1),
        axis.ticks.length = unit(.2, "cm"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.ontop = FALSE,
        panel.background = element_rect(fill = "transparent"),
        axis.text = element_text(color = "#222222"),
        axis.text.x = element_text(color = "#222222",
                                   size = 16, vjust = 0.5, angle = 90),
        axis.text.y = element_text(color = "#222222", size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 26, vjust = 0.5),
        legend.position = "none"
  ) + ylab("Interactions")
}

