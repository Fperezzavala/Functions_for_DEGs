### This script load most of the functions that I made or source for my DEG analysis

### Because Im lazy....
len <- function(x) (length(x)) # len instead of length
#not in function
'%not_in%' <- function(x,y)!('%in%'(x,y))

# Get the row names or a data frame that matches, a value useful for edgeR pipeline as the results from glmLRT (and other tests) are expressed in 
#the 1, 0, and -1, which are upregulated, not significant, or down regulated respectively
get.rownames <- function(x,y) rownames(x[x == y,]) 
# Flat violin plot function taken from https://datavizpyr.com/rain-cloud-plots-using-half-violin-plot-with-jittered-data-points-in-r/
# They take this plot from somewhere else I think.
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R")

# This function uses a matrix and plot its values acording to a provided gene list
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


### This part load some data for ensemble so some functions works

mart <- biomaRt::useMart(biomart = "plants_mart",
                         dataset = "athaliana_eg_gene",
                         host = 'https://plants.ensembl.org')

# Get ensembl gene ids and GO terms
GTOGO <- biomaRt::getBM(attributes = c( "ensembl_gene_id",
                                        "go_id", "name_1006", "tair_symbol"), mart = mart)
#Remove blank entries
GTOGO <- GTOGO[GTOGO$go_id != '',]
# convert from table format to list format
geneID2GO <- by(GTOGO$go_id,
                GTOGO$ensembl_gene_id,
                function(x) as.character(x))
all.genes <- sort(unique(as.character(GTOGO$ensembl_gene_id)))


#### This function get all the genes annotated for a GO, this is essential of many other function I wrote
get_genes_of_a_go <- function(x = "GO:0009862") {
  cbind(GTOGO$ensembl_gene_id[GTOGO$go_id == x], GTOGO$tair_symbol[GTOGO$go_id == x], GTOGO$name_1006[GTOGO$go_id == x]) %>%
    as.data.frame() -> result_table
  colnames(result_table) <- c("ensembl_gene_id", "tair_symbol", "GO_name")
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
GO.heatmap <- function(DE_matrix = mean_z_score_of_CPM, GO = "GO:0009862", option1 = "G", option2 ="A"){
genes_of_go <- get_genes_of_a_go(GO)
genes_of_go <- genes_of_go[order(genes_of_go$ensembl_gene_id), ]
DE_matrix <- DE_matrix[order(rownames(DE_matrix)), ]
genes_of_go_z_score <- DE_matrix[which(rownames(DE_matrix) %in% genes_of_go$ensembl_gene_id), ]
genes_of_go %<>% mutate(ID = paste0(ensembl_gene_id, "_" ,tair_symbol))
pheatmap(DE_matrix[which(rownames(DE_matrix) %in% genes_of_go$ensembl_gene_id), ] %>% t(),
         cluster_rows = FALSE,
         color = c(viridis(12, option = option1, end = 0.8, begin = 0.2, direction = -1), viridis(12, option = option2, end = 0.8, begin = 0.2)),
         main = genes_of_go$GO_name[1],
         border_color = "#00000000",
         labels_col = genes_of_go$ID)
}

GO_PCA <- function( GO = "GO:0009737", DE_matrix = z_score_of_CPM){
  genes_of_go <- get_genes_of_a_go(GO)
  genes_of_go <- genes_of_go[order(genes_of_go$ensembl_gene_id), ]
  DE_matrix <- DE_matrix[order(rownames(DE_matrix)), ]
  genes_of_go_z_score <- DE_matrix[which(rownames(DE_matrix) %in% genes_of_go$ensembl_gene_id), ]
  genes_of_go %<>% mutate(ID = paste0(ensembl_gene_id, "_" ,tair_symbol))
  plotMDS(DE_matrix[which(rownames(DE_matrix) %in% genes_of_go$ensembl_gene_id), ], top = 20000)
}
#the cutoff is the percentage that is discarted ex: 0.75 means that the lower 75% get discarted

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
  results.tab <- mutate(results.tab,
                        "logpval" = -log10(as.numeric(results.tab$elimFisher)),
                        "recall" = results.tab$Significant/results.tab$Annotated)
  results.tab <- results.tab[rev(order(results.tab$logpval)),]
  results.tab$logpval[which(is.na(results.tab$logpval))] <- max(results.tab$logpval[-which(is.na(results.tab$logpval))])+2 
  results.tab$Term <- factor(results.tab$Term, levels=rev(unique(results.tab$Term)))
  return(results.tab)
}

boxplot_of_gene <- function(data = DGE_matrix$counts, gene = "AT1G64280"){
  matrix <- c()
  matrix <- melt(data)
  matrix$X2 %<>% str_extract("[^\\d]{1,6}")
  matrix <- matrix[which(matrix$X1 == gene),]
  if (identical(which(matrix$X1 == gene), integer(0)) == TRUE){
    stop("Gene not found, maybe it was filted out by you analisys?")
  }
  model <- lm(value~X2, data = matrix)
  test <- HSD.test(model,"X2")
  test$groups <- test$groups[c("hp", "hp_phi","lp", "lp_phi"),]
  ggplot(matrix, aes(x = X2, y = value, fill = X2)) +
    geom_boxplot(alpha = 0.5, lwd=1.5) +
    stat_summary(geom = 'text', label = (test$groups$groups), fun.y = max, vjust = -0.4, hjust = 0.5, size = 10) +
    #scale_color_manual(values=c("#050505", "#550000")) +
    scale_fill_manual(values = viridis(len(levels(factor(matrix$X2))), option = "B")) +
    xlab(paste(matrix$Tags[1], "_", GTOGO$tair_symbol[which(GTOGO$ensembl_gene_id == gene)][1])) + ylab("cpm") +
    theme(plot.margin = margin(0.5,0.5,0.5,0, "mm"),
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
      stat_summary(geom = 'text', label = paste0(data$Annotated[Terms:1],"/", data$Significant[Terms:1]),
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



