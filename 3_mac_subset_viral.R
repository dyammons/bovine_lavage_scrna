#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")
library(scProportionTest)

######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   customize variables   ######## <<<<<<<<<<<<<<
######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#set output param
outName <- "mac_viral"
reduction <- "umap.integrated.harmony"
reduction2 <- "integrated.harmony" 
clusMain <- "clusterID_integrated.harmony"
contrast <- c("Diseased", "Healthy") #code will test first vs second

#load in the proprocessed data & subset on pop of interest
seu.obj <- readRDS("../output/s3/20240708_bov_lav_n5n5_dxVSh_viral_S3.rds")
seu.obj <- cleanMeta(seu.obj)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_viral_ID.csv", groupBy = "clusterID_integrated.harmony", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident", metaAdd = "cellSource")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj <- subset(seu.obj, subset = majorID == "Macrophage") #specify subset of interest

######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin preprocessing   ######## <<<<<<<<<<<<<<
######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#integrate the data using all of the four Seurat v5 integration methods
seu.obj <- integrateData(seu.obj = seu.obj, dout = "../output/s2/", outName = outName, vars.to.regress = "percent.mt",
                         runAllMethods = TRUE, normalization.method = "LogNormalize", indReClus = T)

#use clustree to identify clustering parameters that appear most appropriate
clusTree(seu.obj = seu.obj, dout = "../output/clustree/", outName = outName, 
            test_dims = 30, algorithm = 3, prefix = "RNA_snn_res.", reduction = reduction2)

#complete data visualization
for (x in list("integrated.cca", "integrated.harmony", "integrated.joint", "integrated.rcpa")) {
    seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = paste0(outName, "_", x), 
                           final.dims = 30, final.res = 0.6, stashID = "clusterID", algorithm = 3, min.dist = 0.3, n.neighbors = 30,
                           prefix = "RNA_snn_res.", assay = "RNA", reduction = x,
                           saveRDS = F, return_obj = T, returnFeats = T,
                           features = c("PTPRC", "CD3E", "CD163", "CD14", 
                                        "CD40", "EPCAM", "ND5", "HLA-DRA", 
                                        "JCHAIN", "MYO1E", "CXCL2","CD69")
                          )
}

#remove low quality cells (c10; very high MT enrichment) and contaminating B cells (c14) - note c1 in this does look a little susspect
seu.obj <- subset(seu.obj, invert = T, subset = clusterID_integrated.harmony %in% c(10, 14))


#integrate the data using all of the four Seurat v5 integration methods
seu.obj <- integrateData(seu.obj = seu.obj, dout = "../output/s2/", outName = outName, 
                         runAllMethods = TRUE, normalization.method = "LogNormalize", indReClus = T)

#use clustree to identify clustering parameters that appear most appropriate
clusTree(seu.obj = seu.obj, dout = "../output/clustree/", outName = outName, 
            test_dims = 30, algorithm = 3, prefix = "RNA_snn_res.", reduction = "integrated.harmony")

#complete data visualization
# for (x in list("integrated.cca", "integrated.harmony", "integrated.joint", "integrated.rcpa")) {
for (x in list("integrated.harmony")) {
    seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = paste0(outName, "_", x), 
                           final.dims = 30, final.res = 1.2, stashID = "clusterID_sub", algorithm = 3, min.dist = 0.3, n.neighbors = 30,
                           prefix = "RNA_snn_res.", assay = "RNA", reduction = x,
                           saveRDS = F, return_obj = T, returnFeats = T,
                           features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                        "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                        "CD4", "MS4A1", "PPBP","HBM")
                          )
}

#generate viln plots using harmony clusters
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_sub_integrated.harmony", outName = outName,
          outDir = paste0("../output/viln/", outName, "/"), returnViln = T
         )

#identified more low quality cells - remove them - c10 (abnormal gex patterns), c17 (high pct mt enriched) and c20 (contaminating T cells) -- c15 may also need to be removed --
seu.obj <- subset(seu.obj, invert = T, subset = clusterID_sub_integrated.harmony %in% c(10, 17, 20))
seu.obj$clusterID_sub_filterStep2 <- seu.obj$clusterID_sub_integrated.harmony

#integrate the data using all of the four Seurat v5 integration methods
seu.obj <- integrateData(seu.obj = seu.obj, dout = "../output/s2/", outName = outName, runAllMethods = FALSE, 
                         method = "HarmonyIntegration", normalization.method = "LogNormalize", indReClus = T)

#use clustree to identify clustering parameters that appear most appropriate
clusTree(seu.obj = seu.obj, dout = "../output/clustree/", outName = outName, 
            test_dims = 30, algorithm = 3, prefix = "RNA_snn_res.", reduction = "integrated.harmony")

#complete data visualization
# for (x in list("integrated.cca", "integrated.harmony", "integrated.joint", "integrated.rcpa")) {
for (x in list("integrated.harmony")) {
    seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = paste0(outName, "_", x), 
                           final.dims = 30, final.res = 0.6, stashID = "clusterID_sub", algorithm = 3, min.dist = 0.2, n.neighbors = 20,
                           prefix = "RNA_snn_res.", assay = "RNA", reduction = x,
                           saveRDS = F, return_obj = T, returnFeats = T,
                           features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                        "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                        "CD4", "MS4A1", "PPBP","HBM")
                          )
}

saveRDS(seu.obj, paste0("../output/s3/", outName,"CLEAN_S3.rds"))


######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end preprocessing   ######## <<<<<<<<<<<<<<
######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



########################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin subset analysis   ######## <<<<<<<<<<<<<<
########################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in data and metadata
outName <- "mac"
seu.obj <- readRDS(paste0("../output/s3/", outName,"CLEAN_S3.rds"))
reduction <- "umap.integrated.harmony"
clusMain <- "clusterID_sub_integrated.harmony"
contrast <- c("Diseased", "Healthy")
seu.obj$majorID_sub <- droplevels(seu.obj$clusterID_integrated.harmony)
seu.obj$cellSource <- droplevels(seu.obj$cellSource)
seu.obj$name <- droplevels(seu.obj$name)
seu.obj$colz <- droplevels(seu.obj$colz)

# colArray <- read.csv("./metaData/.csv") #load colors if specified

#generate viln plots using harmony clusters
vilnPlots(seu.obj = seu.obj, groupBy = clusMain, outName = outName,
          outDir = paste0("../output/viln/", outName, "/"), returnViln = F
         )

### Export data for interactive cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = paste0(outName, "_CLEAN"), outDir = "../output/cb_input/", 
                markers = paste0("../output/viln/", outName, "/mac_clusterID_sub_integrated.harmony_gene_list.csv"),
                reduction = reduction,  
                colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", 
                                "majorID", clusMain, "name", "cellSource"), 
                skipEXPR = F, test = F,
                feats = c("CD24", "EBF1", "GPNMB", 
                          "RELB", "IL1A", "C1QC", "BOLA-DRA", 
                          "NIPA1", "ISG15", "OAS1Y", "CSF3R", 
                          "SDS", "MEFV", "CD68", "APOE")
              )    

#use singleR to ID cells
singleR(seu.obj = seu.obj, clusters = clusMain, reduction = reduction, 
        outDir = "../output/singleR/", outName = outName)


### Generate dot plots using vilnPlots resuts of majorID
pi <- autoDot(
    seu.integrated.obj = seu.obj, inFile = paste0("../output/viln/", outName, "/", outName, "_", clusMain, "_gene_list.csv"), 
    groupBy = clusMain, MIN_LOGFOLD_CHANGE = 0.5, MIN_PCT_CELLS_EXPR_GENE = 0.1, filterTerm = "ENSBTAG"
) + coord_flip()
ggsave(paste0("../output/", outName, "/", outName, "_autodot_major_clusID.png"), width = 6, height = 15)


### Plot initial cluster UMAP
pi <- DimPlot(seu.obj, 
                reduction = reduction, 
                group.by = clusMain,
#                 cols = colArray$newCol,
                pt.size = 0.25,
                label = TRUE,
                label.box = TRUE
) + NoLegend()
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = "black", smallAxes = F) 
ggsave(paste0("../output/", outName, "/", outName, "_rawUMAP.png"), width = 7, height = 7)


### Key feature plots
features <- c("PTPRC","CD3E","CTSW", 
                "CSF3R","S100A12", 
                "CD68","FLT3","FCER1A", 
                "GPNMB","VEGFB",
                "MS4A1","TOP2A")
p <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 3, title.size = 14, features = features, order = F, legJust = "top", reduction = reduction) 
ggsave(paste0("../output/", outName, "/", outName, "_featPlots.png"), width = 9, height = 15)


### Frequency plots to run stats - cluster
freqy <- freqPlots(seu.obj, method = 1, nrow = 2, 
                   comp = "cellSource", groupBy = "majorID_sub", legTitle = "Cell source", refVal = "name",
                   namez = "name", 
                   colz = "colz"
)
ggsave(paste0("../output/", outName, "/",outName, "_freqPlots_cluster.png"), width = 10, height = 4)



### Premutation testing - only works for 2 groups!!!

# Run scProportionTest::permutation_test
log2FD_threshold <- 1
seu.obj$name <- droplevels(seu.obj$name)
Idents(seu.obj) <- "name"
set.seed(12)
seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$name)))
prop_test <- sc_utils(seu.obj.sub)
prop_test <- permutation_test( prop_test, cluster_identity = "majorID_sub", sample_1 = "Healthy", sample_2 = "Diseased", sample_identity = "cellSource" )

p <- permutation_plot(prop_test)  + theme(axis.title.y = element_blank(),
                                          legend.position = "top") + guides(colour = guide_legend("", nrow = 2, 
                                                                                                  byrow = TRUE)) + coord_flip()
res.df <- prop_test@results$permutation
res.df <- res.df %>% mutate(Significance = as.factor(ifelse(obs_log2FD < -log2FD_threshold & FDR < 0.01,"Down",
                                                            ifelse(obs_log2FD > log2FD_threshold & FDR < 0.01,"Up","n.s.")))
                           ) %>% arrange(obs_log2FD)

res.df$clusters <- factor(res.df$clusters, levels = c(res.df$clusters))
p <- ggplot(res.df, aes(x = clusters, y = obs_log2FD)) + 
geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, 
                    color = Significance)) + theme_bw() + geom_hline(yintercept = log2FD_threshold, 
                                                                     lty = 2) + geom_hline(yintercept = -log2FD_threshold, 
                                                                                           lty = 2) + 
geom_hline(yintercept = 0) + scale_color_manual(values = c("blue", "red","grey")
                                               ) + theme(#axis.text.x = element_blank(),
                                                         axis.title.x = element_blank(),
                                                         legend.position = c(0.2, 0.9),
                                                         legend.title = element_blank(),
                                                         plot.margin = margin(c(7,7,21,7)),
                                                         legend.key = element_rect(fill = 'transparent', colour = NA),
                                                         legend.direction = 'vertical',
                                                         legend.background = element_rect(fill = 'transparent', colour = NA),
                                                         panel.background = element_rect(fill = 'transparent', colour = NA),
                                                         plot.background = element_rect(fill = "transparent", colour = NA)
                                                        ) + ylab(expression(atop(bold("abundance log2(FC)"),atop(italic("(Disease versus healthy)"))))) + guides(color = guide_legend(nrow=1))

ggsave(paste("../output/", outName, "/",outName, "_premTest.png", sep = ""), width = 3, height = 1, scale = 2 )


# Run custom permutation test on each sample
seu.obj$majorID_sub <- droplevels(seu.obj$majorID_sub)
seu.obj$orig.ident <- as.factor(seu.obj$orig.ident)
cluster <- "majorID_sub"
groupBy <- "cellSource"
pairBy <- "orig.ident"
log2FC_threshold <- 1
p.adj_threshold <- 0.01
nPerm <- 10000
set.seed(12)


sim.prop.res <- lapply(levels(seu.obj@meta.data[ , pairBy]), function(pairedName){

    message(paste("Running simmulation for", pairedName))
    
    #extract metadata
    meta.df <- as.data.frame(table(seu.obj@meta.data[ , pairBy], 
                                   seu.obj@meta.data[ , groupBy])) %>% filter(Freq != 0)
    condition <- meta.df %>% filter(Var1 == pairedName) %>% pull(Var2) %>% as.character()
    sampleInclude <-  meta.df %>% filter(Var2 != condition) %>% pull(Var1) %>% as.character()
    sampleInclude <- c(pairedName, sampleInclude)
    
    #remove extra samples and downsample to obtain equal representation
    Idents(seu.obj) <- pairBy
    seu.obj.sub <- subset(seu.obj, idents = sampleInclude)
    seu.obj.sub@meta.data[ , groupBy] <- droplevels(seu.obj.sub@meta.data[ , groupBy])
    
    cellsKeep <- lapply(as.character(meta.df[meta.df$Var2 != condition, ]$Var1), function(x){
        sample(rownames(seu.obj.sub@meta.data[seu.obj.sub$orig.ident == x, ]),
               size = min(meta.df[meta.df$Var2 != condition, ]$Freq), replace = F)        
    })
    cellsKeep <- c(unlist(cellsKeep), rownames(seu.obj.sub@meta.data[seu.obj.sub$orig.ident == pairedName, ]))
    seu.obj.sub <- subset(seu.obj, cells = cellsKeep)
    
    #extract additional metadata
    metaData <- seu.obj.sub@meta.data %>% select(match(c(cluster, groupBy), colnames(.)))
    message(paste("Log fold change will be calcualted by contrasting", levels(metaData[[groupBy]])[2],
                  "to", levels(metaData[[groupBy]])[1]))
    
    actual.prop.df <- as.data.frame(prop.table(table(metaData[ , cluster], droplevels(metaData[ , groupBy])), 2) * 100) %>% 
        pivot_wider(names_from = Var2, values_from = Freq) %>% 
        as.data.frame()
    colnames(actual.prop.df)[1] <- "cluster"
    actual.prop.df$log2FC <- log2(actual.prop.df[ ,3]) - log2(actual.prop.df[ ,2])
    clusToExclude <- actual.prop.df$cluster[is.na(actual.prop.df$log2FC)]

    #if no cells in either condition remove from analysis
    if (length(clusToExclude) > 0) {
        actual.prop.df <- na.omit(actual.prop.df)
        metaData <- metaData[! metaData[[cluster]] %in% clusToExclude, ]
        metaData[[cluster]] <- droplevels(metaData[[cluster]])
    }
    sim.prop <- matrix(NA, nrow(actual.prop.df), nPerm)
    for (i in 1:nPerm) {
        sim.data <- metaData
        sim.data[[groupBy]] <- sample(sim.data[[groupBy]])
        res <- prop.table(table(sim.data[ , cluster], droplevels(sim.data[ , groupBy])), 2) * 100
        sim.prop[ , i] <- log2(unname(res[ , 2])) - log2(unname(res[ , 1]))
    }

    sig.res <- matrix(NA, nrow(actual.prop.df), 1)
    Pvalue.list <- lapply(1:nrow(actual.prop.df), function(i){
        if (actual.prop.df$log2FC[i] > 0 ) {
            Pvalue <- (sum(sim.prop[i, ] >= actual.prop.df$log2FC[i])) / nPerm
        } else{
            Pvalue <- (sum(sim.prop[i, ] <= actual.prop.df$log2FC[i])) / nPerm
        }
        return(Pvalue)
    })
    actual.prop.df <- actual.prop.df %>% 
        mutate(
            "p.value" = unlist(Pvalue.list),
            "p.adj" = p.adjust(p.value),
            "significance" = as.factor(
                ifelse(log2FC < -log2FC_threshold & p.adj < p.adj_threshold, "Down",
                       ifelse(log2FC > log2FC_threshold & p.adj < p.adj_threshold, "Up", "n.s."))
            ),
            "replicate" = pairedName,
            "cohort" = condition,
            "lfc_contrast" = paste0(levels(metaData[[groupBy]])[2], "_vs_", levels(metaData[[groupBy]])[1])
        )
    return(actual.prop.df)
})

res.df <- do.call(rbind, sim.prop.res) 
res.df$cluster <- factor(res.df$cluster, levels = levels(seu.obj$majorID_sub))
res.df$replicate <- factor(res.df$replicate, levels = levels(seu.obj$orig.ident))
res.df$cohort <- factor(res.df$cohort, levels = levels(seu.obj$cellSource))
res.df$significance <- factor(res.df$significance, levels = c("Down", "Up", "n.s."))
res.df <- res.df[res.df$cohort == "Diseased",]

p <- ggplot(res.df, aes(x = cluster, y = log2FC, shape = significance)) + 
    theme_classic() + 
    geom_hline(yintercept = log2FC_threshold, lty = 2) + 
    geom_hline(yintercept = -log2FC_threshold, lty = 2) +
    geom_hline(yintercept = 0) + 
    geom_point(alpha = 0.8) + 
    scale_shape_manual(values = c(25, 24, 21))+
    theme(
        axis.title.x = element_blank(),
        legend.box = 'horizontal',
        legend.position = "top",
#         legend.title = element_blank(),
        plot.margin = margin(c(7,7,21,7)),
        legend.key = element_rect(fill = 'transparent', colour = NA),
        legend.direction = 'horizontal',
        legend.background = element_rect(fill = 'transparent', colour = NA),
        panel.background = element_rect(fill = 'transparent', colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    ) + 
    ylab(expression(atop(bold("abundance log2(FC)"), atop(italic("(Disease versus healthy)"))))) + 
    guides(
        shape = guide_legend(title = "Significance:", nrow = 1)
    ) +
    annotate("rect", xmin = 0, xmax = length(unique(res.df$cluster)) + 1, 
                      ymin = -log2FC_threshold, 
                      ymax = log2FC_threshold, 
                      alpha = 0.1, fill = "grey10")

ggsave(paste("../output/", outName, "/",outName, "_montePerm.png", sep = ""), width = 7, height = 3)


### Highlight c9 to see its spread over the umap
Idents(seu.obj) <- "majorID_sub"
pi <- DimPlot(seu.obj, 
              reduction = reduction, 
              cells.highlight = WhichCells(seu.obj, idents = "9"),
              cols.highlight = "#f22e9d",
              pt.size = 0.5,
              label = F,
              label.box = F,
              shuffle = F
) + ggtitle("c9 highlighted in pink")
ggsave(paste("../output/", outName, "/", "highlight_c9.png", sep = ""), width = 7, height = 7)

### Complete pseudobulk DGE by major cell types
createPB(seu.obj = seu.obj, groupBy = "majorID_sub", comp = "cellSource", biologicalRep = "name", lowFilter = T, dwnSam = F, 
         clusters = NULL, outDir = paste0("../output/", outName, "/pseudoBulk/")
)

pseudoDEG(metaPWD = paste0("../output/", outName, "/pseudoBulk/majorID_sub_deg_metaData.csv"),
          padj_cutoff = 0.05, lfcCut = 1, strict_lfc = F,
          outDir = paste0("../output/", outName, "/pseudoBulk/"), 
          outName = outName, 
          idents.1_NAME = contrast[1], idents.2_NAME = contrast[2],
          inDir = paste0("../output/", outName, "/pseudoBulk/"), title = "Myeloid", 
          filterTerm = "ZZZZ", addLabs = NULL, mkDir = T
)


### Series of plotting depicting sig DEGs
files <- list.files(path = paste0("../output/", outName, "/pseudoBulk/", levels(seu.obj$majorID_sub), "/"), 
           pattern = ".csv", all.files = FALSE, full.names = T)
df.list <- lapply(unlist(files), read.csv, header = T)
res.df <- do.call(rbind, df.list)


# Plot number of degs by cell type
cnts_mat <- res.df %>% 
    mutate(
        direction = ifelse(log2FoldChange > 0, "Up", "Down")
    ) %>% 
    group_by(gs_base,direction) %>% 
    summarize(nRow = n()) %>% 
    pivot_wider(names_from = gs_base, values_from = nRow) %>% 
    as.matrix() %>% t()

colnames(cnts_mat) <- cnts_mat[1,]
cnts_mat <- cnts_mat[-c(1),]
class(cnts_mat) <- "numeric"
cnts_mat[is.na(cnts_mat)] <- 0

#order by number of total # of DEGs
orderList <- rev(rownames(cnts_mat)[order(rowSums(cnts_mat))])
cnts_mat <- cnts_mat[match(orderList, rownames(cnts_mat)),]        

png(file = paste0("../output/", outName, "/", outName, "_deg_heat.png"), width=1000, height=2000, res=400)
par(mfcol=c(1,1))         
ht <- Heatmap(cnts_mat,#name = "mat", #col = col_fun,
              name = "# of DEGs",
              cluster_rows = F,
              row_title = "Cell type",
              col = circlize ::colorRamp2(c(0,max(cnts_mat)), colors = c("white","red")),
              cluster_columns = F,
              column_title = "# of DEGs",
              show_column_names = TRUE,
              column_title_side = "top",
              column_names_rot = 0,
              column_names_centered = TRUE,
              heatmap_legend_param = list(legend_direction = "horizontal", title_position = "topleft",  title_gp = gpar(fontsize = 16), 
                                          labels_gp = gpar(fontsize = 8), legend_width = unit(6, "cm")),
              cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(sprintf("%.0f", as.matrix(cnts_mat)[i, j]), x, y, gp = gpar(fontsize = 14, col = "black"))
              })
draw(ht, padding = unit(c(2, 12, 2, 5), "mm"),show_heatmap_legend = FALSE)
dev.off()

# Plot top 10 degs of each cell type (very crowded plating) 
res.df <- res.df %>% 
    filter(!grepl("^ENS", gene)) %>%
    group_by(gs_base) %>%
    top_n(n = -10, wt = padj)

cond_colz <- c("mediumseagreen","mediumpurple1")
names(cond_colz) <- c("Healthy","Diseased")

res.df2 <- res.df[!duplicated(res.df$gene), ] %>% filter(!grepl("^ENS", gene)) %>% top_n(n = -10, wt = padj)
p <- splitDot(
    seu.obj = seu.obj,
    groupBy = "majorID_sub",
    splitBy = "cellSource",
    namedColz = cond_colz,
    left_margin = 35,
    geneList_UP = filter(res.df2, log2FoldChange > 0) %>% pull(gene),
    geneList_DWN = filter(res.df2, log2FoldChange < 0) %>% pull(gene),
    geneColz = c("red", "blue")
    ) + theme(legend.direction = "horizontal", legend.position = "top", legend.box = "horizontal")
ggsave(paste0("../output/", outName, "/", "splitDot.png"), width = 16, height = 8)

# Plot heatmap with significance indicator
clus <- gg_color_hue(length(levels(seu.obj$majorID_sub)))
names(clus) <- paste0("(c", ((1:length(levels(seu.obj$majorID_sub))) - 1), ") ", levels(seu.obj$majorID_sub))
clus_colz <- clus

res.df$gs_base <- toupper(res.df$gs_base)
ht <- sigDEG_heatmap(
    seu.obj = seu.obj, groupBy = "majorID_sub", splitBy = "cellSource", forceCleanPlot = F, 
    dge_res = res.df, lfc_thres = 1, cond_colz = cond_colz, clus_colz = clus,
    saveName = paste0("../output/", outName, "/", "splitHeat.png"),
    ht_height = 4500, ht_width = 3000
)


### Complete GSEA 
lapply(levels(seu.obj$majorID_sub), function(group){

    df <- res.df %>% filter(gs_base == group) %>% arrange(padj)

    upGenes <- df %>% filter(log2FoldChange > 0) %>% pull(gene)
    dwnGenes <- df %>% filter(log2FoldChange < 0) %>% pull(gene)
        skip <- FALSE
        if (length(upGenes) > 5 & length(dwnGenes) > 5) {
            p <- plotGSEA(geneList = upGenes, geneListDwn = dwnGenes, category = "C5",
                          subcategory = "GO:BP", upCol = "red", dwnCol = "blue", size = 3.5,
                          saveRes = paste0("../output/", outName, "/gsea/", group, "_gseaRes.csv")
                         )
        } else if (length(upGenes) > 5 & length(dwnGenes) <= 5) {
            p <- plotGSEA(geneList = upGenes, upOnly = T, category = "C5",
                          subcategory = "GO:BP", upCol = "red", dwnCol = "blue", size = 3.5,
                          saveRes = paste0("../output/", outName, "/gsea/", group, "_gseaRes.csv")
                         )
        } else if (length(upGenes) <= 5 & length(dwnGenes) > 5) {
            p <- plotGSEA(geneListDwn = dwnGenes, dwnOnly = T, category = "C5",
                          subcategory = "GO:BP", upCol = "red", dwnCol = "blue", size = 3.5,
                          saveRes = paste0("../output/", outName, "/gsea/", group, "_gseaRes.csv")
                         )
        } else {
            skip <- TRUE
        }
        if (nrow(p$data) > 0 & ! skip) {
            minVal <- -15
            maxVal <- 15
            pi <- p + scale_x_continuous(limits = c(minVal, maxVal), name = "Signed log10(padj)")
            ggsave(paste("../output/", outName, "/gsea/", group, "_gsea.png", sep = ""), width = 10, height = 5)
        } else{
            message(paste("No pathways enriched for", group, "comparison."))
        }
})


#organize gsea results
df.list <- lapply(c("gsea"), function(timepoint){
    res.df.list <- lapply(levels(seu.obj$majorID_sub), function(group){
        df <- try(read.csv(paste0("../output/", outName, "/", timepoint, "/", group, "_gseaRes.csv"), header = T), silent = T)
        if (class(df) != "try-error") {
            if (nrow(df) > 0) {
                df$CellType <- group
                df$TimePoint <- timepoint
                message(paste("Saving", group, "and", timepoint))
                makePseudoDf <- FALSE
            } else { makePseudoDf <- TRUE }
        } else { makePseudoDf <- TRUE }
        
        if (makePseudoDf) {
            df <- t(as.data.frame(
                setNames(
                    c(rep("NO_TERMS_ENRICHED", 10), 0, "NONE", group, timepoint),
                    c("X", "ID", "Description", "GeneRatio", "BgRatio", "pvalue", 
                      "p.adjust", "qvalue", "geneID", "Count", "x_axis", "direction", 
                      "CellType", "TimePoint")
                )
            ))
            rownames(df) <- 1
                
#                 "CellType" = group,
#                 "TimePoint" = timepoint,
#                 "ID" = "NO_TERMS_ENRICHED",
#                 "x_axis" = 0
        }
        return(df)
    })
    res.df <- do.call(rbind, res.df.list)
    return(res.df)
})
res.df <- do.call(rbind, df.list)

#clean up the gsea terms
res.df$ID <- gsub("GOBP_", "", res.df$ID)
res.df$ID <- gsub("ANTIGEN_PROCESSING_AND_PRESENTATION_", "ANTIGEN_PRESENTATION_", res.df$ID)

#data cleaning and conversion to matrix
df <- res.df %>% 
    mutate(
        row = row_number(),
        colName = paste0(CellType, "-_-", TimePoint)
    ) %>% 
    select(colName, ID, x_axis, row) %>%
    pivot_wider(names_from = colName, values_from = x_axis) %>% 
    as.data.frame()

mat <- aggregate(df[ , 3:ncol(df)], by = df[1], FUN = function(x){sum(as.numeric(x), na.rm = TRUE)}) %>%
    column_to_rownames("ID") %>%
    as.matrix()

mat <- mat[rownames(mat) != "NO_TERMS_ENRICHED", ]

#set annotations
cond_colz <- setNames(
    c("#C47AEA"),
    c("Disease_VS_Healthy")
)
clus_colz <- setNames(
    gg_color_hue(length(unique(unlist(lapply(colnames(mat), function(x){strsplit(x, "-_-")[[1]][1]}))))),
    unique(unlist(lapply(colnames(mat), function(x){strsplit(x, "-_-")[[1]][1]})))
)
ha <- HeatmapAnnotation(
    Cluster = unlist(lapply(colnames(mat), function(x){strsplit(x, "-_-")[[1]][1]})),
    border = TRUE,
    col = list(Cluster = clus_colz),
    annotation_legend_param = list(
        Cluster = list(direction = "horizontal", nrow = 1)
    ),
    show_annotation_name = FALSE,
    show_legend = FALSE
)
lgd1 <- Legend(labels = unique(unlist(lapply(colnames(mat), function(x){strsplit(x, "-_-")[[1]][1]}))),
               legend_gp = gpar(fill = clus_colz), 
               title = "Cluster", 
               direction = "vertical",
               nrow = 3, 
               gap = unit(0.6, "cm")
              )
lgd2 <- Legend(labels = names(cond_colz),
               legend_gp = gpar(fill = cond_colz), 
               title = "Condition", 
               direction = "vertical",
               nrow = 2
              )
pd <- packLegend(lgd1, lgd2, max_width = unit(45, "cm"), 
    direction = "horizontal", column_gap = unit(5, "mm"), row_gap = unit(0.5, "cm"))

#complete hc and extract most sig genes to label
M <- (1- cor(t(mat), method = "pearson"))/2
hc_row <- hclust(as.dist(M), method = "complete")

# M <- (1- cor(mat, method = "pearson"))/2
# hc_col <- hclust(as.dist(M), method = "complete")

topn <- rownames(mat)[rev(order(rowSums(abs(mat))))][1:50]
topn <- topn[nchar(topn) <= 75]
position <- match(topn, hc_row$labels[hc_row$order])
ra_right <- rowAnnotation(key_feats = anno_mark(at = hc_row$order[position], 
                                                labels = hc_row$labels[hc_row$order][position],
                                                labels_gp = gpar(fontsize = 6)))

#create the heatmap
ht <- Heatmap(
    mat,
    name = "Signed log10(FDR)",
    col = circlize::colorRamp2(c(-4, 0, 4), c("blue", "#EEEEEE", "red")),
    cluster_rows = hc_row,
    row_title_gp = gpar(fontsize = 24),
    show_row_names = F,
    cluster_columns = F,
    show_column_names = F,
    top_annotation = ha, 
    column_names_gp = gpar(fontsize = 12),
    right_annotation = ra_right,
    column_split = unlist(lapply(colnames(mat), function(x){strsplit(x, "-_-")[[1]][2]})),
    heatmap_legend_param = list(
        direction = "horizontal",
        legend_width = unit(6, "cm")
    )
)
png(file = paste0("../output/", outName, "/",  outName, "_heatTest.png"), width = 3000, height = 4000, res = 400)
par(mfcol = c(1, 1))   
draw(
    ht, padding = unit(c(2, 10, 2, 2), "mm"), heatmap_legend_side = "bottom",
    annotation_legend_list = pd, annotation_legend_side = "top"
)
dev.off()
