#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")

######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   customize variables   ######## <<<<<<<<<<<<<<
######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#set output param
outName <- "tcell"
reduction <- "umap.integrated.harmony"
reduction2 <- "integrated.harmony" 
clusMain <- "clusterID_integrated.harmony"
contrast <- c("Diseased", "Healthy") #code will test first vs second

#load in the proprocessed data & subset on pop of interest
seu.obj <- readRDS("../output/s3/20240313_bov_lav_n5n5_dxVSh_S3.rds")
seu.obj <- cleanMeta(seu.obj)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID.csv", groupBy = "clusterID_integrated.harmony", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident", metaAdd = "cellSource")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident", metaAdd = "name")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj <- subset(seu.obj, subset = majorID == "tcell")


######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin preprocessing   ######## <<<<<<<<<<<<<<
######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#integrate the data using all of the four Seurat v5 integration methods
seu.obj <- integrateData(seu.obj = seu.obj, dout = "../output/s2/", outName = outName, 
                         runAllMethods = TRUE, normalization.method = "LogNormalize", indReClus = T,
                         min.cell = 100) # this exludes bov_lav_453212 and bov_lav_46425 from analysis

#use clustree to identify clustering parameters that appear most appropriate
clusTree(seu.obj = seu.obj, dout = "../output/clustree/", outName = outName, 
            test_dims = 30, algorithm = 3, prefix = "RNA_snn_res.", reduction = reduction2)

#complete data visualization
for (x in list("integrated.cca", "integrated.harmony", "integrated.joint", "integrated.rcpa")) {
    seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = paste0(outName, "_", x), 
                           final.dims = 30, final.res = 0.6, stashID = "clusterID", algorithm = 3, min.dist = 0.3, n.neighbors = 30,
                           prefix = "RNA_snn_res.", assay = "RNA", reduction = x,
                           saveRDS = F, return_obj = T, returnFeats = T,
                           features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                        "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                        "CD4", "MS4A1", "PPBP","HBM")
                          )
}
saveRDS(seu.obj, paste0("../output/s3/", outName,"_S3.rds"))

######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end preprocessing   ######## <<<<<<<<<<<<<<
######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



########################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin subset analysis   ######## <<<<<<<<<<<<<<
########################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in data and metadata
outName <- "tcell"
reduction <- "umap.integrated.harmony"
clusMain <- "clusterID_integrated.harmony"
contrast <- c("Diseased", "Healthy") #code will test first vs second
seu.obj <- readRDS(paste0("../output/s3/", outName,"_S3.rds"))
# colArray <- read.csv("./metaData/.csv") #load colors if specified

#generate viln plots using harmony clusters
vilnPlots(seu.obj = seu.obj, groupBy = clusMain, outName = outName,
          outDir = paste0("../output/viln/", outName, "/"), returnViln = T 
         )

### Export data for interactive cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, outDir = "../output/cb_input/", 
                markers = paste0("../output/viln/", outName, "/", outName, "_", clusMain,"_gene_list.csv"),
                reduction = reduction,  
                colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", 
                                "majorID", clusMain, "name", "cellSource"), 
                skipEXPR = F, test = F,
                feats = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                          "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                          "CD4", "MS4A1", "PPBP", "HBM")
)    

#use singleR to ID cells
singleR(seu.obj = seu.obj, clusters = clusMain, reduction = reduction, 
        outDir = "../output/singleR/", outName = outName)


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


### UMAP by sample -- if unequal sample size downsample by cellSource
Idents(seu.obj) <- "orig.ident"
set.seed(12)
seu.obj.ds <- subset(x = seu.obj, downsample = min(table(seu.obj$orig.ident)))
pi <- DimPlot(seu.obj.ds, 
                reduction = reduction, 
                group.by = "name",
                cols = levels(seu.obj.ds$colz),
                pt.size = 0.25,
                label = FALSE,
                shuffle = TRUE
)
p <- formatUMAP(pi) + labs(colour="Cell source:") + theme(legend.position = "top", 
                                                            legend.direction = "horizontal",
                                                            legend.title=element_text(size=12)
                                                            ) + guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4)))
ggsave(paste0("../output/", outName, "/", outName, "_umap_bySample.png"), width =7, height = 7)


### Frequency plots to run stats - cluster
freqy <- freqPlots(seu.obj, method = 1, nrow = 3, 
                   comp = "cellSource", groupBy = clusMain, legTitle = "Cell source", refVal = "name",
                   namez = "name", 
                   colz = "colz"
)
ggsave(paste0("../output/", outName, "/",outName, "_freqPlots_cluster.png"), width = 12, height = 8)


### Stacked bar graph by clusMain
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "name", groupBy = "name", clusters = clusMain) + scale_x_discrete(expand = c(0, 0)) +
scale_fill_manual(labels = levels(seu.obj$name), 
                    values = levels(seu.obj$colz)) + 
    theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12))
ggsave(paste0("../output/", outName, "/", outName, "_stackedBar_cluster.png"), width =7, height = 5)


### Complete pseudobulk DGE by all cells
createPB(seu.obj = seu.obj, groupBy = "allCells", comp = "cellSource", biologicalRep = "name", lowFilter = T, dwnSam = F, 
         clusters = NULL, outDir = paste0("../output/", outName, "/pseudoBulk/")
)

pseudoDEG(metaPWD = paste0("../output/", outName, "/pseudoBulk/allCells_deg_metaData.csv"),
          padj_cutoff = 0.05, lfcCut = 0.58, outDir = paste0("../output/", outName, "/pseudoBulk/"), 
          outName = outName, 
          idents.1_NAME = contrast[1], idents.2_NAME = contrast[2],
          inDir = paste0("../output/", outName, "/pseudoBulk/"), title = "T cells", 
          filterTerm = "ZZZZ", addLabs = NULL, mkDir = T
)


### Fig - Split UMAP of selected DEGs
Idents(seu.obj) <- "cellSource"
seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$cellSource)))

df <- read.csv(paste0("../output/", outName, "/pseudoBulk/allCells/", outName, "_cluster_allCells_all_genes.csv")) %>% arrange(padj)
upGenes <- df %>% filter(log2FoldChange > 0) %>% pull(gene)

features <- head(upGenes, 10)

p <- FeaturePlot(seu.obj.sub,features = features, pt.size = 0.1, reduction = reduction, split.by = "cellSource", order = T, cols = c("lightgrey","darkblue"), by.col = F) + labs(x = "UMAP1", y = "UMAP2") & theme(axis.text = element_blank(),
                                                                                                                                                axis.title.y.right = element_text(size = 11),
                                                                                                                                            
                                                                                                                           axis.ticks = element_blank(),
                                                                                                                           axis.title = element_blank(),
                                                                                                                           axis.line = element_blank(),
                                                                                                                           plot.title = element_text(size=11),
                                                                                                                           title = element_blank(),
                                                                                                                           plot.margin = unit(c(0, 0, 0, 0), "cm")
                                                                                                                          ) 
ggsave(paste("../output/", outName, "/", "splitFeats.png", sep = ""), width = 12, height = 4)





### Complete pseudobulk DGE by all cells
createPB(seu.obj = seu.obj, groupBy = clusMain, comp = "cellSource", biologicalRep = "name", lowFilter = T, dwnSam = F, 
         clusters = NULL, outDir = paste0("../output/", outName, "/pseudoBulk/")
)

pseudoDEG(metaPWD = paste0("../output/", outName, "/pseudoBulk/clusterID_integrated.harmony_deg_metaData.csv"),
          padj_cutoff = 0.05, lfcCut = 0.58, outDir = paste0("../output/", outName, "/pseudoBulk/"), 
          outName = outName, 
          idents.1_NAME = contrast[1], idents.2_NAME = contrast[2],
          inDir = paste0("../output/", outName, "/pseudoBulk/"), title = "T cells", 
          filterTerm = "ZZZZ", addLabs = NULL, mkDir = T
)

### Supp fig xx -- heatmap of sig DEGs
files <- list.files(path = paste0("../output/", outName, "/pseudoBulk/", levels(seu.obj$clusterID_integrated.harmony), "/"), 
           pattern = ".csv", all.files = FALSE, full.names = T)
df.list <- lapply(unlist(files), read.csv, header = T)
res.df <- do.call(rbind, df.list)

res.df <- res.df %>% 
    filter(!grepl("^ENSCAFG", gene)) %>%
    group_by(gs_base) %>%
    top_n(n = -20, wt = padj)

cond_colz <- c("mediumseagreen","mediumpurple1")
names(cond_colz) <- c("Healthy","Diseased")

res.df2 <- res.df[!duplicated(res.df$gene), ] %>% filter(!grepl("^ENSCAFG", gene)) %>% top_n(n = -20, wt = padj)
p <- splitDot(
    seu.obj = seu.obj,
    groupBy = "clusterID_integrated.harmony",
    splitBy = "cellSource",
    namedColz = cond_colz,
    geneList_UP = filter(res.df2, log2FoldChange > 0) %>% pull(gene),
    geneList_DWN = filter(res.df2, log2FoldChange < 0) %>% pull(gene),
    geneColz = c("red", "blue")
    )
ggsave(paste0("../output/", outName, "/", "splitDot.png"), width = 14, height = 7)
    

clus <- gg_color_hue(length(levels(seu.obj$clusterID_integrated.harmony)))
names(clus) <- paste0("(c", ((1:length(levels(seu.obj$clusterID_integrated.harmony))) - 1), ") ", levels(seu.obj$clusterID_integrated.harmony))
clus_colz <- clus

res.df$gs_base <- toupper(res.df$gs_base)
ht <- sigDEG_heatmap(
    seu.obj = seu.obj, groupBy = "clusterID_integrated.harmony", splitBy = "cellSource", forceCleanPlot = F, 
    dge_res = res.df, lfc_thres = 1, cond_colz = cond_colz, clus_colz = clus,
    saveName = paste0("../output/", outName, "/", "splitHeat.png"),
    ht_height = 4500, ht_width = 3000
)


### Complete GSEA 
lapply(levels(seu.obj$clusterID_integrated.harmony), function(group){

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
                          saveRes = ppaste0("../output/", outName, "/gsea/", group, "_gseaRes.csv")
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
    res.df.list <- lapply(levels(seu.obj$clusterID_integrated.harmony), function(group){
        df <- try(read.csv(paste0("../output/", outName, "/", timepoint, "/", group, "_gseaRes.csv"), header = T))
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