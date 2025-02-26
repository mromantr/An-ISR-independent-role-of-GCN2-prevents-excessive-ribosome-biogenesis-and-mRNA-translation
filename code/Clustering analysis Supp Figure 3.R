
#########################################
# Log2FC heatmap with previous study data
#########################################

  require(DESeq2)
  require(ComplexHeatmap)
  require(circlize)
  require(fgsea)
  require(tidyverse)
  require(msigdbr)

  options(scipen = 99)

  # read in previous data
    studyid <- 'Cell lines RNAseq'

    dds <- readRDS(paste0('../../HolgerImperial/18/', studyid, '/12.dds.Rds'))
    idx <- which(colData(dds)$Treatment == 'GiB')
    dds <- dds[,idx]
    colData(dds)$Dependency <- ifelse(colData(dds)$Line %in% 'A375', 'MostDependent',
      ifelse(colData(dds)$Line %in% c('A375MA2', 'MPCa', 'A549'), 'DependentSomewhat',
        ifelse(colData(dds)$Line %in% c('IPC', 'HepG2'), 'Independent', NA)))
    dds <- DESeqDataSetFromMatrix(counts(dds, normalized = FALSE), colData = colData(dds), design = ~ Dependency)
    dds <- DESeq(dds)
    rld <- rlog(dds, blind = FALSE)
    vst <- vst(dds, blind = FALSE)

    colData(dds)$Line
    a375 <- data.frame(readxl::read_xlsx('../../HolgerImperial/18/Cell lines RNAseq/A375 GiB48h_vs_DMSO48h.xlsx')[,c('symbol','log2FoldChange')])
    a375ma2 <- data.frame(readxl::read_xlsx('../../HolgerImperial/18/Cell lines RNAseq/A375MA2 GiB48h_vs_DMSO48h.xlsx')[,c('symbol','log2FoldChange')])
    a549 <- data.frame(readxl::read_xlsx('../../HolgerImperial/18/Cell lines RNAseq/A549 GiB48h_vs_DMSO48h.xlsx')[,c('symbol','log2FoldChange')])
    hepg2 <- data.frame(readxl::read_xlsx('../../HolgerImperial/18/Cell lines RNAseq/HepG2 GiB48h_vs_DMSO48h.xlsx')[,c('symbol','log2FoldChange')])
    ipc <- data.frame(readxl::read_xlsx('../../HolgerImperial/18/Cell lines RNAseq/IPC298 GiB48h_vs_DMSO48h.xlsx')[,c('symbol','log2FoldChange')])
    mpca <- data.frame(readxl::read_xlsx('../../HolgerImperial/18/Cell lines RNAseq/MiaPaCa2 GiB48h_vs_DMSO48h.xlsx')[,c('symbol','log2FoldChange')])
    symbols <- sort(unique(c(a375$symbol, a375ma2$symbol, a549$symbol, hepg2$symbol, ipc$symbol, mpca$symbol)))
    l2fc <- data.frame(
      A375 = a375[match(symbols, a375$symbol),'log2FoldChange'],
      A375MA2 = a375ma2[match(symbols, a375ma2$symbol),'log2FoldChange'],
      A549 = a549[match(symbols, a549$symbol),'log2FoldChange'],
      HepG2 = hepg2[match(symbols, hepg2$symbol),'log2FoldChange'],
      IPC = ipc[match(symbols, ipc$symbol),'log2FoldChange'],
      MPCa = mpca[match(symbols, mpca$symbol),'log2FoldChange'],
      row.names = symbols)
    l2fcmeta <- colData(dds)[!duplicated(colData(dds)$Line),][,c('Line','Dependency')]
    l2fcmeta <- l2fcmeta[match(colnames(l2fc), l2fcmeta$Line),]

  # now add in l2fc data for the current GiB vs DMSO at 48h
    studyid <- 'GCN2 Revisions'

    setwd('out/deseq2 Less Outlier/protein_coding/')
    load(paste0(studyid, ' savepoint2.rdata'))

    huvec <- read.table('Exports/GCN2 Revisions DEA GiB vs DM 48h 2025-01-04.tsv', header = TRUE, sep = '\t')
    l2fc <- l2fc[which(rownames(l2fc) %in% huvec$symbol),]
    huvec <- huvec[which(huvec$symbol %in% rownames(l2fc)),]
    huvec <- huvec[match(rownames(l2fc), huvec$symbol),]
    all(rownames(l2fc) == huvec$symbol)

    l2fc$HUVEC <- huvec$log2FoldChange
    l2fcmeta <- data.frame(
      Line = c(as.character(l2fcmeta$Line), 'HUVEC'),
      Dependency = c(as.character(l2fcmeta$Dependency), 'HUVEC'),
      row.names = c(rownames(l2fcmeta), 'HUVEC'))

  # all DE genes
    MostDependent_vs_Independent <- read.table(
      '../../../../../HolgerImperial/18/Cell lines RNAseq/out/DEA MostDependent vs Independent 2024-05-11.tsv', header = TRUE)
    MostDependent_vs_DependentSomewhat <- read.table(
      '../../../../../HolgerImperial/18/Cell lines RNAseq/out/DEA MostDependent vs DependentSomewhat 2024-05-11.tsv', header = TRUE)
    DependentSomewhat_vs_Independent <- read.table(
      '../../../../../HolgerImperial/18/Cell lines RNAseq/out/DEA DependentSomewhat vs Independent 2024-05-11.tsv', header = TRUE)

    myCol <- colorRampPalette(c('blue2', 'white', 'red2'))(100)
    myBreaks <- seq(-1, 1, length.out = 100)
    log2cutoff <- 2
    qvaluecutoff <- 0.01
    sigGenes <- unique(c(
      subset(MostDependent_vs_Independent, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff)$symbol,
      subset(MostDependent_vs_DependentSomewhat, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff)$symbol,
      subset(DependentSomewhat_vs_Independent, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff)$symbol,
      subset(huvec, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff)$symbol))
    length(sigGenes)
    heat <- l2fc[which(rownames(l2fc) %in% sigGenes),]
    heat <- data.matrix(heat)
    ann <- data.frame(Dependency = l2fcmeta$Dependency)
    colours <- list('Dependency' = c('MostDependent' = 'red', 'DependentSomewhat' = 'forestgreen',
      'Independent' = 'royalblue', 'HUVEC' = 'black'))
    colAnn <- HeatmapAnnotation(df = ann, which = 'col', col = colours,
      annotation_width = unit(2, 'cm'), gap = unit(2, 'mm'),
      annotation_legend_param = list(Dependency = list(nrow = 4, title = '',
        title_position = 'leftcenter', legend_direction = 'vertical',
        title_gp = gpar(fontsize = 12, fontface = 'bold'),
        labels_gp = gpar(fontsize = 12, fontface = 'bold'))))
    boxplotCol <- HeatmapAnnotation(
      boxplot = anno_boxplot(
        data.matrix(heat),
        border = FALSE,
        gp = gpar(fill = '#CCCCCC'),
        pch = '.',
        size = unit(2, 'mm'),
        axis = TRUE,
        axis_param = list(gp = gpar(fontsize = 12), side = 'left')),
      annotation_width = unit(c(2.0), 'cm'), which = 'col')
    boxplotRow <- HeatmapAnnotation(
      boxplot = row_anno_boxplot(
        data.matrix(heat),
        border = FALSE,
        gp = gpar(fill = '#CCCCCC'),
        pch = '.',
        size = unit(2, 'mm'),
        axis = TRUE,
        axis_param = list(gp = gpar(fontsize = 12), side = 'top')),
      annotation_width = unit(c(2.0), 'cm'), which = 'row')
    genelabels <- rowAnnotation(
      Genes = anno_mark(
        at = seq(1, nrow(heat), 70),
        labels = rownames(heat)[seq(1, nrow(heat), 70)],
        labels_gp = gpar(fontsize = 14, fontface = 'bold'),
        padding = 0.5,
        link_width = unit(15, 'mm')),
      width = unit(2.0, 'cm') +
      max_text_width(
        rownames(heat)[seq(1, nrow(heat), 70)],
        gp = gpar(fontsize = 14,  fontface = 'bold')))
    png(paste0('Plots/ComplexHeatmap Supervised Log2FC ', Sys.Date(), '.png'), units = 'in', res = 300,
      width = 6.5, height = 15)
      hmap1 <- Heatmap(heat, name = 'Log [base 2] fold-change', col = colorRamp2(myBreaks, myCol), row_km = 1,
        #column_split = colData(dds)$Dependency,
        heatmap_legend_param = list(color_bar = 'continuous', legend_direction = 'horizontal',
          legend_width = unit(8,'cm'), legend_height = unit(5.5,'cm'), title_position = 'topcenter',
            title_gp = gpar(fontsize = 18, fontface = 'bold'), labels_gp = gpar(fontsize = 16, fontface = 'bold')),
        cluster_rows = TRUE, show_row_dend = FALSE, row_title = 'Statistically significant genes',
        row_title_side = 'left', row_title_gp = gpar(fontsize = 18,  fontface = 'bold'),
        row_title_rot = 90, show_row_names = FALSE, row_names_gp = gpar(fontsize = 16, fontface = 'plain'),
        row_names_side = 'left', row_dend_width = unit(25,'mm'),
        cluster_columns = TRUE, show_column_dend = TRUE, column_title = '', column_title_side = 'bottom',
        column_title_gp = gpar(fontsize = 18, fontface = 'bold'), column_title_rot = 0,
        show_column_names = TRUE, column_names_gp = gpar(fontsize = 16, fontface = 'plain'),
        #column_names_max_height = unit(16, 'cm'),
        column_dend_height = unit(25,'mm'),
        #clustering_distance_columns = function(x) as.dist(1-cor(t(x))),
        clustering_method_columns = 'ward.D2',
        #clustering_distance_rows = function(x) as.dist(1-cor(t(x))),
        clustering_method_rows = 'ward.D2',
        top_annotation = colAnn,
        bottom_annotation = boxplotCol)
      pushViewport(viewport(layout = grid.layout(nr = 1, nc = 1)))
        pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
          draw(hmap1 + genelabels, heatmap_legend_side = 'bottom',
            annotation_legend_side = 'top', row_sub_title_side = 'left',
            newpage = FALSE)
        popViewport()
      popViewport()
    dev.off()

  # signature / pathway genes
    plotheatmap <- function(name, width, height, genes, labelskip) {
      genes <- genes[!duplicated(genes)]
      myCol <- colorRampPalette(c('blue2', 'white', 'red2'))(100)
      myBreaks <- seq(-1, 1, length.out = 100)
      heat <- l2fc[which(rownames(l2fc) %in% genes),]
      heat <- data.matrix(heat)
      ann <- data.frame(Dependency = l2fcmeta$Dependency)
      colours <- list('Dependency' = c('MostDependent' = 'red', 'DependentSomewhat' = 'forestgreen',
        'Independent' = 'royalblue', 'HUVEC' = 'black'))
      colAnn <- HeatmapAnnotation(df = ann, which = 'col', col = colours,
        annotation_width = unit(2, 'cm'), gap = unit(2, 'mm'),
        annotation_legend_param = list(Dependency = list(nrow = 4, title = '',
          title_position = 'leftcenter', legend_direction = 'vertical',
          title_gp = gpar(fontsize = 12, fontface = 'bold'),
          labels_gp = gpar(fontsize = 12, fontface = 'bold'))))
      boxplotCol <- HeatmapAnnotation(
        boxplot = anno_boxplot(
          data.matrix(heat),
          border = FALSE,
          gp = gpar(fill = '#CCCCCC'),
          pch = '.',
          size = unit(2, 'mm'),
          axis = TRUE,
          axis_param = list(gp = gpar(fontsize = 12), side = 'right')),
        annotation_width = unit(c(2.0), 'cm'), which = 'col')
      genelabels <- rowAnnotation(
        Genes = anno_mark(
          at = seq(1, nrow(heat), labelskip),
          labels = rownames(heat)[seq(1, nrow(heat), labelskip)],
          labels_gp = gpar(fontsize = 14, fontface = 'bold'),
          padding = 0.5,
          link_width = unit(15, 'mm')),
        width = unit(2.0, 'cm') +
        max_text_width(
          rownames(heat)[seq(1, nrow(heat), labelskip)],
          gp = gpar(fontsize = 14,  fontface = 'bold')))
      png(paste0('Plots/ComplexHeatmap Supervised Log2FC ', name, ' ', Sys.Date(), '.png'), units = 'in', res = 300,
        width = width, height = height)
        hmap1 <- Heatmap(heat, name = 'Log [base 2] fold-change', col = colorRamp2(myBreaks, myCol), row_km = 1,
          #column_split = colData(dds)$Dependency,
          heatmap_legend_param = list(color_bar = 'continuous', legend_direction = 'horizontal',
            legend_width = unit(8,'cm'), legend_height = unit(5.5,'cm'), title_position = 'topcenter',
              title_gp = gpar(fontsize = 18, fontface = 'bold'), labels_gp = gpar(fontsize = 16, fontface = 'bold')),
          cluster_rows = TRUE, show_row_dend = FALSE, row_title = 'Pathway / Signature Genes',
          row_title_side = 'left', row_title_gp = gpar(fontsize = 18,  fontface = 'bold'),
          row_title_rot = 90, show_row_names = FALSE, row_names_gp = gpar(fontsize = 16, fontface = 'plain'),
          row_names_side = 'left', row_dend_width = unit(25,'mm'),
          cluster_columns = TRUE, show_column_dend = TRUE, column_title = '', column_title_side = 'bottom',
          column_title_gp = gpar(fontsize = 18, fontface = 'bold'), column_title_rot = 0,
          show_column_names = TRUE, column_names_gp = gpar(fontsize = 16, fontface = 'plain'),
          #column_names_max_height = unit(16, 'cm'),
          column_dend_height = unit(25,'mm'),
          #clustering_distance_columns = function(x) as.dist(1-cor(t(x))),
          clustering_method_columns = 'ward.D2',
          #clustering_distance_rows = function(x) as.dist(1-cor(t(x))),
          clustering_method_rows = 'ward.D2',
          top_annotation = colAnn,
          bottom_annotation = boxplotCol)
        pushViewport(viewport(layout = grid.layout(nr = 1, nc = 1)))
          pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
            draw(hmap1 + genelabels, heatmap_legend_side = 'bottom',
              annotation_legend_side = 'top', row_sub_title_side = 'left',
              newpage = FALSE)
          popViewport()
        popViewport()
      dev.off()
    }

    genes <- data.frame(readxl::read_xlsx('../../../../../HolgerImperial/18/Ribosomal proteins and translation factors.xlsx',
      sheet = 'Cytoplasm Ribosomal Proteins'))[,1]
    plotheatmap('Cytoplasm Ribosomal Proteins', width = 6.5, height = 12, genes, labelskip = 3)

    genes <- data.frame(readxl::read_xlsx('../../../../../HolgerImperial/18/Ribosomal proteins and translation factors.xlsx',
      sheet = 'Mit Ribosomal Prot'))[,1]
    plotheatmap('Mit Ribosomal Prot', width = 6.5, height = 12, genes, labelskip = 3)

    genes <- data.frame(readxl::read_xlsx('../../../../../HolgerImperial/18/Ribosomal proteins and translation factors.xlsx',
      sheet = 'Initiation Factors'))[,1]
    plotheatmap('Initiation Factors', width = 6.5, height = 12, genes, labelskip = 2)

    genes <- data.frame(readxl::read_xlsx('../../../../../HolgerImperial/18/Ribosomal proteins and translation factors.xlsx',
      sheet = 'Elongation factors'))[,1]
    plotheatmap('Elongation factors', width = 6.5, height = 8, genes, labelskip = 1)

    genes <- data.frame(readxl::read_xlsx('../../../../../HolgerImperial/18/Ribosomal proteins and translation factors.xlsx',
      sheet = 'Termination factors'))[,1]
    plotheatmap('Termination factors', width = 6.5, height = 6, genes, labelskip = 1)

    gsets_ref <- msigdbr(species = 'Homo sapiens') %>% split(x = .$gene_symbol, f = .$gs_name)
    gsets_ref_hallmark <- gsets_ref[grep('HALLMARK', names(gsets_ref))]
    genes <- gsets_ref[['HALLMARK_MYC_TARGETS_V1']]
    plotheatmap('HALLMARK_MYC_TARGETS_V1', width = 6.5, height = 12, genes, labelskip = 6)

    genes <- gsets_ref[['HALLMARK_MYC_TARGETS_V2']]
    plotheatmap('HALLMARK_MYC_TARGETS_V2', width = 6.5, height = 12, genes, labelskip = 2)



#####
#END#
#####
