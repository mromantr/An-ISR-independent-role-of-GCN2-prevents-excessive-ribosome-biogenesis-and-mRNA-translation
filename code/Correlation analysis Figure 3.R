
  dir.create('Exports', showWarnings = FALSE)
  dir.create('Plots', showWarnings = FALSE)

  rnaseq <- read.table('MetadataDoc/GiB.48h_vs_DMSO.48h_stats.tsv', header = TRUE, sep = '\t')
  protein <- read.table('MetadataDoc/proteome-GiB-vs-DMSO-ttest.csv', header = TRUE, sep = ',')

  # remove duplicates
    rnaseq <- rnaseq[!duplicated(rnaseq$symbol),]
    protein <- protein[!duplicated(protein$symbols),]

  require(ggplot2)
  require(ggrepel)
  require(grid)
  require(gridExtra)
  require(cowplot)
  require(dplyr)
  require(msigdbr)
  require(ComplexHeatmap)
  require(circlize)
  require(DESeq2)

  mytheme <- theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'plain', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'plain', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'plain', vjust = 1),
    axis.text.x = element_text(angle = 0, size = 12, face = 'plain', vjust = 0.5),
    axis.text.y = element_text(angle = 0, size = 12, face = 'plain', vjust = 0.5),
    axis.title = element_text(angle = 0, size = 14, face = 'plain'),
    legend.key = element_blank(),
    legend.key.size = unit(0.75, 'cm'),
    legend.text = element_text(angle = 0, size = 12, face = 'plain'),
    title = element_text(angle = 0, size = 12, face = 'plain'))

  metadata <- data.frame(readxl::read_xlsx('MetadataDoc/1260419__excel_tabless1-s18 1.xlsx',
    sheet = 'S1. All genes'))
    # has a 1 (or >1) in column SP and a 0 in TM --> it is SP
      sp_ens <- subset(metadata, SP >= 1 & TM == 0)$ensg_id
      sp <- subset(metadata, SP >= 1 & TM == 0)$display_name

    # has a 1 (or >1) in column TM and 0 in SP --> it is TM
      tm_ens <- subset(metadata, SP == 0 & TM >= 1)$ensg_id
      tm <- subset(metadata, SP == 0 & TM >= 1)$display_name

    # has a 1 (or >1) in columns SP and a 1 (or >1) in TM --> it is SP+TM
      sptm_ens <- subset(metadata, SP >= 1 & TM >= 1)$ensg_id
      sptm <- subset(metadata, SP >= 1 & TM >= 1)$display_name

    # has a 0 in column SP and a 0 in TM --> it is IC
      ic_ens <- subset(metadata, SP == 0 & TM == 0)$ensg_id
      ic <- subset(metadata, SP == 0 & TM == 0)$display_name

  # main results
      # focused analysis on SP, TM, SP+TM, IC
        # scatter plot of fold changes
          rnaseq.subset <- rnaseq
          protein.subset <- protein
          vec <- intersect(rnaseq.subset$symbol, protein.subset$symbols)
          vec <- vec[which(vec %in% c(sp, tm, sptm, ic))]
          rnaseq.subset <- subset(rnaseq.subset, symbol %in% vec)
          protein.subset <- subset(protein.subset, symbols %in% vec)
          protein.subset <- protein.subset[match(rnaseq.subset$symbol, protein.subset$symbols),]
          all(rnaseq.subset$symbol == protein.subset$symbols)

          ns <- rnaseq.subset[which(rnaseq.subset$padj >= 0.05 & protein.subset$padj >= 0.05),'symbol']

          ggdata <- data.frame(
            Symbol = rnaseq.subset$symbol,
            RNA = rnaseq.subset$log2FoldChange,
            Protein = protein.subset$log2FoldChange)
          ggdata$Category <- ggdata$Symbol
          ggdata$Category[ggdata$Category %in% sp] <- 'SP'
          ggdata$Category[ggdata$Category %in% tm] <- 'TM'
          ggdata$Category[ggdata$Category %in% sptm] <- 'SP+TM'
          ggdata$Category[ggdata$Category %in% ic] <- 'IC'
          ggdata$Category[ggdata$Symbol %in% ns] <- 'NS'
          # re-order to plot NS first
            ggdata <- rbind(
              ggdata[which(ggdata$Category %in% 'NS'),],
              ggdata[-which(ggdata$Category %in% 'NS'),])
          r <- round(cor(ggdata$RNA, ggdata$Protein), digits = 3)
          p <- cor.test(ggdata$RNA, ggdata$Protein)$p.value
          p <- ifelse(p<0.0001, 'p < 0.0001', paste0('p = ', round(p, digits = 4)))

          # export lists of proteins
            subset(ggdata, RNA == 0 & Protein == 0) # check nothing is zero

            # All
              write.table(
                data.frame(UpperRight = sort(subset(ggdata, RNA > 0 & Protein > 0)$Symbol)),
                paste0('Exports/ABCDE mRNA X protein XY Plot v1 All Upper-Right ', Sys.Date(), '.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)
              write.table(
                data.frame(UpperLeft = sort(subset(ggdata, RNA < 0 & Protein > 0)$Symbol)),
                paste0('Exports/ABCDE mRNA X protein XY Plot v1 All Upper-Left ', Sys.Date(), '.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)
              write.table(
                data.frame(BottomRight = sort(subset(ggdata, RNA > 0 & Protein < 0)$Symbol)),
                paste0('Exports/ABCDE mRNA X protein XY Plot v1 All Bottom-Right ', Sys.Date(), '.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)
              write.table(
                data.frame(BottomLeft = sort(subset(ggdata, RNA < 0 & Protein < 0)$Symbol)),
                paste0('Exports/ABCDE mRNA X protein XY Plot v1 All Bottom-Left ', Sys.Date(), '.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)

            # IC
              write.table(
                data.frame(UpperRight = sort(subset(ggdata, RNA > 0 & Protein > 0 & Category == 'IC')$Symbol)),
                paste0('Exports/ABCDE mRNA X protein XY Plot v1 IC Upper-Right ', Sys.Date(), '.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)
              write.table(
                data.frame(UpperLeft = sort(subset(ggdata, RNA < 0 & Protein > 0 & Category == 'IC')$Symbol)),
                paste0('Exports/ABCDE mRNA X protein XY Plot v1 IC Upper-Left ', Sys.Date(), '.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)
              write.table(
                data.frame(BottomRight = sort(subset(ggdata, RNA > 0 & Protein < 0 & Category == 'IC')$Symbol)),
                paste0('Exports/ABCDE mRNA X protein XY Plot v1 IC Bottom-Right ', Sys.Date(), '.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)
              write.table(
                data.frame(BottomLeft = sort(subset(ggdata, RNA < 0 & Protein < 0 & Category == 'IC')$Symbol)),
                paste0('Exports/ABCDE mRNA X protein XY Plot v1 IC Bottom-Left ', Sys.Date(), '.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)

            # SP
              write.table(
                data.frame(UpperRight = sort(subset(ggdata, RNA > 0 & Protein > 0 & Category == 'SP')$Symbol)),
                paste0('Exports/ABCDE mRNA X protein XY Plot v1 SP Upper-Right ', Sys.Date(), '.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)
              write.table(
                data.frame(UpperLeft = sort(subset(ggdata, RNA < 0 & Protein > 0 & Category == 'SP')$Symbol)),
                paste0('Exports/ABCDE mRNA X protein XY Plot v1 SP Upper-Left ', Sys.Date(), '.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)
              write.table(
                data.frame(BottomRight = sort(subset(ggdata, RNA > 0 & Protein < 0 & Category == 'SP')$Symbol)),
                paste0('Exports/ABCDE mRNA X protein XY Plot v1 SP Bottom-Right ', Sys.Date(), '.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)
              write.table(
                data.frame(BottomLeft = sort(subset(ggdata, RNA < 0 & Protein < 0 & Category == 'SP')$Symbol)),
                paste0('Exports/ABCDE mRNA X protein XY Plot v1 SP Bottom-Left ', Sys.Date(), '.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)

            # SP+TM
              write.table(
                data.frame(UpperRight = sort(subset(ggdata, RNA > 0 & Protein > 0 & Category == 'SP+TM')$Symbol)),
                paste0('Exports/ABCDE mRNA X protein XY Plot v1 SP+TM Upper-Right ', Sys.Date(), '.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)
              write.table(
                data.frame(UpperLeft = sort(subset(ggdata, RNA < 0 & Protein > 0 & Category == 'SP+TM')$Symbol)),
                paste0('Exports/ABCDE mRNA X protein XY Plot v1 SP+TM Upper-Left ', Sys.Date(), '.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)
              write.table(
                data.frame(BottomRight = sort(subset(ggdata, RNA > 0 & Protein < 0 & Category == 'SP+TM')$Symbol)),
                paste0('Exports/ABCDE mRNA X protein XY Plot v1 SP+TM Bottom-Right ', Sys.Date(), '.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)
              write.table(
                data.frame(BottomLeft = sort(subset(ggdata, RNA < 0 & Protein < 0 & Category == 'SP+TM')$Symbol)),
                paste0('Exports/ABCDE mRNA X protein XY Plot v1 SP+TM Bottom-Left ', Sys.Date(), '.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)

            # TM
              write.table(
                data.frame(UpperRight = sort(subset(ggdata, RNA > 0 & Protein > 0 & Category == 'TM')$Symbol)),
                paste0('Exports/ABCDE mRNA X protein XY Plot v1 TM Upper-Right ', Sys.Date(), '.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)
              write.table(
                data.frame(UpperLeft = sort(subset(ggdata, RNA < 0 & Protein > 0 & Category == 'TM')$Symbol)),
                paste0('Exports/ABCDE mRNA X protein XY Plot v1 TM Upper-Left ', Sys.Date(), '.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)
              write.table(
                data.frame(BottomRight = sort(subset(ggdata, RNA > 0 & Protein < 0 & Category == 'TM')$Symbol)),
                paste0('Exports/ABCDE mRNA X protein XY Plot v1 TM Bottom-Right ', Sys.Date(), '.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)
              write.table(
                data.frame(BottomLeft = sort(subset(ggdata, RNA < 0 & Protein < 0 & Category == 'TM')$Symbol)),
                paste0('Exports/ABCDE mRNA X protein XY Plot v1 TM Bottom-Left ', Sys.Date(), '.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)

          df <- data.frame(rbind(
            ggdata[-which(ggdata$Category == 'NS'),],
            ggdata[which(ggdata$Category == 'NS'),]))
          df$RNA <- as.numeric(df$RNA)
          df$Protein <- as.numeric(df$Protein)
          df$Category <- factor(df$Category)
          p <- ggplot(df, aes(x = RNA, y = Protein)) +
            geom_point(aes(colour = Category, fill = Category), size = 0.75) +
            geom_smooth(method = 'lm', formula = y~x, level = 0.95) +
            stat_smooth(method = 'lm', fullrange = TRUE, level = 0.95, colour = 'red2') +
            scale_color_manual(
              values = c('IC' = 'tomato1', 'SP' = 'forestgreen', 'TM' = 'purple', 'SP+TM' = 'turquoise2', 'NS' = 'grey50')) +
            scale_fill_manual(
              values = c('IC' = 'tomato1', 'SP' = 'forestgreen', 'TM' = 'purple', 'SP+TM' = 'turquoise2', 'NS' = 'grey50')) +
            theme_bw(base_size = 24) + mytheme + theme(
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank()) +
            guides(colour = guide_legend(override.aes = list(size = 2.5))) +
            ylab(bquote(Proteomics~log[2]~FC)) +
            xlab(bquote(RNA-seq~log[2]~FC)) +
            geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
            geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
            labs(title = NULL,
              subtitle = NULL,
              caption = paste('Pearson r = ', r, ' | ', p, sep = ''))
          pdf(paste0('Plots/ABCDE mRNA X protein XY Plot ', Sys.Date(), '.pdf'), width = 7, height = 5.5)
            cowplot::plot_grid(p)
          dev.off()
          png(paste0('Plots/ABCDE mRNA X protein XY Plot ', Sys.Date(), '.png'),
            units = 'in', res = 300, width = 7, height = 5.5)
            cowplot::plot_grid(p)
          dev.off()

          p1 <- ggplot(ggdata, aes(x = RNA, y = Protein)) +
            geom_point(aes(colour = Category, fill = Category), size = 0.75) +
            geom_smooth(method = 'lm', formula = y~x, level = 0.95) +
            stat_smooth(method = 'lm', fullrange = TRUE, level = 0.95, colour = 'red2') +
            scale_color_manual(
              values = c('IC' = 'grey60', 'SP' = 'blue2', 'TM' = 'red2', 'SP+TM' = 'forestgreen', 'NS' = 'grey80')) +
            scale_fill_manual(
              values = c('IC' = 'grey60', 'SP' = 'blue2', 'TM' = 'red2', 'SP+TM' = 'forestgreen', 'NS' = 'grey80')) +
            theme_bw(base_size = 24) + mytheme + theme(
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank()) +
            guides(colour = guide_legend(override.aes = list(size = 2.5))) +
            ylab(bquote(Proteomics~log[2]~FC)) +
            xlab(bquote(RNA-seq~log[2]~FC)) +
            geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
            geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
            labs(title = 'E all mRNA and proteins',
              subtitle = NULL,
              caption = paste('Pearson r = ', r, ' | ', p, sep = '')) +
            geom_text(data = data.frame(
              xpos = c(-Inf, -Inf, Inf, Inf),
              ypos =  c(-Inf, Inf, -Inf, Inf),
              annotateText = c('E', 'D', 'C', 'B'),
                hjustvar = c(-0.5, -0.5, 1.5, 1.5) ,
                vjustvar = c(-1, 2, -1, 2)),
              aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText),
              size = 9)
          pdf(paste0('Plots/ABCDE mRNA X protein XY Plot v1', Sys.Date(), '.pdf'), width = 7, height = 6)
            cowplot::plot_grid(p1)
          dev.off()
          pdf(paste0('Plots/ABCDE mRNA X protein XY Plot v1 Stratified ', Sys.Date(), '.pdf'), width = 16, height = 5)
            p1 + facet_grid(. ~ Category) + labs(caption = '')
          dev.off()

          p1 <- ggplot(ggdata, aes(x = RNA, y = Protein)) +
            geom_point(aes(colour = Category, fill = Category), size = 0.75) +
            scale_color_manual(
              values = c('IC' = 'grey60', 'SP' = 'blue2', 'TM' = 'red2', 'SP+TM' = 'forestgreen', 'NS' = 'grey80')) +
            scale_fill_manual(
              values = c('IC' = 'grey60', 'SP' = 'blue2', 'TM' = 'red2', 'SP+TM' = 'forestgreen', 'NS' = 'grey80')) +
            theme_bw(base_size = 24) + mytheme + theme(
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank()) +
            guides(colour = guide_legend(override.aes = list(size = 2.5))) +
            ylab(bquote(Proteomics~log[2]~FC)) +
            xlab(bquote(RNA-seq~log[2]~FC)) +
            geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
            geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
            labs(title = 'E all mRNA and proteins',
              subtitle = NULL,
              caption = paste('Pearson r = ', r, ' | ', p, sep = '')) +
            geom_text(data = data.frame(
              xpos = c(-Inf, -Inf, Inf, Inf),
              ypos =  c(-Inf, Inf, -Inf, Inf),
              annotateText = c('E', 'D', 'C', 'B'),
                hjustvar = c(-0.5, -0.5, 1.5, 1.5) ,
                vjustvar = c(-1, 2, -1, 2)),
              aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText),
              size = 9)
          pdf(paste0('Plots/ABCDE mRNA X protein XY Plot v2 ', Sys.Date(), '.pdf'), width = 7, height = 6)
            cowplot::plot_grid(p1)
          dev.off()
          pdf(paste0('Plots/ABCDE mRNA X protein XY Plot v2 Stratified ', Sys.Date(), '.pdf'), width = 16, height = 5)
            p1 + facet_grid(. ~ Category) + labs(caption = '')
          dev.off()

          p1 <- ggplot(ggdata, aes(x = RNA, y = Protein)) +
            geom_point(aes(colour = Category, fill = Category), size = 0.75) +
            scale_color_manual(
              values = c('IC' = 'grey60', 'SP' = 'blue2', 'TM' = 'red2', 'SP+TM' = 'forestgreen', 'NS' = 'grey80')) +
            scale_fill_manual(
              values = c('IC' = 'grey60', 'SP' = 'blue2', 'TM' = 'red2', 'SP+TM' = 'forestgreen', 'NS' = 'grey80')) +
            # top 5 from each category (B-E)
              geom_label_repel(data = subset(ggdata, Symbol %in% c(top5B, top5C, top5D, top5E)),
                aes(x = RNA, y = Protein, label = Symbol),
                size = 3,  max.overlaps = Inf) +
            theme_bw(base_size = 24) + mytheme + theme(
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank()) +
            guides(colour = guide_legend(override.aes = list(size = 2.5))) +
            ylab(bquote(Proteomics~log[2]~FC)) +
            xlab(bquote(RNA-seq~log[2]~FC)) +
            geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
            geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
            labs(title = 'E all mRNA and proteins',
              subtitle = NULL,
              caption = paste('Pearson r = ', r, ' | ', p, sep = '')) +
            geom_text(data = data.frame(
              xpos = c(-Inf, -Inf, Inf, Inf),
              ypos =  c(-Inf, Inf, -Inf, Inf),
              annotateText = c('E', 'D', 'C', 'B'),
                hjustvar = c(-0.5, -0.5, 1.5, 1.5) ,
                vjustvar = c(-1, 2, -1, 2)),
              aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText),
              size = 9)
          pdf(paste0('Plots/ABCDE mRNA X protein XY Plot v3 ', Sys.Date(), '.pdf'), width = 7, height = 6)
            cowplot::plot_grid(p1)
          dev.off()
          pdf(paste0('Plots/ABCDE mRNA X protein XY Plot v3 Stratified ', Sys.Date(), '.pdf'), width = 16, height = 5)
            p1 + facet_grid(. ~ Category) + labs(caption = '')
          dev.off()



  # a set of correlation plots, no grid, no correlation line, all genes in medium-light grey. Per plot, highlight one set of genes of the following sets (each in medium purple):  
    # GO:0042254 Ribosome biogenesis
      entrez <- gsets_ref[grep(toupper('GO_Ribosome_biogenesis$'), names(gsets_ref))]
      symbol <- annotMaster[match(entrez[[1]], annotMaster$entrezgene_id),'hgnc_symbol']

      rnaseq.subset <- rnaseq
      protein.subset <- protein
      vec <- intersect(rnaseq.subset$symbol, protein.subset$symbols)
      vec <- vec[which(vec %in% c(sp, tm, sptm, ic))]
      rnaseq.subset <- subset(rnaseq.subset, symbol %in% vec)
      protein.subset <- subset(protein.subset, symbols %in% vec)
      protein.subset <- protein.subset[match(rnaseq.subset$symbol, protein.subset$symbols),]
      all(rnaseq.subset$symbol == protein.subset$symbols)

      ggdata <- data.frame(
        Symbol = rnaseq.subset$symbol,
        RNA = rnaseq.subset$log2FoldChange,
        Protein = protein.subset$log2FoldChange)
      ggdata$Category <- ggdata$Symbol
      ggdata$Category[ggdata$Category %in% symbol] <- 'Ribosome biogenesis'
      ggdata$Category[ggdata$Category != 'Ribosome biogenesis'] <- 'non-Ribosome\nbiogenesis'
      # re-order to plot nonGO first
        ggdata <- rbind(
          ggdata[which(ggdata$Category %in% 'non-Ribosome\nbiogenesis'),],
          ggdata[-which(ggdata$Category %in% 'non-Ribosome\nbiogenesis'),])
      p1 <- ggplot(ggdata, aes(x = RNA, y = Protein)) +
        geom_point(aes(colour = Category, fill = Category)) +
        scale_color_manual(
          values = c('non-Ribosome\nbiogenesis' = 'grey50', 'Ribosome biogenesis' = 'purple')) +
        scale_fill_manual(
          values = c('non-Ribosome\nbiogenesis' = 'grey50', 'Ribosome biogenesis' = 'purple')) +
        theme_bw(base_size = 24) + mytheme + theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 2.5))) +
        ylab(bquote(Proteomics~log[2]~FC)) +
        xlab(bquote(RNA-seq~log[2]~FC)) +
        geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
        geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
        labs(title = 'E all mRNA and proteins',
          subtitle = 'GO:0042254 Ribosome biogenesis',
          caption = NULL)
      pdf(paste0('Plots/ABCDE GO0042254 Ribosome biogenesis ', Sys.Date(), '.pdf'), width = 7, height = 5)
        cowplot::plot_grid(p1)
      dev.off()

    # GO:0022613 Ribonucleoprotein complex biogenesis
      entrez <- gsets_ref[grep(toupper('GO_Ribonucleoprotein_complex_biogenesis$'), names(gsets_ref))]
      symbol <- annotMaster[match(entrez[[1]], annotMaster$entrezgene_id),'hgnc_symbol']

      rnaseq.subset <- rnaseq
      protein.subset <- protein
      vec <- intersect(rnaseq.subset$symbol, protein.subset$symbols)
      vec <- vec[which(vec %in% c(sp, tm, sptm, ic))]
      rnaseq.subset <- subset(rnaseq.subset, symbol %in% vec)
      protein.subset <- subset(protein.subset, symbols %in% vec)
      protein.subset <- protein.subset[match(rnaseq.subset$symbol, protein.subset$symbols),]
      all(rnaseq.subset$symbol == protein.subset$symbols)

      ggdata <- data.frame(
        Symbol = rnaseq.subset$symbol,
        RNA = rnaseq.subset$log2FoldChange,
        Protein = protein.subset$log2FoldChange)
      ggdata$Category <- ggdata$Symbol
      ggdata$Category[ggdata$Category %in% symbol] <- 'Ribonucleoprotein\ncomplex biogenesis'
      ggdata$Category[ggdata$Category != 'Ribonucleoprotein\ncomplex biogenesis'] <- 'non-Ribonucleoprotein\ncomplex biogenesis'
      # re-order to plot nonGO first
        ggdata <- rbind(
          ggdata[which(ggdata$Category %in% 'non-Ribonucleoprotein\ncomplex biogenesis'),],
          ggdata[-which(ggdata$Category %in% 'non-Ribonucleoprotein\ncomplex biogenesis'),])
      p1 <- ggplot(ggdata, aes(x = RNA, y = Protein)) +
        geom_point(aes(colour = Category, fill = Category)) +
        scale_color_manual(
          values = c('non-Ribonucleoprotein\ncomplex biogenesis' = 'grey50', 'Ribonucleoprotein\ncomplex biogenesis' = 'purple')) +
        scale_fill_manual(
          values = c('non-Ribonucleoprotein\ncomplex biogenesis' = 'grey50', 'Ribonucleoprotein\ncomplex biogenesis' = 'purple')) +
        theme_bw(base_size = 24) + mytheme + theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 2.5))) +
        ylab(bquote(Proteomics~log[2]~FC)) +
        xlab(bquote(RNA-seq~log[2]~FC)) +
        geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
        geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
        labs(title = 'E all mRNA and proteins',
          subtitle = 'GO:0022613 Ribonucleoprotein complex biogenesis',
          caption = NULL)
      pdf(paste0('Plots/ABCDE GO0022613 Ribonucleoprotein complex biogenesis ', Sys.Date(), '.pdf'), width = 7.5, height = 5.5)
        cowplot::plot_grid(p1)
      dev.off()

    # GO:0005730 Nucleolus
      entrez <- gsets_ref[grep(toupper('GO_Nucleolus$'), names(gsets_ref))]
      symbol <- annotMaster[match(entrez[[1]], annotMaster$entrezgene_id),'hgnc_symbol']

      rnaseq.subset <- rnaseq
      protein.subset <- protein
      vec <- intersect(rnaseq.subset$symbol, protein.subset$symbols)
      vec <- vec[which(vec %in% c(sp, tm, sptm, ic))]
      rnaseq.subset <- subset(rnaseq.subset, symbol %in% vec)
      protein.subset <- subset(protein.subset, symbols %in% vec)
      protein.subset <- protein.subset[match(rnaseq.subset$symbol, protein.subset$symbols),]
      all(rnaseq.subset$symbol == protein.subset$symbols)

      ggdata <- data.frame(
        Symbol = rnaseq.subset$symbol,
        RNA = rnaseq.subset$log2FoldChange,
        Protein = protein.subset$log2FoldChange)
      ggdata$Category <- ggdata$Symbol
      ggdata$Category[ggdata$Category %in% symbol] <- 'Nucleolus'
      ggdata$Category[ggdata$Category != 'Nucleolus'] <- 'non-Nucleolus'
      # re-order to plot nonGO first
        ggdata <- rbind(
          ggdata[which(ggdata$Category %in% 'non-Nucleolus'),],
          ggdata[-which(ggdata$Category %in% 'non-Nucleolus'),])
      p1 <- ggplot(ggdata, aes(x = RNA, y = Protein)) +
        geom_point(aes(colour = Category, fill = Category)) +
        scale_color_manual(
          values = c('non-Nucleolus' = 'grey50', 'Nucleolus' = 'purple')) +
        scale_fill_manual(
          values = c('non-Nucleolus' = 'grey50', 'Nucleolus' = 'purple')) +
        theme_bw(base_size = 24) + mytheme + theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 2.5))) +
        ylab(bquote(Proteomics~log[2]~FC)) +
        xlab(bquote(RNA-seq~log[2]~FC)) +
        geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
        geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
        labs(title = 'E all mRNA and proteins',
          subtitle = 'GO:0005730 Nucleolus',
          caption = NULL)
      pdf(paste0('Plots/ABCDE GO0005730 Nucleolus ', Sys.Date(), '.pdf'), width = 7, height = 5.5)
        cowplot::plot_grid(p1)
      dev.off()

    # GO:0030684 Preribosome
      entrez <- gsets_ref[grep(toupper('GO_Preribosome$'), names(gsets_ref))]
      symbol <- annotMaster[match(entrez[[1]], annotMaster$entrezgene_id),'hgnc_symbol']

      rnaseq.subset <- rnaseq
      protein.subset <- protein
      vec <- intersect(rnaseq.subset$symbol, protein.subset$symbols)
      vec <- vec[which(vec %in% c(sp, tm, sptm, ic))]
      rnaseq.subset <- subset(rnaseq.subset, symbol %in% vec)
      protein.subset <- subset(protein.subset, symbols %in% vec)
      protein.subset <- protein.subset[match(rnaseq.subset$symbol, protein.subset$symbols),]
      all(rnaseq.subset$symbol == protein.subset$symbols)

      ggdata <- data.frame(
        Symbol = rnaseq.subset$symbol,
        RNA = rnaseq.subset$log2FoldChange,
        Protein = protein.subset$log2FoldChange)
      ggdata$Category <- ggdata$Symbol
      ggdata$Category[ggdata$Category %in% symbol] <- 'Preribosome'
      ggdata$Category[ggdata$Category != 'Preribosome'] <- 'non-Preribosome'
      # re-order to plot nonGO first
        ggdata <- rbind(
          ggdata[which(ggdata$Category %in% 'non-Preribosome'),],
          ggdata[-which(ggdata$Category %in% 'non-Preribosome'),])
      p1 <- ggplot(ggdata, aes(x = RNA, y = Protein)) +
        geom_point(aes(colour = Category, fill = Category)) +
        scale_color_manual(
          values = c('non-Preribosome' = 'grey50', 'Preribosome' = 'purple')) +
        scale_fill_manual(
          values = c('non-Preribosome' = 'grey50', 'Preribosome' = 'purple')) +
        theme_bw(base_size = 24) + mytheme + theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 2.5))) +
        ylab(bquote(Proteomics~log[2]~FC)) +
        xlab(bquote(RNA-seq~log[2]~FC)) +
        geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
        geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
        labs(title = 'E all mRNA and proteins',
          subtitle = 'GO:0030684 Preribosome',
          caption = NULL)
      pdf(paste0('Plots/ABCDE GO0030684 Preribosome ', Sys.Date(), '.pdf'), width = 7, height = 5.5)
        cowplot::plot_grid(p1)
      dev.off()

    # GO:0000502 Proteasome complex
      # gsets_ref[grep(toupper('GO_Proteasome_complex$'), names(gsets_ref))]
      symbol <- read.table('MetadataDoc/GO0000502_proteasome_complex.list', header = FALSE)[,1]

      rnaseq.subset <- rnaseq
      protein.subset <- protein
      vec <- intersect(rnaseq.subset$symbol, protein.subset$symbols)
      vec <- vec[which(vec %in% c(sp, tm, sptm, ic))]
      rnaseq.subset <- subset(rnaseq.subset, symbol %in% vec)
      protein.subset <- subset(protein.subset, symbols %in% vec)
      protein.subset <- protein.subset[match(rnaseq.subset$symbol, protein.subset$symbols),]
      all(rnaseq.subset$symbol == protein.subset$symbols)

      ggdata <- data.frame(
        Symbol = rnaseq.subset$symbol,
        RNA = rnaseq.subset$log2FoldChange,
        Protein = protein.subset$log2FoldChange)
      ggdata$Category <- ggdata$Symbol
      ggdata$Category[ggdata$Category %in% symbol] <- 'Proteasome complex'
      ggdata$Category[ggdata$Category != 'Proteasome complex'] <- 'non-Proteasome\ncomplex'
      # re-order to plot nonGO first
        ggdata <- rbind(
          ggdata[which(ggdata$Category %in% 'non-Proteasome\ncomplex'),],
          ggdata[-which(ggdata$Category %in% 'non-Proteasome\ncomplex'),])
      p1 <- ggplot(ggdata, aes(x = RNA, y = Protein)) +
        geom_point(aes(colour = Category, fill = Category)) +
        scale_color_manual(
          values = c('non-Proteasome\ncomplex' = 'grey50', 'Proteasome complex' = 'purple')) +
        scale_fill_manual(
          values = c('non-Proteasome\ncomplex' = 'grey50', 'Proteasome complex' = 'purple')) +
        theme_bw(base_size = 24) + mytheme + theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 2.5))) +
        ylab(bquote(Proteomics~log[2]~FC)) +
        xlab(bquote(RNA-seq~log[2]~FC)) +
        geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
        geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
        labs(title = 'E all mRNA and proteins',
          subtitle = 'GO:0000502 Proteasome complex',
          caption = NULL)
      pdf(paste0('Plots/ABCDE GO0000502 Proteasome complex ', Sys.Date(), '.pdf'), width = 7, height = 5)
        cowplot::plot_grid(p1)
      dev.off()

    # GO:0005839 Proteasome core complex
      entrez <- gsets_ref[grep(toupper('GO_Proteasome_core_complex$'), names(gsets_ref))]
      symbol <- annotMaster[match(entrez[[1]], annotMaster$entrezgene_id),'hgnc_symbol']

      rnaseq.subset <- rnaseq
      protein.subset <- protein
      vec <- intersect(rnaseq.subset$symbol, protein.subset$symbols)
      vec <- vec[which(vec %in% c(sp, tm, sptm, ic))]
      rnaseq.subset <- subset(rnaseq.subset, symbol %in% vec)
      protein.subset <- subset(protein.subset, symbols %in% vec)
      protein.subset <- protein.subset[match(rnaseq.subset$symbol, protein.subset$symbols),]
      all(rnaseq.subset$symbol == protein.subset$symbols)

      ggdata <- data.frame(
        Symbol = rnaseq.subset$symbol,
        RNA = rnaseq.subset$log2FoldChange,
        Protein = protein.subset$log2FoldChange)
      ggdata$Category <- ggdata$Symbol
      ggdata$Category[ggdata$Category %in% symbol] <- 'Proteasome\ncore complex'
      ggdata$Category[ggdata$Category != 'Proteasome\ncore complex'] <- 'non-Proteasome\ncore complex'
      # re-order to plot nonGO first
        ggdata <- rbind(
          ggdata[which(ggdata$Category %in% 'non-Proteasome\ncore complex'),],
          ggdata[-which(ggdata$Category %in% 'non-Proteasome\ncore complex'),])
      p1 <- ggplot(ggdata, aes(x = RNA, y = Protein)) +
        geom_point(aes(colour = Category, fill = Category)) +
        scale_color_manual(
          values = c('non-Proteasome\ncore complex' = 'grey50', 'Proteasome\ncore complex' = 'purple')) +
        scale_fill_manual(
          values = c('non-Proteasome\ncore complex' = 'grey50', 'Proteasome\ncore complex' = 'purple')) +
        theme_bw(base_size = 24) + mytheme + theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 2.5))) +
        ylab(bquote(Proteomics~log[2]~FC)) +
        xlab(bquote(RNA-seq~log[2]~FC)) +
        geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
        geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
        labs(title = 'E all mRNA and proteins',
          subtitle = 'GO:0005839 Proteasome core complex',
          caption = NULL)
      pdf(paste0('Plots/ABCDE GO0005839 Proteasome core complex ', Sys.Date(), '.pdf'), width = 7.5, height = 5.5)
        cowplot::plot_grid(p1)
      dev.off()

    # GO:0061630 Ubiquitin protein ligase activity
      # gsets_ref[grep(toupper('GO_Ubiquitin_protein_ligase_activity$'), names(gsets_ref))]
      symbol <- read.table('MetadataDoc/GO0061630_ubiquitin_protein_ligase_activity.list', header = FALSE)[,1]

      rnaseq.subset <- rnaseq
      protein.subset <- protein
      vec <- intersect(rnaseq.subset$symbol, protein.subset$symbols)
      vec <- vec[which(vec %in% c(sp, tm, sptm, ic))]
      rnaseq.subset <- subset(rnaseq.subset, symbol %in% vec)
      protein.subset <- subset(protein.subset, symbols %in% vec)
      protein.subset <- protein.subset[match(rnaseq.subset$symbol, protein.subset$symbols),]
      all(rnaseq.subset$symbol == protein.subset$symbols)

      ggdata <- data.frame(
        Symbol = rnaseq.subset$symbol,
        RNA = rnaseq.subset$log2FoldChange,
        Protein = protein.subset$log2FoldChange)
      ggdata$Category <- ggdata$Symbol
      ggdata$Category[ggdata$Category %in% symbol] <- 'Ubiquitin protein\nligase activity'
      ggdata$Category[ggdata$Category != 'Ubiquitin protein\nligase activity'] <- 'non-Ubiquitin protein\nligase activity'
      # re-order to plot nonGO first
        ggdata <- rbind(
          ggdata[which(ggdata$Category %in% 'non-Ubiquitin protein\nligase activity'),],
          ggdata[-which(ggdata$Category %in% 'non-Ubiquitin protein\nligase activity'),])
      p1 <- ggplot(ggdata, aes(x = RNA, y = Protein)) +
        geom_point(aes(colour = Category, fill = Category)) +
        scale_color_manual(
          values = c('non-Ubiquitin protein\nligase activity' = 'grey50', 'Ubiquitin protein\nligase activity' = 'purple')) +
        scale_fill_manual(
          values = c('non-Ubiquitin protein\nligase activity' = 'grey50', 'Ubiquitin protein\nligase activity' = 'purple')) +
        theme_bw(base_size = 24) + mytheme + theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 2.5))) +
        ylab(bquote(Proteomics~log[2]~FC)) +
        xlab(bquote(RNA-seq~log[2]~FC)) +
        geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
        geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
        labs(title = 'E all mRNA and proteins',
          subtitle = 'GO:0061630 Ubiquitin protein ligase activity',
          caption = NULL)
      pdf(paste0('Plots/ABCDE GO0061630 Ubiquitin protein ligase activity ', Sys.Date(), '.pdf'), width = 7.5, height = 5.5)
        cowplot::plot_grid(p1)
      dev.off()

    # GO:0030198 Extracellular matrix organization
      # gsets_ref[grep(toupper('GO_Extracellular_matrix_organization$'), names(gsets_ref))]
      symbol <- read.table('MetadataDoc/GO0030198_extracellular_matrix_organization.list', header = FALSE)[,1]

      rnaseq.subset <- rnaseq
      protein.subset <- protein
      vec <- intersect(rnaseq.subset$symbol, protein.subset$symbols)
      vec <- vec[which(vec %in% c(sp, tm, sptm, ic))]
      rnaseq.subset <- subset(rnaseq.subset, symbol %in% vec)
      protein.subset <- subset(protein.subset, symbols %in% vec)
      protein.subset <- protein.subset[match(rnaseq.subset$symbol, protein.subset$symbols),]
      all(rnaseq.subset$symbol == protein.subset$symbols)

      ggdata <- data.frame(
        Symbol = rnaseq.subset$symbol,
        RNA = rnaseq.subset$log2FoldChange,
        Protein = protein.subset$log2FoldChange)
      ggdata$Category <- ggdata$Symbol
      ggdata$Category[ggdata$Category %in% symbol] <- 'Extracellular\nmatrix organization'
      ggdata$Category[ggdata$Category != 'Extracellular\nmatrix organization'] <- 'non-Extracellular\nmatrix organization'
      # re-order to plot nonGO first
        ggdata <- rbind(
          ggdata[which(ggdata$Category %in% 'non-Extracellular\nmatrix organization'),],
          ggdata[-which(ggdata$Category %in% 'non-Extracellular\nmatrix organization'),])
      p1 <- ggplot(ggdata, aes(x = RNA, y = Protein)) +
        geom_point(aes(colour = Category, fill = Category)) +
        scale_color_manual(
          values = c('non-Extracellular\nmatrix organization' = 'grey50', 'Extracellular\nmatrix organization' = 'purple')) +
        scale_fill_manual(
          values = c('non-Extracellular\nmatrix organization' = 'grey50', 'Extracellular\nmatrix organization' = 'purple')) +
        theme_bw(base_size = 24) + mytheme + theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 2.5))) +
        ylab(bquote(Proteomics~log[2]~FC)) +
        xlab(bquote(RNA-seq~log[2]~FC)) +
        geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
        geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
        labs(title = 'E all mRNA and proteins',
          subtitle = 'GO:0030198 Extracellular matrix organization',
          caption = NULL)
      pdf(paste0('Plots/ABCDE GO0030198 Extracellular matrix organization ', Sys.Date(), '.pdf'), width = 7.5, height = 5.5)
        cowplot::plot_grid(p1)
      dev.off()


    # GO:0043062 Extracellular structure organization
      entrez <- gsets_ref[grep(toupper('GO_Extracellular_structure_organization$'), names(gsets_ref))]
      symbol <- annotMaster[match(entrez[[1]], annotMaster$entrezgene_id),'hgnc_symbol']

      rnaseq.subset <- rnaseq
      protein.subset <- protein
      vec <- intersect(rnaseq.subset$symbol, protein.subset$symbols)
      vec <- vec[which(vec %in% c(sp, tm, sptm, ic))]
      rnaseq.subset <- subset(rnaseq.subset, symbol %in% vec)
      protein.subset <- subset(protein.subset, symbols %in% vec)
      protein.subset <- protein.subset[match(rnaseq.subset$symbol, protein.subset$symbols),]
      all(rnaseq.subset$symbol == protein.subset$symbols)

      ggdata <- data.frame(
        Symbol = rnaseq.subset$symbol,
        RNA = rnaseq.subset$log2FoldChange,
        Protein = protein.subset$log2FoldChange)
      ggdata$Category <- ggdata$Symbol
      ggdata$Category[ggdata$Category %in% symbol] <- 'Extracellular\nstructure organization'
      ggdata$Category[ggdata$Category != 'Extracellular\nstructure organization'] <- 'non-Extracellular\nstructure organization'
      # re-order to plot nonGO first
        ggdata <- rbind(
          ggdata[which(ggdata$Category %in% 'non-Extracellular\nstructure organization'),],
          ggdata[-which(ggdata$Category %in% 'non-Extracellular\nstructure organization'),])
      p1 <- ggplot(ggdata, aes(x = RNA, y = Protein)) +
        geom_point(aes(colour = Category, fill = Category)) +
        scale_color_manual(
          values = c('non-Extracellular\nstructure organization' = 'grey50', 'Extracellular\nstructure organization' = 'purple')) +
        scale_fill_manual(
          values = c('non-Extracellular\nstructure organization' = 'grey50', 'Extracellular\nstructure organization' = 'purple')) +
        theme_bw(base_size = 24) + mytheme + theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 2.5))) +
        ylab(bquote(Proteomics~log[2]~FC)) +
        xlab(bquote(RNA-seq~log[2]~FC)) +
        geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
        geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
        labs(title = 'E all mRNA and proteins',
          subtitle = 'GO:0043062 Extracellular structure organization',
          caption = NULL)
      pdf(paste0('Plots/ABCDE GO0043062 Extracellular structure organization ', Sys.Date(), '.pdf'), width = 7.5, height = 5.5)
        cowplot::plot_grid(p1)
      dev.off()

    # GO:0005615 Extracellular space
      # gsets_ref[grep(toupper('GO_Extracellular_space$'), names(gsets_ref))]
      symbol <- read.table('MetadataDoc/GO0005615_extracellular_space.list', header = FALSE)[,1]

      rnaseq.subset <- rnaseq
      protein.subset <- protein
      vec <- intersect(rnaseq.subset$symbol, protein.subset$symbols)
      vec <- vec[which(vec %in% c(sp, tm, sptm, ic))]
      rnaseq.subset <- subset(rnaseq.subset, symbol %in% vec)
      protein.subset <- subset(protein.subset, symbols %in% vec)
      protein.subset <- protein.subset[match(rnaseq.subset$symbol, protein.subset$symbols),]
      all(rnaseq.subset$symbol == protein.subset$symbols)

      ggdata <- data.frame(
        Symbol = rnaseq.subset$symbol,
        RNA = rnaseq.subset$log2FoldChange,
        Protein = protein.subset$log2FoldChange)
      ggdata$Category <- ggdata$Symbol
      ggdata$Category[ggdata$Category %in% symbol] <- 'Extracellular\nspace'
      ggdata$Category[ggdata$Category != 'Extracellular\nspace'] <- 'non-Extracellular\nspace'
      # re-order to plot nonGO first
        ggdata <- rbind(
          ggdata[which(ggdata$Category %in% 'non-Extracellular\nspace'),],
          ggdata[-which(ggdata$Category %in% 'non-Extracellular\nspace'),])
      p1 <- ggplot(ggdata, aes(x = RNA, y = Protein)) +
        geom_point(aes(colour = Category, fill = Category)) +
        scale_color_manual(
          values = c('non-Extracellular\nspace' = 'grey50', 'Extracellular\nspace' = 'purple')) +
        scale_fill_manual(
          values = c('non-Extracellular\nspace' = 'grey50', 'Extracellular\nspace' = 'purple')) +
        theme_bw(base_size = 24) + mytheme + theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 2.5))) +
        ylab(bquote(Proteomics~log[2]~FC)) +
        xlab(bquote(RNA-seq~log[2]~FC)) +
        geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
        geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
        labs(title = 'E all mRNA and proteins',
          subtitle = 'GO:0005615 Extracellular space',
          caption = NULL)
      pdf(paste0('Plots/ABCDE GO0005615 Extracellular space ', Sys.Date(), '.pdf'), width = 7.5, height = 5.5)
        cowplot::plot_grid(p1)
      dev.off()
