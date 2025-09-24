# DESCRIPTION ----------------------------------------------------------------------------------
# Author: Elizabeth Hutchins
# Date: January 26th, 2023

#Purpose: gtex summary - genes detected

#Imports:
# aggregated count table (GTeX)
# biotype files

#Exports: 
#

# import packages -------------------------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(edgeR)
#toolkitPath <- "~/Dropbox/repos/toolkit_ehutchins/"
#source(paste(toolkitPath, "R-functions/createColorPalettes.R", sep = ""))
#source(paste(toolkitPath, "R-functions/biotypePlotter.R", sep = ""))
library(rtracklayer)
library(devtools)

#R package to calculate tau
#https://rdrr.io/github/roonysgalbi/tispec/f/vignettes/UserGuide.Rmd
#install.packages("remotes")
#remotes::install_github("roonysgalbi/tispec")
library(tispec)

source_gist("c579c4ddc06fd2ceb121e690dd9df186") #color palettes
myPalette <- getPalette(4, "sequential")
dev.off()

#function updates -----
getDist <- function (x, y) 
{
  tissues <- names(x[, c(-1, -2, -3)])
  dist <- data.frame()
  for (i in tissues) {
    asg <- subset(x, x[, i] >= y)
    t <- length(row.names(asg))
    p <- length(row.names(subset(asg, asg$gene_biotype == 
                                   "protein_coding")))
    m <- length(row.names(subset(asg, asg$gene_biotype == 
                                   "long_noncoding")))
    l <- length(row.names(subset(asg, asg$gene_biotype == 
                                   "short_noncoding")))
    n <- length(row.names(subset(asg, asg$gene_biotype == 
                                   "pseudogene")))
    q <- length(row.names(subset(asg, asg$gene_biotype == 
                                   "TEC")))
    newrow = c(p, m, l, n, q, t)
    dist <- rbind(dist, newrow)
  }
  colnames(dist) <- c("protein_coding", "long_noncoding", "short_noncoding", "pseudogene", "TEC",
                      "total")
  row.names(dist) <- tissues
  return(dist)
}

# load gtex count table with metadata ----------------------------------------------------------------------------------
gtex.cnts <- read_tsv("data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", skip = 2)
attributes <- read_tsv("data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", col_names = TRUE)


# incorporate biotypes and gene names ----------------------------------------------------------
#load annotation info for GENCODE v19 used for gtex
#generated like so:
gencode.v26.gtf <- "data/gencode.v26.GRCh38.genes.gtf"
gtf <- rtracklayer::import(gencode.v26.gtf, format = "gtf")
tmp <- as.data.frame(gtf)

tmp <- subset(tmp, type == "gene")

gencode.anno <- tmp[,c("gene_id", "gene_name", "gene_type")]
gencode.anno <- unique(gencode.anno)

#consolidate biotypes
pcBiotypes <- c("protein_coding", "IG_V_gene", "IG_C_gene", "IG_J_gene",
                "TR_C_gene", "TR_J_gene", "TR_V_gene", "TR_D_gene",
                "IG_D_gene", "nonsense_mediated_decay", "non_stop_decay")
lncBiotypes <- c("lincRNA", "antisense",
                 "non_coding",
                 "processed_transcript",
                 "3prime_overlapping_ncRNA",
                 "bidirectional_promoter_lncRNA", "macro_lncRNA",
                 "sense_intronic", "sense_overlapping",
                 "retained_intron", "lncipedia_lncRNA"
)
short.ncBiotypes <- c("snRNA", "snoRNA", "misc_RNA", "miRNA", "rRNA",
                      "Mt_tRNA", "Mt_rRNA", "ribozyme", "scaRNA",
                      "vaultRNA", "sRNA", "scRNA")
pseudogeneBiotypes <- c("processed_pseudogene", "transcribed_unprocessed_pseudogene",
                        "unprocessed_pseudogene", "transcribed_processed_pseudogene",
                        "transcribed_unitary_pseudogene", "rRNA_pseudogene",
                        "unitary_pseudogene", "polymorphic_pseudogene",
                        "pseudogene", "IG_V_pseudogene",
                        "translated_processed_pseudogene", "TR_V_pseudogene",
                        "IG_C_pseudogene", "TR_J_pseudogene",
                        "IG_J_pseudogene", "IG_pseudogene")
TECbiotypes <- c("TEC")

gencode.anno$gene_biotype == NULL
gencode.anno$gene_biotype[gencode.anno$gene_type %in% pcBiotypes] <- "protein_coding"
gencode.anno$gene_biotype[gencode.anno$gene_type %in% lncBiotypes] <- "long_noncoding"
gencode.anno$gene_biotype[gencode.anno$gene_type %in% short.ncBiotypes] <- "short_noncoding"
gencode.anno$gene_biotype[gencode.anno$gene_type %in% pseudogeneBiotypes] <- "pseudogene"
gencode.anno$gene_biotype[gencode.anno$gene_type %in% TECbiotypes] <- "TEC"
unique(gencode.anno$gene_biotype)

#countMatGtex.wb.anno <- right_join(gencode.anno, countMatGtex.wb, by = c("gene_id" = "Name"))
#head(countMatGtex.wb.anno[,c(1:10)])


#CPM normalization
countTable <- gtex.cnts %>% select(-Description) %>%
  column_to_rownames(var = "Name")
cpmTable <- cpm(countTable)
cpmTable.longer <- cpmTable %>% as.data.frame() %>% rownames_to_column(var = "gene_id") %>%
  pivot_longer(-gene_id, values_to = "CPM", names_to = "SAMPID") %>%
  left_join(attributes %>% select(SAMPID, SMTS, SMTSD))

cpmTable.longer.filt <- cpmTable.longer %>% filter(!SMTSD %in% c("Cells - EBV-transformed lymphocytes",
                                                                 "Cells - Leukemia cell line (CML)",
                                                                 "Cells - Cultured fibroblasts"))
#dim(cpmTable.longer)
#dim(cpmTable.longer.filt)
#unique(cpmTable.longer.filt$SMTSD)

cpmTable.meanTissue <- cpmTable.longer.filt %>%
  group_by(gene_id, SMTS) %>%
  summarise(meanCPM = mean(CPM)) %>%
  pivot_wider(names_from = "SMTS", values_from = "meanCPM")

#normalization
#meanExp <- cpmTable.meanTissue %>% column_to_rownames(var = "gene_id")
meanExp <- cpmTable.meanTissue %>%
  select(-Breast, -`Cervix Uteri`, -`Fallopian Tube`, -Ovary, -Prostate, -Testis, -Uterus, -Vagina) %>%
  column_to_rownames(var = "gene_id")
log2Exp <- log2Tran(meanExp) 
head(log2Exp[,1:5], n=5)
qnExp <- quantNorm(log2Exp)
head(qnExp[,1:5], n=5)

#calculate tau
tauExp <- calcTau(qnExp) 
head(tauExp[,1:5], n=5)

tauAnno <- right_join(gencode.anno, rownames_to_column(tauExp, var = "gene_id")) %>%
  column_to_rownames(var = "gene_id") %>%
  rename("gene_name" = "external_gene_name") %>%
  select(-gene_type)

#tau distribution
plotDensity(tauAnno) 
ggsave("tau_density_noReproductive.png", width = 8, height = 6)

asg <- getDist(tauAnno, 1)
head(asg, n = 5)

hsg <- getDist(tauAnno, 0.95)
head(hsg, n = 5)
plotDist(tauAnno)

#plot distribution - number of genes at different cutoffs ---------
plotDistBiotype <- function(tauAnno, biotype){
  tau1.00 <- getDist(tauAnno, 1)[, biotype]
  tau0.99 <- getDist(tauAnno, 0.99)[, biotype]
  tau0.95 <- getDist(tauAnno, 0.95)[, biotype]
  tau0.90 <- getDist(tauAnno, 0.90)[, biotype]
  dist <- data.frame(cbind(tau1.00, tau0.99, tau0.95, tau0.90))
  tissues <- names(tauAnno[, c(-1, -2, -3)])
  dist$tissues <- tissues
  dist <- tidyr::gather(dist, "spec", "count", -tissues)
  ggplot2::ggplot(dist, ggplot2::aes_string(x = "tissues", 
                                            y = "count", fill = "spec")) +
    ggplot2::geom_bar(position = "dodge", 
                      stat = "identity", colour = "black") +
    ggplot2::geom_text(ggplot2::aes_string(label = "count"), 
                       size = 3, vjust = -0.5, hjust = -0.5, angle = 35,
                       position = ggplot2::position_dodge(width = 1)) + 
    ggplot2::labs(title = paste("Distribution:", biotype)) + 
    ggplot2::scale_y_continuous(name = "Number of Genes", 
                                limits = c(0, 3000), breaks = round(seq(min(0), max(3000), by = 200), 1)) +
    ggplot2::scale_fill_manual(values = myPalette, name = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 12), 
                   plot.title = ggplot2::element_text(size = 20, face = "bold", 
                                                      hjust = 0.5), axis.title.x = ggplot2::element_blank(), 
                   axis.text.x = ggplot2::element_text(size = 12, angle = 35, 
                                                       hjust = 1, face = "bold"), axis.title.y = ggplot2::element_text(size = 12, 
                                                                                                                       face = "bold"), axis.text.y = ggplot2::element_text(size = 12), 
                   panel.border = ggplot2::element_rect(colour = "BLACK", 
                                                        size = 0.5), panel.grid.major = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank(), strip.background = ggplot2::element_rect())
  ggsave(paste("tau_tissue_specificity_GTeX_numberGenes_", biotype, "_noReproductive.png", sep = ""), width = 22, height = 8)
  
}
plotDistBiotype(tauAnno, "protein_coding")
plotDistBiotype(tauAnno, "long_noncoding")
plotDistBiotype(tauAnno, "short_noncoding")
plotDistBiotype(tauAnno, "pseudogene")
plotDistBiotype(tauAnno, "TEC")

#plot distribution - number of genes at different cutoffs ---------
tau1.00 <- getDist(tauAnno, 1)[, "total"]
tau0.99 <- getDist(tauAnno, 0.99)[, "total"]
tau0.95 <- getDist(tauAnno, 0.95)[, "total"]
tau0.90 <- getDist(tauAnno, 0.90)[, "total"]
dist <- data.frame(cbind(tau1.00, tau0.99, tau0.95, tau0.90))
tissues <- names(tauAnno[, c(-1, -2, -3)])
dist$tissues <- tissues
dist <- tidyr::gather(dist, "spec", "count", -tissues)
ggplot2::ggplot(dist, ggplot2::aes_string(x = "tissues", 
                                          y = "count", fill = "spec")) +
  ggplot2::geom_bar(position = "dodge", 
                    stat = "identity", colour = "black") +
  ggplot2::geom_text(ggplot2::aes_string(label = "count"), 
                     size = 3, vjust = -0.5, hjust = -0.5, angle = 35,
                     position = ggplot2::position_dodge(width = 1)) + 
  ggplot2::labs(title = "Distribution") + 
  ggplot2::scale_y_continuous(name = "Number of Genes", 
                              limits = c(0, 1000), breaks = round(seq(min(0), max(1000), by = 100), 1)) +
  ggplot2::scale_fill_manual(values = myPalette, name = "") +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 12), 
                 plot.title = ggplot2::element_text(size = 20, face = "bold", 
                                                    hjust = 0.5), axis.title.x = ggplot2::element_blank(), 
                 axis.text.x = ggplot2::element_text(size = 12, angle = 35, 
                                                     hjust = 1, face = "bold"), axis.title.y = ggplot2::element_text(size = 12, 
                                                                                                                     face = "bold"), axis.text.y = ggplot2::element_text(size = 12), 
                 panel.border = ggplot2::element_rect(colour = "BLACK", 
                                                      size = 0.5), panel.grid.major = ggplot2::element_blank(), 
                 panel.grid.minor = ggplot2::element_blank(), strip.background = ggplot2::element_rect())
ggsave("tau_tissue_specificity_GTeX_numberGenes_noReproductive.png", width = 22, height = 8)

#focus on brain
tissueA <- getTissue('Brain', qnExp, tauAnno)
head(tissueA[, -2], n = 5)

optimum <- getOptimum(tauAnno, tissueA, 10)
optimum$dataframe[, -2]

optimum$barplot

meanExp[c(rownames(optimum$dataframe)),]

corrPlots <- plotCorr(tissueA, c("CDR1", "AVP", "HCRT",
                                 "TPH2", "BARHL1", "TMEM88B",
                                 "GPR101", "OTP", "GABRA6",
                                 "MC3R", "F8A2", "MEPE"))
#corrPlots$tauPlot
corrPlots$fracPlot
ggsave("HPA_top12_vs_GTeX_tau_noReproductive.png", width = 8, height = 8)

plotGene(tauAnno, 'TPH2')
plotGene(tauAnno, 'CDR1')
plotGene(tauAnno, 'MOG')
plotGene(tauAnno, 'POPDC2')
plotGene(tauAnno, 'MOBP')
plotGene(tauAnno, 'MYL7')
plotGene(tauAnno, 'NPPA')
plotGene(tauAnno, 'NPPB')


id <- rownames(subset(tauAnno, tauAnno$external_gene_name == 'TPH2'))
round(meanExp[id, ], digits = 3)


#output list of tissue specific genes
myTauAnno <- right_join(gencode.anno, rownames_to_column(tauExp, var = "gene_id"))
write_csv(myTauAnno, "tau_tissue_specificity_gtex_noReproductive.csv")

tau.all <- read_csv("tau_tissue_specificity_gtex.csv")
tau.norepro <- read_csv("tau_tissue_specificity_gtex_noReproductive.csv")

plot(tau.all$Brain, tau.norepro$Brain)
