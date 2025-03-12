##HEADER------------------------------------------------------------------------

##SCRIPT FINAL-CMV MAIN FIGURES##
#Authors: Gusztav Milotay, Martin Little, Robert Watson & Benjamin Fairfax
#Last update: 15/01/2025
#Fully annotated: Yes
#Reviewed: Yes

##START-------------------------------------------------------------------------
##LOAD PACKAGES-----------------------------------------------------------------

#install.packages("survminer")
library(tidyverse)
library(ggpubr)
library(rstatix)
library(survival)
library(survminer)
library(DESeq2)
library(limma)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(org.Hs.eg.db) 
library(AnnotationDbi)
library(gplots)
library(ashr)
#library(XGR)
library(decoupleR)
library(grid)
library(broom)
library(grid)
library(gt)
library(patchwork)
library(tidycmprsk)
library(ggsurvfit)
library(ggbeeswarm)
library(gridExtra)
library(cowplot)
library(pROC)
library(MatchIt)
library(coin)
library(ggridges)
library(data.table)



##LOAD DATA---------------------------------------------------------------------

#bulk TCR data

tcr_processed <- readRDS("/data/TCRdata.rds")

#CMV data 

cmv <- fread("/data/CMVdata.txt", sep = "\t")

cmv <- cmv %>% 
  separate(`PtName ResultTrans`, c("PtName", "ResultTrans"), sep = " ")

#Clinical data

clinicalData <- read.csv("/data/clinicalData.csv")

#Bulk RNA-seq counts data 

counts <- readRDS("/data/BulkData.rds")

#Bulk RNA-seq batch data

batch_data <- read.csv("/data/BatchData.csv")

#T-cell flow data

flow <- readRDS("/data/FlowData.rds")

#Toxicity metadata 

tox <- read.csv("/data/AIpatient411complete.csv")

#Haematological data 

bloods <- read.csv("/data/HaemData.csv")

#UK Biobank data 

Biobank <- readRDS("/data/UKB.Rds")

#Single cell data for CD8+ T cells  

cd8 <- readRDS("/data/scData.rds")

#Gene list used for calculation of cytotoxicity score- same as used in Watson et al. 2021 

cytotoxicity <- read.csv("/data/BulkScores.csv")

clinicalData <- clinicalData %>%
  mutate(PtName = as.character(ID)) %>% 
  dplyr::select(-X)

#Corresponding CMV serology for most up to date patient list 

cmv$PtName <- gsub(c("R|PM|CS"), "", cmv$PtName)
clinicalData <- left_join(clinicalData, cmv)
clinicalData <- clinicalData %>% 
  mutate(cmv = ifelse(ResultTrans == "EQUIVOCAL", NA, 
                      ifelse(ResultTrans == "DETECTED", "CMV+",
                             ifelse(ResultTrans == "NOT_DETECTED", "CMV-", NA)))) 

survival_mm <- clinicalData %>% 
  filter(cancer == "Melanoma") %>% 
  filter(Rx != "AdjPembro" & Rx != "AdjNivo") 

#Metastatic melanoma dataset 

mmdata <- survival_mm %>% 
  mutate(patient = ID) %>% 
  mutate(C2 = "C2", C1 = "C1") %>% 
  pivot_longer(names_to = "cycle", cols = c(C2, C1), values_to = "cycle2") %>% 
  dplyr::select(-cycle2) %>% 
  mutate(sample = paste(patient, cycle, sep = "_"))
  
#Set theme 

theme_set(theme_classic()+
            theme(axis.ticks.length = unit(0.15, "cm"),
                  axis.text.x = element_text(angle=45, hjust=1, size=14),
                  axis.text.y = element_text(size=14),
                  axis.line = element_line(colour = "grey30", size = 0.2),
                  axis.ticks = element_line(color= "grey30", size=0.2), 
                  strip.background =element_rect(fill= "grey95"),
                  strip.text = element_text(colour = "grey30"),
                  axis.title.x = element_text(size=18),
                  axis.title.y = element_text(size=18),
                  strip.text.x = element_text(size = 18), 
                  strip.text.y = element_text(size = 18)))


##FUNCTIONS---------------------------------------------------------------------

#Extract normalised counts from DESEQ

baselineNormaliseTransform<-function(files,samples){
  allddsHTSeq<-DESeqDataSetFromMatrix(countData = files,
                                      colData = samples,
                                      design=~1)#No design
  allddsHTSeq<-estimateSizeFactors(allddsHTSeq)
  allddsHTSeq<-estimateDispersions(allddsHTSeq)
  data  <- counts(allddsHTSeq, normalized = TRUE)
  vsd<-getVarianceStabilizedData(allddsHTSeq)
  data <- as.data.frame(vsd)
  data$X <- rownames(data)
  
  return(norm_counts = data)
}

#Geometric mean calculator 

GeometricMean <- function(x) {
  x <- x[x > 0]  
  if (length(x) == 0) {
    return(NA) 
  }
  exp(mean(log(x)))  
}

################################################################################
##FIGURE I- BASELINE CONSEQUENCES OF CMV INFECTION##
################################################################################

bloods <- bloods %>% 
  mutate(patient = ID)
BloodsCMV <- left_join(survival_mm, bloods)

################################################################################
#FIGURE 1a- Baseline Lymphocyte count according to CMV status 
# 
# BloodsCMV %>% 
#   filter(!is.na(cmv)) %>% 
#   filter(!is.na(Lymph_pre)) %>% 
#   dplyr::mutate(sex = ifelse(sex == 1, "Male", "Female")) %>% 
#   mutate(across(where(is.character), as.factor)) %>% 
#   dplyr::select(age, sex, BRAF_status, cmv) %>% 
#   gtsummary::tbl_summary(by = cmv) %>% 
#   gtsummary::add_p() %>% 
#   modify_spanning_header(all_stat_cols() ~ "**Baseline Haematological data**")

#Check baseline lymphocyte count in CMV

#total 229, 124 CMV- and 105 CMV+
#Median CMV+ = 1.81
#Median CMV- = 1.48
# data:  Lymph_pre by cmv
# W = 4554.5, p-value = 9.09e-05
# alternative hypothesis: true location shift is not equal to 0
# 95 percent confidence interval:
#   -0.4599776 -0.1500544
# sample estimates:
#   difference in location 
# -0.300001

x <- ggplot(subset(BloodsCMV, !is.na(cmv)), aes(x = cmv, y = Lymph_pre))+
  ggbeeswarm::geom_quasirandom(aes(fill=cmv),
                               colour='grey30',
                               size=1.5,
                               pch=21,
                               alpha=0.6,
                               width=0.075)+ 
  geom_boxplot(colour = "grey30", aes(fill = cmv), alpha = 0.5, width = 0.5)+
  scale_fill_brewer(palette = 1)+
  theme(legend.position = "none")+
  xlab("CMV status")+
  ylab(paste("Baseline Lymphocyte count"))+
  # stat_compare_means(method = "wilcox", label.x.npc = "center", comparisons = list(c("CMV+", "CMV-")), size = 5)
  stat_compare_means(comparisons = list(c("CMV+", "CMV-")), label = "p.signif", size = 7)
  ggsave("/results/Fig1a.pdf", plot = x, width = 5, height = 8)

################################################################################
#FIGURE 1b- NLR according to CMV status 

  #229 patients, 124 CMV- and 105 CMV+
  #Median CMV+ = 2.59
  #Median CMV- = 3.26
  # data:  NLR_pre by cmv
  # W = 8092, p-value = 0.001546
  # alternative hypothesis: true location shift is not equal to 0
  # 95 percent confidence interval:
  #   0.2640733 1.0960376
  # sample estimates:
  #   difference in location 
  # 0.6783185

#Look at NLR which is associated with negative prognosis on immunotherapy 

x <- ggplot(subset(BloodsCMV, !is.na(cmv)), aes(x = cmv, y = NLR_pre))+
  ggbeeswarm::geom_quasirandom(aes(fill = cmv),
                               colour='grey30',
                               size=1.5,
                               pch=21,
                               alpha=0.6,
                               width=0.075)+
  geom_boxplot(colour = "grey30", aes(fill = cmv), alpha = 0.5, width = 0.5)+
  scale_fill_brewer(palette = 1)+
  theme(legend.position = "none")+
  xlab("CMV status")+
  ylab(paste("Baseline Neutrophil/Lymphocyte ratio"))+
  #stat_compare_means(method = "wilcox", label.x.npc = "center", comparisons = list(c("CMV+", "CMV-")), size = 5)
  stat_compare_means(comparisons = list(c("CMV+", "CMV-")), label = "p.signif", size = 7)
ggsave("/results/Fig1b.pdf", plot = x, width = 5, height = 8)

################################################################################

#Create a sample column in the flow data that will allow for merging with the metadata

flow <- flow %>% 
  mutate(ID = as.numeric(ID)) %>% 
  mutate(patient = ID) %>% 
  dplyr::select(-ID) %>% 
  mutate(sample = paste(patient, cycle, sep = "_"))

#Join flow data with metadata 

FlowMeta <- left_join(flow, mmdata)

#Identify major subsets of interest in CD8+ T cells

SubsetsCD8 <- c("TEM_CD8pos_CD3pos", "TNaive_CD8pos_CD3pos", "TCM_CD8pos_CD3pos", "TEMRA_CD8pos_CD3pos")
FlowMetaSubsetCD8 <- FlowMeta %>% 
  filter(cycle == "C1") %>% 
  filter(!is.na(cancer)) %>% 
  filter(!is.na(cmv)) %>% 
  filter(metric %in% SubsetsCD8) %>% 
  mutate(metric = ifelse(metric == "TNaive_CD8pos_CD3pos", "T Naive",
                         ifelse(metric == "TCM_CD8pos_CD3pos", "TCM",
                                ifelse(metric == "TEMRA_CD8pos_CD3pos", "TEMRA", "TEM"))))

FlowMetaSubsetCD8$metric <- factor(FlowMetaSubsetCD8$metric, levels = c("T Naive", "TCM", "TEM", "TEMRA"))

#74 patients 

#Repeat for CD4+ T cells 

SubsetsCD4 <- c("TEM_CD4pos_CD3pos", "TNaive_CD4pos_CD3pos", "TCM_CD4pos_CD3pos", "TEMRA_CD4pos_CD3pos")
FlowMetaSubsetCD4 <- FlowMeta %>% 
  filter(cycle == "C1") %>% 
  filter(!is.na(cancer)) %>% 
  filter(!is.na(cmv)) %>% 
  filter(metric %in% SubsetsCD4) %>% 
  mutate(metric = ifelse(metric == "TNaive_CD4pos_CD3pos", "T Naive",
                         ifelse(metric == "TCM_CD4pos_CD3pos", "TCM",
                                ifelse(metric == "TEMRA_CD4pos_CD3pos", "TEMRA", "TEM"))))

FlowMetaSubsetCD4$metric <- factor(FlowMetaSubsetCD4$metric, levels = c("T Naive", "TCM", "TEM", "TEMRA"))
  
################################################################################
#FIGURE 1c- CD4+ T cell subset changes in CMV 

# FlowMetaSubsetCD4 %>% 
#   filter(!is.na(cmv)) %>% 
#   dplyr::mutate(sex = ifelse(sex == 1, "Male", "Female")) %>% 
#   mutate(across(where(is.character), as.factor)) %>% 
#   dplyr::select(age, sex, BRAF_status, cmv, patient) %>% 
#   unique() %>% 
#   dplyr::select(-patient) %>% 
#   gtsummary::tbl_summary(by = cmv) %>% 
#   gtsummary::add_p() %>% 
#   modify_spanning_header(all_stat_cols() ~ "**Baseline Flow cytometry data**")

x <- ggplot(FlowMetaSubsetCD4, aes(x = cmv, y = value))+
  ggbeeswarm::geom_quasirandom(aes(fill= cmv),
                               colour='grey30',
                               size=1.5,
                               pch=21,
                               alpha=0.6,
                               width=0.075)+   
  geom_boxplot(colour = "grey30", aes(fill = cmv), alpha=0.5, width = 0.5)+
  scale_fill_brewer(palette = 2)+
  facet_wrap(~metric, scales = "free", nrow = 1)+
  theme(legend.position = "none")+
  xlab("CMV status")+
  ylab("Baseline proportion of CD4+ T cells (%)")+
  #stat_compare_means(method = "wilcox", label.x.npc = "center", comparisons = list(c("CMV+", "CMV-")), size = 5)
  stat_compare_means(comparisons = list(c("CMV+", "CMV-")), label = "p.signif", size = 7)
  ggsave("/results/Fig1c.pdf", plot = x, width = 12, height = 8)

################################################################################
#Figure 1d- CD8+ T cells subset changed in CMV 

x <- ggplot(FlowMetaSubsetCD8, aes(x = cmv, y = value))+
  ggbeeswarm::geom_quasirandom(aes(fill=cmv),
                               colour='grey30',
                               size=1.5,
                               pch=21,
                               alpha=0.6,
                               width=0.075)+   
  geom_boxplot(colour = "grey30", aes(fill = cmv), alpha=0.5, width = 0.5)+
  scale_fill_brewer(palette = 2)+
  facet_wrap(~metric, scales = "free", nrow = 1)+
  xlab("CMV status")+
  ylab("Baseline proportion of CD8+ T cells (%)")+
  theme(legend.position = "none")+
  #stat_compare_means(method = "wilcox", label.x.npc = "center", comparisons = list(c("CMV+", "CMV-")), size = 5)
  stat_compare_means(comparisons = list(c("CMV+", "CMV-")), label = "p.signif", size = 7)
ggsave("/results/Fig1d.pdf", plot = x, width = 12, height = 8)

################################################################################
#Figure 1e- CD25+FOXP3+ Treg subsets in CMV 

FlowTreg <- c("TregCD25_CD4pos_CD3pos")

FlowTregSubsets <- FlowMeta %>%
  filter(cycle == "C1") %>%
  filter(!is.na(cancer)) %>%
  filter(!is.na(cmv)) %>%
  filter(metric %in% FlowTreg) %>%
  mutate(metric = ifelse(metric == "TregCD25_CD4pos_CD3pos", "CD25+ Treg", NA))

x <- ggplot(subset(FlowTregSubsets, metric == "CD25+ Treg"), aes(x = cmv, y = value))+
  ggbeeswarm::geom_quasirandom(aes(fill=cmv),
                               colour='grey30',
                               size=1.5,
                               pch=21,
                               alpha=0.6,
                               width=0.075)+   
  geom_boxplot(colour = "grey30", aes(fill = cmv), alpha=0.5, width = 0.5)+
  scale_fill_brewer(palette = 2)+
  xlab("CMV status")+
  ylab("Baseline proportion of CD25+ Tregs (%)")+
  theme(legend.position = "none")+
  # stat_compare_means(method = "wilcox", label.x.npc = "center", comparisons = list(c("CMV+", "CMV-")), size = 5)
  stat_compare_means(comparisons = list(c("CMV+", "CMV-")), label = "p.signif", size = 7)
  ggsave("/results/Fig1e.pdf", plot = x, width = 5, height = 8)


  # data:  value by cmv
  # W = 986.5, p-value = 0.000763
  # alternative hypothesis: true location shift is not equal to 0
  # 95 percent confidence interval:
  #   0.4000571 1.4699570
  # sample estimates:
  #   difference in location 
  # 0.9399625 

################################################################################

#Filter metastatic melanoma data for C1 expression and create a new sample variable 
#Baseline Ipi sample removed and patient 25 who was never treated 

metadata <- mmdata %>% 
  filter(patient != 25 & Rx != "Ipilimumab") %>% 
  filter(!is.na(cmv)) %>% 
  filter(cycle == "C1") %>% 
  mutate(sample = paste0(sample, "_CD8_None")) %>% 
  mutate(sample = as.factor(sample)) %>% 
  filter(!is.na(Rx)) %>% 
  mutate(combination = ifelse(Rx == "IpiNivo", "cICB", "sICB")) 

#For the batch_data, dplyr::select CD8 T cell samples taken at baseline (format C1_CD8_)

batch <- batch_data %>% 
  mutate(patient = as.numeric(patient)) %>% 
  mutate(sample = as.factor(sample)) %>% 
  filter(grepl(paste("C1_CD8_None", sep = "|"), sample))

#Merge batch data with metadata which can be used in DESeq2 for batch correction 
#merge using patient, cycle, and sample and discard anything that does not match (unmatched = "drop")

MetaBatch <- left_join(metadata, batch, by = c("patient", "cycle", "sample") , unmatched = "drop", keep = F)

NewSamples <- c(262,264,271, 274, 275, 276, 279, 281, 286, 291, 297, 302,265, 269, 273, 284, 288, 292)
MissingBatches2022.Sep <- c(233, 239, 243, 254, 261)
MissingBatchesMay2021 <- c(224)
MissingBatchesDec2020 <- c(161)
MissingBeforeSep2019 <- c(2, 26, 41, 44, 45, 81, 86)
MissingSep2019 <- c(153, 155)

FilteredMeta <- MetaBatch %>%
  mutate_if(is.character, as.factor) %>% 
  mutate(cmv = ifelse(cmv == "CMV+", "CMVpositive",
                      ifelse(cmv == "CMV-", "CMVnegative", NA))) %>%
  mutate(combination = ifelse(Rx == "IpiNivo", "cICB", "sICB")) %>% 
  mutate(origBatch = case_when(
    patient %in% NewSamples ~ "Jan2024",
    TRUE ~ origBatch)) %>% 
  mutate(origBatch = case_when(patient %in% MissingBatches2022.Sep ~ "2022.Sep", TRUE ~ origBatch)) %>% 
  mutate(origBatch = case_when(patient %in% MissingBatchesMay2021 ~ "May_2021", TRUE ~ origBatch)) %>% 
  mutate(origBatch = case_when(patient %in% MissingBatchesDec2020 ~ "Dec_2020", TRUE ~ origBatch)) %>% 
  mutate(origBatch = case_when(patient %in% MissingBeforeSep2019 ~ "Before_SEP2019", TRUE ~ origBatch)) %>% 
  mutate(origBatch = case_when(patient %in% MissingSep2019 ~ "Sep_19", TRUE ~ origBatch)) %>% 
  mutate(origBatch = ifelse(patient > 302, "June_2024", origBatch)) 

#Ensure the rownames still correspond to the sample ids

rownames(FilteredMeta) <- FilteredMeta$sample

#Subset counts data for baseline CD8+ T cell samples (C1_CD8_None)
#Add Ensembl ids as rownames in this new dataset

#rownames(counts) <- counts$Gene
SubsetCounts <- counts %>%
  dplyr::select(matches(c("C1_CD8_None")))
rownames(SubsetCounts) <- rownames(counts)

#Metadata and counts data do not have matching samples. Match them and subset the matches 
#Ensure that the counts column names are in the rownames of the metadata
#Find where the two dataframes intersect and filter according to the matches

all(colnames(SubsetCounts) %in% rownames(FilteredMeta))
matching <- intersect(colnames(SubsetCounts), rownames(FilteredMeta))
FilteredMeta <- subset(FilteredMeta, (sample %in% matching))
rownames(FilteredMeta) <- FilteredMeta$sample
SubsetCounts <- SubsetCounts %>% 
  dplyr::select(contains(matching))
rownames(SubsetCounts) <- rownames(counts)

#Not only do all the samples have to match in the metadata and the counts, but they have to be in the same order. Reorder the counts data in accordance with the row order in the metadata 
#dplyr::select metadata columns of interest: sample, sex, cycle, batch, age as continuous variable and cmv status- remove columns that are not required 

CountsReordered <- SubsetCounts[, match(rownames(FilteredMeta), colnames(SubsetCounts))]
all(colnames(CountsReordered) == rownames(FilteredMeta))
rownames(CountsReordered) = rownames(counts) 
CountsReordered$id <- rownames(CountsReordered)
FinalMeta <- FilteredMeta %>% 
  dplyr::select(sample, patient, sex, age, cmv, origBatch, Rx, BRAF_status) %>% 
  mutate(sex = as.factor(sex)) %>% 
  mutate(Age = scale(age, center = T))
rownames(FinalMeta) <- FinalMeta$sample

#Check that columns and rows match again 

all(colnames(CountsReordered) %in% rownames(FinalMeta))
matching <- intersect(rownames(FinalMeta), colnames(CountsReordered))
CountsReordered2 <- CountsReordered %>% 
  dplyr::select(all_of(matching))
CountsReordered2$id <- CountsReordered$id

#Filter out genes which have low reads as these could easily be classed as differentially expressed. Vignette proposes using a count of 10 as a cutoff 
#Do a final check that the samples in the counts data are matching and in the same order as the samples in the metadata 

keep <- rowMeans(CountsReordered2[,1:(ncol(CountsReordered2)-1)]) >= 10
CountsReordered2 <- CountsReordered2[keep, ] 
CountsFinal <- dplyr::select(CountsReordered2, -id)
rownames(CountsFinal) <- CountsReordered2$id
all(colnames(CountsFinal) == rownames(FinalMeta))

#Create a DESEQ object with age, sex and batch as covariates  

#206 samples from metastatic melanoma patients with CMV data
#95 CMV+ samples
#111 CMV- samples

dds <- DESeqDataSetFromMatrix(countData = CountsFinal, 
                              colData = FinalMeta, 
                              design = ~ Age  + sex + origBatch + cmv)

#Set reference CMV to CMV- 

dds$cmv <- relevel(dds$cmv, ref = "CMVnegative")

#Run DESeq

dds <- DESeq(dds)

#Extract results for CMV 

ResultsDDS <- results(dds, name = "cmv_CMVpositive_vs_CMVnegative", alpha = 0.05)

summary(ResultsDDS)

#Reorder the genes by significance 

ResultsOrdered <- ResultsDDS[order(ResultsDDS$pvalue),]

#Save results in a dataframe 

ResultsOrdered <- as.data.frame(ResultsOrdered)

#Map ENSEMBL ids 

ResultsOrdered$symbol <- mapIds(org.Hs.eg.db, keys = rownames(ResultsOrdered), keytype = "ENSEMBL", column = "SYMBOL")

#Label significantly differentially expressed genes (according to the adjusted P-value) as being up or down-regulated depending on the log2fc 

ResultsOrdered$Response <- "None"
ResultsOrdered$Response[ResultsOrdered$log2FoldChange > 0 & ResultsOrdered$padj < 0.05] <- "CMV+"
ResultsOrdered$Response[ResultsOrdered$log2FoldChange < 0 & ResultsOrdered$padj < 0.05] <- "CMV-"

################################################################################
#TABLE 1- Differentially expressed genes in CD8+ T cells in CMV seropositivity 

Table1 <- ResultsOrdered %>% 
  filter(padj < 0.05)

#7576 genes out of 17,308 were differentially expressed
#19% upregulated (3224) and 25% downregulated (4352)

################################################################################

#Find top 25 up and down-regulated genes and combine them in a dataframe of the top 50 DE genes

Downregulated <- ResultsOrdered %>% 
  subset(Response == "CMV-") %>% 
  arrange(padj) %>% 
  head(n = 25)
Upregulated <- ResultsOrdered %>% 
  subset(Response == "CMV+") %>% 
  arrange(padj) %>% 
  head(n = 25)
Labels <- rbind(Upregulated, Downregulated)

################################################################################
#FIGURE 1f- Differentially expressed genes by CMV in CD8+ T cells 
# 
# FinalMeta %>% 
#   mutate(CMV = ifelse(cmv == "CMVpositive", "CMV+", "CMV-")) %>% 
#   dplyr::mutate(sex = ifelse(sex == 1, "Male", "Female")) %>% 
#   mutate(across(where(is.character), as.factor)) %>% 
#   dplyr::select(age, sex, CMV, origBatch) %>% 
#   gtsummary::tbl_summary(by = CMV) %>% 
#   gtsummary::add_p() %>% 
#   modify_spanning_header(all_stat_cols() ~ "**Baseline Bulk CD8+ RNA-seq data**")

ylab <- expression(-log[10]~ (adjusted ~ P-value))
xlab <- expression(log[2]~FC )
x <- ggplot(ResultsOrdered,
       aes(x = log2FoldChange,
           y = -log10(padj),
           col = Response,
           label = symbol))+
  geom_point(alpha = 0.2, size = 2)+
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", alpha = 0.5)+
  geom_vline(xintercept = c(1,-1), linetype = "dotted", alpha = 0.5)+
  xlab(xlab)+
  ylab(ylab)+
  theme_classic()+
  theme(axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(hjust=1, size=14),
        axis.text.y = element_text(size=14),
        axis.line = element_line(colour = "grey30", size = 0.2),
        axis.ticks = element_line(color= "grey30", size=0.2), 
        strip.background =element_rect(fill= "grey95"),
        strip.text = element_text(colour = "grey30"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        strip.text.x = element_text(size = 18), 
        strip.text.y = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18))+
  geom_text_repel(data = Labels, colour = "grey30", size = 3.75, max.overlaps = 10, alpha = 0.6)+
  scale_color_manual(values=c("#2171B5", "#FFBB44", 'grey80'))

ggsave("/results/Fig1f.pdf", plot = x, width = 8, height = 8)

################################################################################

#Loading ontology from server for XGR pathway analysis 

RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
GOBP <- xRDataLoader(RData = "org.Hs.egGOBP", RData.location = RData.location)

#Obtain the background which will be all genes- filtering required to obtain genes with symbol annotations

ResultsFiltered <- ResultsOrdered %>%
  filter(!is.na(symbol))
Background <- ResultsFiltered$symbol

#14614 genes used as background

#Filter for genes which are significantly upregulated in CMV

CMVseropositive <- ResultsOrdered %>%
  filter(log2FoldChange > 0 & padj < 0.05)
CMVseropositive <- CMVseropositive$symbol

#Filter for genes which are downregulated in CMV

CMVseronegative <- ResultsOrdered %>%
  filter(log2FoldChange < 0 & padj < 0.05)
CMVseronegative <- CMVseronegative$symbol

#Run GSEA using GOBP ontology for up- and downregulated genes separately

Induced <- xEnricherGenes(data = CMVseropositive,
                          background = Background,
                          ontology = "GOBP",
                          test = "hypergeo",
                          min.overlap = 10,
                          elim.pvalue  = 0.01,
                          ontology.algorithm = "elim")
Table2Induced <- as.data.frame(xEnrichViewer(Induced, top_num = 1000))
Table2Induced$direction <- "Induced"
Induced <- as.data.frame(xEnrichViewer(Induced, top_num = 25))

Downregulated <- xEnricherGenes(data = CMVseronegative,
                                background = Background,
                                ontology = "GOBP",
                                test = "hypergeo",
                                min.overlap = 10,
                                elim.pvalue  = 0.01,
                                ontology.algorithm = "elim")
Table2Downregulated <- as.data.frame(xEnrichViewer(Downregulated, top_num = 1000))
Table2Downregulated$direction <- "Downregulated"
Downregulated <- as.data.frame(xEnrichViewer(Downregulated, top_num = 25))

################################################################################
#TABLE 2- Pathways induced and downregulated in CD8+ T cells in CMV seropositivity

Table2 <- rbind(Table2Induced, Table2Downregulated)
Table2 <- Table2 %>%
  filter(adjp < 0.05)

#181 differentially regulated pathways
#141 induced and 40 downregulated

#TCR signalling 8.2e-8
#Interferon gamma: 1.8e-5
#Antigen presentation: 1.1e-5

################################################################################

#Indicate which genes have been induced and which have been suppressed

Induced$direction <- "induced"
Downregulated$direction <- "suppressed"

#Bind the results dataframes and multiple fc by -1 for suppressed pathways

Pathways <- rbind(Induced, Downregulated)
Pathways <- Pathways %>%
  mutate(fc = ifelse(direction == "suppressed", fc* (-1), fc)) %>%
  filter(fc > 1.8 | fc < -1.8) %>%
  filter(or > 2.3) %>%
  filter(adjp <= 1.3e-3)

#Create a new dataframe to displace lables from points

Pathways$GOBP <- row.names(Pathways)
Pathways$GOBP <- factor(Pathways$GOBP, levels = c(row.names(Pathways[order(Pathways$fc), ])))
Pathways <- Pathways %>%
  mutate(`adjP` = adjp) %>%
  mutate(`Overlap (n)` = nOverlap)
Naming <- Pathways[Pathways$adjp < 0.05, ]
Naming$fc <-
  ifelse(Naming$fc > 0,
         Naming$fc - 2.75,
         Naming$fc + 2)
Naming <- Naming %>%
  mutate(name = ifelse(name == "antigen processing and presentation of peptide antigen via MHC class I",
                       "MHC class I peptide presentation", name))

################################################################################
#FIGURE 1g- Differential pathway analysis according to CMV status in CD8+ T cells

ggplot(Pathways, aes(x = fc, y = GOBP)) +
  geom_point(aes(fill = -log10(`adjP`), size = `Overlap (n)`), pch = 21) +
  scale_fill_distiller(palette = "Spectral") +
  geom_text(data = Naming,
            aes(label = name, colour = direction),
            size = 6, show.legend = F) +
  xlab("Fold Change")+
  theme_classic()+
  theme(axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(hjust=1, size=12),
        axis.text.y = element_text(size=12),
        axis.line = element_line(colour = "grey30", size = 0.2),
        axis.ticks = element_line(color= "grey30", size=0.2),
        strip.background =element_rect(fill= "grey95"),
        strip.text = element_text(colour = "grey30"),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))+
  scale_color_manual(values = c("#EB7926", "#134B73"))

################################################################################
##FIGURE II- CMV INTERACTION WITH TREATMENT##
################################################################################

#Filter metadata for baseline and C2 samples- sample 260 was incorrectly assigned as a Melanoma 
#patient 25 never received treatment
#Analysis restricted to anti-PD-1 or anti-PD-1/anti-CTLA-4

metadata <- mmdata %>% 
  filter(patient != 25 & Rx != "Ipilimumab" & sample != "260_C2") %>% 
  filter(cycle == "C1" | cycle == "C2") %>% 
  mutate(sample = paste0(sample, "_CD8_None")) %>% 
  mutate(sample = as.factor(sample)) 

#For the batch_data, dplyr::select cells of interest in the format 'CX_CD8' referring to cycle and cell marker

batch <- batch_data %>% 
  mutate(patient = as.numeric(patient)) %>% 
  mutate(sample = as.factor(sample)) %>% 
  filter(grepl(paste("C1_CD8_None|C2_CD8_None", collapse = "|"), sample))

#Merge batch data with metadata which can be used in DESeq2 for batch correction 
#merge using patient, cycle, and sample and discard anything that does not match

MetadataBatch <- left_join(metadata, batch, by = c("patient", "cycle", "sample") , unmatched = "drop", keep = F)

FilteredMeta <- MetadataBatch %>% 
  mutate(over68 = ifelse(age >= 68, "68+", "<68")) %>% 
  mutate(cmv = ifelse(cmv == "CMV+", "CMVpositive",
                      ifelse(cmv == "CMV-", "CMVnegative", NA))) #this is because DESEQ2 does not like non-conventional characters 

#Ensure that sample names are rownames 

rownames(FilteredMeta) <- FilteredMeta$sample

#Subset counts data for cells of interest in the format '_CX_CD8_None'
#Add Ensembl ids as rownames in this new dataset- still unsure why R always removes rownames when you modify a dataframe

SubsetCounts <- counts %>% 
  dplyr::select(matches(c("C2_CD8_None", "C1_CD8_None")))
rownames(SubsetCounts) <- rownames(counts)

#Metadata and counts data do not have matching samples. Match them and subset the matches 
#see if the counts column names are in the rownames of the metadata
#Find where the two dataframes intersect and filter according to the matches

all(colnames(SubsetCounts) %in% rownames(FilteredMeta))
matching <- intersect(colnames(SubsetCounts), rownames(FilteredMeta))
FilteredMeta <- subset(FilteredMeta, (sample %in% matching))
rownames(FilteredMeta) <- FilteredMeta$sample
SubsetCounts <- SubsetCounts %>% 
  dplyr::select(contains(matching))
rownames(SubsetCounts) <- rownames(counts)

#Not only do all the samples have to match in the metadata and the counts, but they have to be in the same order. Reorder the counts data 
#dplyr::select metadata columns of interest: sample, sex, cycle, batch, age as continuous variable, cmv sereostatus, etc...

SubsetReordered <- SubsetCounts[, match(rownames(FilteredMeta), colnames(SubsetCounts))]
all(colnames(SubsetReordered) == rownames(FilteredMeta))
rownames(SubsetReordered) <- rownames(counts) 
SubsetReordered$id <- rownames(SubsetReordered)
FinalMeta <- FilteredMeta %>% 
dplyr::select(sample, sex, cycle, age, cmv, Rx, patient, BRAF_status)

#Filter out genes which have low reads as these could easily be classed as differentially expressed. 10 is the count threshold suggested by the vignette.
#Do a final check that the samples in the counts data are matching and in the same order as the samples in the metadata 

keep <- rowMeans(SubsetReordered[,1:(ncol(SubsetReordered)-1)]) >= 10
SubsetReordered <- SubsetReordered[keep, ] 
rownames(SubsetReordered) <- SubsetReordered$id
CountsFinal <- dplyr::select(SubsetReordered, -id)
rownames(CountsFinal) <- rownames(SubsetReordered)
rownames(FinalMeta) <- FinalMeta$sample
all(colnames(CountsFinal) == rownames(FinalMeta))

#Filter for DEGs by CMV that are below the optimal threshold of significance- see Extended Figure 2a to see ROC analysis for determining optimal threshold

CMVgenes <- ResultsOrdered %>% 
  filter(padj < 1.581139e-14)

#75 genes

NormCounts <- baselineNormaliseTransform(files = CountsFinal, samples = FinalMeta)
NormCounts$symbol <- mapIds(org.Hs.eg.db, keys = rownames(NormCounts), keytype = "ENSEMBL", column = "SYMBOL")
CMVinduced <- NormCounts %>% 
  dplyr::select(-X) %>% 
  filter(rownames(NormCounts) %in% rownames(CMVgenes)) %>% 
  dplyr::select(-c(symbol))
CMVinducedWide <- t(CMVinduced)

#Add genes to metadata 

rownames(CMVinducedWide) <- FinalMeta$sample
CMVinducedWide <- as.data.frame(CMVinducedWide)
CMVinducedWide$sample <- rownames(CMVinducedWide)
CountsMeta <- left_join(FinalMeta, CMVinducedWide)

#Perform PCA analysis to create a score for each individual on the basis of CMV-induced and downregulated genes- this will be PC1

GeneCounts <- CMVinducedWide %>% 
  dplyr::select(-sample)
PCA <- prcomp(x = GeneCounts)
PCAdata <- as.data.frame(PCA$x)
CMVscore <- PCAdata %>% 
  dplyr::select(PC1) 
CMVscore$sample <- rownames(CMVscore)
PCAmeta <- left_join(FinalMeta, CMVscore)

################################################################################
#FIGURE 2a- CMV score according to treatment

#Identify paired samples
#151 paired samples with serology 
#108 cICB patients (64 CMV- and 44 CMV+) and 43 sICB patients (15 CMV- and 28 CMV+) 

PlotData <- PCAmeta %>% 
  mutate(Treatment = ifelse(Rx == "IpiNivo", "cICB", "sICB")) %>% 
  group_by(patient) %>% 
  filter(n_distinct(cycle) > 1) %>% 
  ungroup() %>% 
  filter(!is.na(cmv)) %>% 
  mutate(cmv = ifelse(cmv == "CMVpositive", "CMV+", "CMV-")) %>% 
  mutate(cycle = ifelse(cycle == "C1", "Baseline", "Post cycle 1"))

# PlotData %>% 
#   dplyr::mutate(sex = ifelse(sex == 1, "Male", "Female")) %>% 
#   mutate(across(where(is.character), as.factor)) %>% 
#   filter(cycle == "Baseline") %>% 
#   dplyr::select(age, sex, BRAF_status, cmv, Treatment) %>% 
#   gtsummary::tbl_summary(by = cmv) %>% 
#   gtsummary::add_p() %>% 
#   modify_spanning_header(all_stat_cols() ~ "**Paired CD8+ Bulk RNA-seq data**")

x <- ggplot(PlotData, aes(x = cycle, y = PC1))+
  ggbeeswarm::geom_quasirandom(aes(fill = cycle),
                               colour='grey30',
                               size=1.5,
                               pch=21,
                               alpha=0.6,
                               width=0.075)+   
  geom_line(aes(group = patient), color = "grey", alpha = 0.5, linetype = "dashed") +
  geom_boxplot(colour = "grey30", aes(fill = cycle), alpha=0.5, width = 0.5)+
  scale_fill_brewer(palette = 7)+
  facet_grid(vars(Treatment), vars(cmv))+
  xlab("Cycle")+
  ylab("CMV score")+
  theme(legend.position = "none")+
  labs(fill = "Cycle", color = "Cycle")+
  #stat_compare_means(method = "wilcox", label.x.npc = "center", paired = T, comparisons = list(c("Baseline", "Post cycle 1")), size = 5)
  stat_compare_means(comparisons = list(c("Baseline", "Post cycle 1")), label = "p.signif", paired = T,  size = 6)

ggsave("/results/Fig2a.pdf", plot = x, width = 5, height = 12)

#CMV- cICB: p = 2.8e-9 
#CMV- sICB: p = 0.28 
#CMV+ cICB: p = 0.41 
#CMV+ sICB: p = 0.54 
  
################################################################################

#Explore Gene expression scores that were previously associated with survival- Cytotoxicity score  

metadata <- mmdata %>% 
  filter(patient != 25 & Rx != "Ipilimumab" & sample != "260_C2") %>% 
  filter(cycle == "C1" | cycle == "C2") %>% 
  mutate(sample = paste0(sample, "_CD8_None")) %>% 
  mutate(sample = as.factor(sample)) 

#For the batch_data 

batch <- batch_data %>% 
  mutate(patient = as.numeric(patient)) %>% 
  mutate(sample = as.factor(sample)) %>% 
  filter(grepl(paste("C1_CD8_None|C2_CD8_None", collapse = "|"), sample))

#Merge batch data with metadata which can be used in DESeq2 for batch correction 
#merge using patient, cycle, and sample and discard anything that does not match

MetaBatch <- left_join(metadata, batch, by = c("patient", "cycle", "sample") , unmatched = "drop", keep = F)
FilteredMeta <- MetaBatch %>% 
  mutate(cmv = ifelse(cmv == "CMV+", "CMVpositive",
                      ifelse(cmv == "CMV-", "CMVnegative", NA)))#this is because DESEQ2 does not like non-conventional characters 

#Ensure that sample names are rownames 

rownames(FilteredMeta) <- FilteredMeta$sample

#Subset counts data for baseline and C2 CD8+ T cells

SubsetCounts <- counts %>% 
  dplyr::select(matches(c("C2_CD8_None", "C1_CD8_None")))
#rownames(SubsetCounts) <- counts$Gene

#Metadata and counts data do not have matching samples. Match them and subset the matches 
#see if the counts column names are in the rownames of the metadata
#Find where the two dataframes intersect and filter according to the matches

all(colnames(SubsetCounts) %in% rownames(FilteredMeta))
Matching <- intersect(colnames(SubsetCounts), rownames(FilteredMeta))
FilteredMeta <- subset(FilteredMeta, (sample %in% Matching))
rownames(FilteredMeta) <- FilteredMeta$sample
SubsetCounts <- SubsetCounts %>% 
  dplyr::select(contains(Matching))
#rownames(SubsetCounts) <- counts$Gene

#Not only do all the samples have to match in the metadata and the counts, but they have to be in the same order. Reorder the counts data 
#dplyr::select metadata columns of interest

SubsetReordered <- SubsetCounts[, match(rownames(FilteredMeta), colnames(SubsetCounts))]
all(colnames(SubsetReordered) == rownames(FilteredMeta))
rownames(SubsetReordered) = rownames(SubsetCounts) 
SubsetReordered$id <- rownames(SubsetReordered)
FinalMeta <- FilteredMeta %>% 
  dplyr::select(sample, sex, cycle, age, cmv, Rx, sample, origBatch, patient, BRAF_status, cancer)

#Filter out genes which have low reads as these could easily be classed as differentially expressed. 10 is the count threshold suggested by the vignette.
#Do a final check that the samples in the counts data are matching and in the same order as the samples in the metadata 

keep <- rowMeans(SubsetReordered[,1:(ncol(SubsetReordered)-1)]) >= 10
SubsetReordered <- SubsetReordered[keep, ] 
rownames(SubsetReordered) <- SubsetReordered$id
CountsFinal <- dplyr::select(SubsetReordered, -id)
rownames(CountsFinal) <- rownames(SubsetReordered)
rownames(FinalMeta) <- FinalMeta$sample
all(colnames(CountsFinal) == rownames(FinalMeta))

#Perform count normalisation by passing counts and metadata arguments 

NormCounts <- baselineNormaliseTransform(files = CountsFinal, samples = FinalMeta)

#Ensure the rownames match the ENSEMBL IDs

NormCounts2 <- as.data.frame(NormCounts)

#Map the ENSEMBL ids to symbols using mapIds()

NormCounts2$symbol <- mapIds(org.Hs.eg.db, keys = rownames(NormCounts2), keytype = "ENSEMBL", column = "SYMBOL")

#Subset the counts to only include positively associated genes in the Cytotoxicity score  
#Merge normalised counts for genes of interest with FinalMeta

CytotoxicityGenes <- NormCounts2 %>% 
  dplyr::select(-X) %>% 
  filter(NormCounts2$symbol %in% cytotoxicity$Cytox.score) %>% 
  dplyr::select(-c(symbol))
CytotoxicityWide <- t(CytotoxicityGenes)
rownames(CytotoxicityWide) <- FinalMeta$sample
CytotoxicityWide <- as.data.frame(CytotoxicityWide)
CytotoxicityWide$sample <- rownames(CytotoxicityWide)
FinalMeta <- left_join(FinalMeta, CytotoxicityWide, by = "sample")

#Compute the geometric mean for genes associated with cytotoxicity 

ScoreMeta <- FinalMeta %>% 
  rowwise() %>%
  mutate(Cytotoxicity = GeometricMean(c_across(11:61))) %>% 
  ungroup()

################################################################################
#FIGURE 2b- Cytotoxicity score pre- and post-treatment according to treatment type and CMV status

#Identify paired samples 

PlotData <- ScoreMeta %>% 
  group_by(patient) %>% 
  filter(n_distinct(cycle) > 1) %>% 
  ungroup() %>% 
  filter(!is.na(cmv)) %>% 
  mutate(combination = ifelse(Rx == "IpiNivo", "cICB", "sICB")) %>% 
  mutate(cmv = ifelse(cmv == "CMVpositive", "CMV+", "CMV-")) %>% 
  mutate(cycle = ifelse(cycle == "C1", "Baseline", "Post cycle 1"))

#Plot 

x <- ggplot(subset(PlotData, combination == "cICB"), aes(x = cycle, y = Cytotoxicity))+
  ggbeeswarm::geom_quasirandom(aes(fill = cycle),
                               colour='grey30',
                               size=1.5,
                               pch=21,
                               alpha=0.6,
                               width=0.075)+   
  geom_line(aes(group = patient), color = "grey", alpha = 0.5, linetype = "dashed") +
  geom_boxplot(colour = "grey30", aes(fill = cycle), alpha=0.5, width = 0.5)+
  scale_fill_brewer(palette = 13)+
  facet_wrap(vars(cmv))+
  theme(legend.position = "None")+
  xlab("Cycle")+
  ylab("Cytotoxicity score")+
  #stat_compare_means(method = "wilcox", label.x.npc = "center", paired = T, comparisons = list(c("Baseline", "Post cycle 1")), size = 5)
  stat_compare_means(comparisons = list(c("Baseline", "Post cycle 1")), label = "p.signif", paired = T,  size = 6)

ggsave("/results/Fig2b.pdf", plot = x, width = 5, height = 8)

  #CMV- cICB: p = 2.1e-7 
  #CMV+ cICB: p = 0.0013

################################################################################

#Survival data

surv_last_update <- survival_mm %>% 
  filter(progression == 1 | months_progression > 6) %>% 
  mutate(months_death = as.numeric(months_death)) %>% 
  mutate(months_progression = as.numeric(months_progression)) %>% 
  mutate(months_death = ifelse(months_death > 90, 90, months_death)) %>% 
  mutate(months_progression = ifelse(months_progression > 90, 90, months_progression))

#302 Metastatic melanoma patients with 6 months followup 

surv_last_update <- surv_last_update %>% 
  mutate(combination = ifelse(Rx == "IpiNivo", "cICB", "sICB")) %>% 
  mutate(`CMV-Treatment` = paste(cmv, combination))

#Filter for anti-PD1 treated 
#105 patients 
#30 CMV- and 45 CMV+

sICB_surv <- surv_last_update %>% 
  filter(combination == "sICB") %>% 
  filter(Rx != "Ipilimumab" & Rx != "RelatNivo")

#Filter for anti-PD1/ anti-CTLA4 treated 
#194 patients 
#113 CMV- and 78 CMV+ 

cICB_surv <- surv_last_update %>% 
  filter(combination == "cICB") 

################################################################################
#FIGURE 2c- Kaplan Meier of CMV in cICB-treated cohort 

# cICB_surv %>%
#   dplyr::mutate(sex = ifelse(sex == 1, "Male", "Female")) %>%
#   mutate(across(where(is.character), as.factor)) %>%
#   dplyr::select(age, sex, BRAF_status, cmv) %>%
#   gtsummary::tbl_summary(by = cmv) %>%
#   gtsummary::add_p() %>%
#   modify_spanning_header(all_stat_cols() ~ "**cICB survival data**")

fit <- survfit(Surv(months_death, death_status) ~ cmv, data = cICB_surv)
summary(fit)$table
log_rank_test <- survdiff(Surv(months_death, death_status) ~ cmv, data = cICB_surv)
p_value <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)
cox <- coxph(Surv(months_death, death_status) ~ cmv, data = cICB_surv)
cICBsummary <- summary(cox)
coefcICB <- as.data.frame(cICBsummary$coef[, 1])
cicICB <- as.data.frame(cICBsummary$conf.int)
PlotcICB <- data.frame(
  covariate = c("CMV"),
  HR = cicICB[, "exp(coef)"],
  lower = cicICB[,"lower .95"],
  upper = cicICB[,"upper .95"],
  p_value = cICBsummary$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)

time_point_1yr <- 12
time_point_5yr <- 60
surv_prob_1yr <- summary(fit, times = time_point_1yr)$surv
surv_prob_5yr <- summary(fit, times = time_point_5yr)$surv
km <- ggsurvplot(fit,
               # pval = TRUE,  # Add log rank p-value
               # pval.method = TRUE,  # Add label for method used,
               pval.size=5,
               pval.coord = c(35, 0.70), pval.method.coord = c(35, 0.78),
               xlab = "Months", ylab = "Probability of overall survival", legend.title = "CMV status", legend.labs = c("CMV-", "CMV+"),
               font.main = 16,font.x =  16, font.y = 16, font.tickslab = 16, font.legend = 16, font.title = c(12, "bold"), break.x.by = 6, censor.shape = 124, panel.border = element_rect(colour = "black", fill=NA, size=4),
               risk.table = "nrisk_cumevents",risk.table.legend = FALSE, risk.table.y.text.col= FALSE, risk.table.font = 3.5, risk.table.order = c(0,1),
               ggtheme = theme_bw() + theme(panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(),
                                            plot.title = element_text(hjust = 0.5)),
               conf.int = F,
               palette = c("#2171B5","#FFBB44"),
               xlim = c(0, 90)
               )

km$plot<- km$plot +
  geom_segment(x = time_point_1yr, xend = time_point_1yr, y = -Inf, yend = 0.85 , color = "grey60", lty="dashed")+
  geom_segment(x = time_point_5yr, xend = time_point_5yr, y = -Inf, yend = 0.65, color = "grey60", lty="dashed")+
  annotate("text", x=time_point_1yr+5, y=surv_prob_1yr[1] -0.1, label=paste0(round(surv_prob_1yr[[1]],3)*100,"%"), size=5, colour="#2171B5")+
  annotate("text", x=time_point_1yr+5, y=surv_prob_1yr[2]+0.05, label=paste0(round(surv_prob_1yr[[2]],3)*100,"%"), size=5, colour="#FFBB44")+
  annotate("text", x=time_point_1yr, y=surv_prob_1yr[1]+0.2, label="1 year", size=5, colour="grey60")+
  annotate("text", x=time_point_5yr+5, y=surv_prob_5yr[1]+0.07, label=paste0(round(surv_prob_5yr[[1]],3)*100,"%"), size=5, colour="#2171B5")+
  annotate("text", x=time_point_5yr+5, y=surv_prob_5yr[2]-0.1, label=paste0(round(surv_prob_5yr[[2]],3)*100,"%"), size=5, colour="#FFBB44")+
  annotate("text", x=time_point_5yr, y=surv_prob_5yr[1]+0.2, label="5 years", size=5, colour="grey60")+
  annotate("text", x=time_point_5yr-20, y = 0.85, label=paste0("HR for death: ",
                                                               round(PlotcICB$HR, digits = 2),"\n (95% CI: ", round(PlotcICB$lower, digits = 2),
                                                               "-",round(PlotcICB$upper, digits = 2), ") \n Log-rank \n", round(p_value, digits = 2), sep = ""), size=5, fontface="plain")
km

ggsave("/results/Fig2c.pdf", plot = km$plot, width = 8, height = 8)

#p = 0.92
#HR = 1.02
#CMV+ median = 54 months
#CMV- median = 78.2 months 

################################################################################

#Prepare data for multivariate survival analysis 

surv_last_update <- surv_last_update %>% 
  mutate(Age = age) %>% 
  mutate(Sex = ifelse(sex == 2, "Female", "Male")) %>% 
  mutate(Treatment = Rx) %>% 
  mutate(CMV = cmv) %>% 
  mutate(`BRAF status` = BRAF_status) 

################################################################################

#FIGURE 2d- Kaplan Meier of CMV in sICB-treated cohort 

# sICB_surv %>%
#   dplyr::mutate(sex = ifelse(sex == 1, "Male", "Female")) %>%
#   mutate(across(where(is.character), as.factor)) %>%
#   dplyr::select(age, sex, BRAF_status, cmv) %>%
#   gtsummary::tbl_summary(by = cmv) %>%
#   gtsummary::add_p() %>%
#   modify_spanning_header(all_stat_cols() ~ "**sICB survival data**")


fit <- survfit(Surv(months_death, death_status) ~ cmv, data = sICB_surv)
summary(fit)$table
log_rank_test <- survdiff(Surv(months_death, death_status) ~ cmv, data = sICB_surv)
p_value <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)
cox <- coxph(Surv(months_death, death_status) ~ cmv, data = sICB_surv)
sICBsummary <- summary(cox)
coefsICB <- as.data.frame(sICBsummary$coef[, 1])
cisICB <- as.data.frame(sICBsummary$conf.int)
PlotsICB <- data.frame(
  covariate = c("CMV"),
  HR = cisICB[, "exp(coef)"],
  lower = cisICB[,"lower .95"],
  upper = cisICB[,"upper .95"],
  p_value = sICBsummary$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)

time_point_1yr <- 12
time_point_5yr <- 60
surv_prob_1yr <- summary(fit, times = time_point_1yr)$surv
surv_prob_5yr <- summary(fit, times = time_point_5yr)$surv
km <- ggsurvplot(fit,
                 # pval = TRUE,  # Add log rank p-value
                 # pval.method = TRUE,  # Add label for method used,
                 pval.size=5,
                 pval.coord = c(35, 0.70), pval.method.coord = c(35, 0.78),
                 xlab = "Months", ylab = "Probability of overall survival", legend.title = "CMV status", legend.labs = c("CMV-", "CMV+"),
                 font.main = 16,font.x =  16, font.y = 16, font.tickslab = 16, font.legend = 16, font.title = c(12, "bold"), break.x.by = 6, censor.shape = 124, panel.border = element_rect(colour = "black", fill=NA, size=4),
                 risk.table = "nrisk_cumevents",risk.table.legend = FALSE, risk.table.y.text.col= FALSE, risk.table.font = 3.5, risk.table.order = c(0,1),
                 ggtheme = theme_bw() + theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(),
                                              plot.title = element_text(hjust = 0.5)),
                 conf.int = F,
                 palette = c("#2171B5","#FFBB44"),
                 xlim = c(0, 90)
)

km$plot<- km$plot +
  geom_segment(x = time_point_1yr, xend = time_point_1yr, y = -Inf, yend = 0.85 , color = "grey60", lty="dashed")+
  geom_segment(x = time_point_5yr, xend = time_point_5yr, y = -Inf, yend = 0.70, color = "grey60", lty="dashed")+
  annotate("text", x=time_point_1yr+5, y=surv_prob_1yr[1] -0.15, label=paste0(round(surv_prob_1yr[[1]],3)*100,"%"), size=5, colour="#2171B5")+
  annotate("text", x=time_point_1yr+5, y=surv_prob_1yr[2]+0.05, label=paste0(round(surv_prob_1yr[[2]],3)*100,"%"), size=5, colour="#FFBB44")+
  annotate("text", x=time_point_1yr, y=surv_prob_1yr[1]+0.3, label="1 year", size=5, colour="grey60")+
  annotate("text", x=time_point_5yr+5, y=surv_prob_5yr[1]-0.1, label=paste0(round(surv_prob_5yr[[1]],3)*100,"%"), size=5, colour="#2171B5")+
  annotate("text", x=time_point_5yr+5, y=surv_prob_5yr[2]+0.05, label=paste0(round(surv_prob_5yr[[2]],3)*100,"%"), size=5, colour="#FFBB44")+
  annotate("text", x=time_point_5yr, y=surv_prob_5yr[1]+0.4, label="5 years", size=5, colour="grey60")+
  annotate("text", x=time_point_5yr-20, y = 0.85, label=paste0("HR for death: ",
                                                               round(PlotsICB$HR, digits = 2),"\n (95% CI: ", round(PlotsICB$lower, digits = 2),
                                                               "-",round(PlotsICB$upper, digits = 2), ") \n Log-rank \n", round(p_value, digits = 3), sep = ""), size=5, fontface="plain")
km

ggsave("/results/Fig2d.pdf", plot = km$plot, width = 8, height = 8)

#p = 0.039
#HR = 0.51
#CMV+ median = NA
#CMV- median = 34.4

################################################################################

#Perform a cox regression for OS according to CMV, correcting for age, sex, BRAF status 

sICB_cox <- coxph(Surv(months_death, death_status) ~ Age + Sex + `BRAF status`+ `CMV`,
                  data = subset(surv_last_update, combination == "sICB" & Rx != "Ipilimumab" & Rx != "RelatNivo"))
sICBsummary <- summary(sICB_cox)
coefsICB <- sICBsummary$coef[, 1]
cisICB <- as.data.frame(sICBsummary$conf.int)

#Create a data frame for plotting

PlotsICB <- data.frame(
  covariate = c("Age", "Female \n vs Male", "BRAF Wild-Type \n vs Mutant", "CMV+ \n vs CMV-"),
  HR = cisICB[, "exp(coef)"],
  lower = cisICB[,"lower .95"],
  upper = cisICB[,"upper .95"],
  p_value = sICBsummary$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)

PlotsICB <- PlotsICB %>% 
  mutate(Significant = ifelse(p_value < 0.05, "Significant", "Not significant")) %>% 
  mutate(Alpha = ifelse(Significant == "Significant", 1, 0.5))

################################################################################
#FIGURE 2e- Multivariate survival analysis in sICB-treated patients 

x <- ggplot(PlotsICB, aes(x = HR, y = covariate)) +
  geom_errorbar(data = subset(PlotsICB, Significant == "Significant"),
                aes(xmax = upper, xmin = lower), 
                width = 0, linetype = "longdash", colour = "#FFBB44", alpha = 1.2, size = 0.5) +
  geom_point(data = subset(PlotsICB, Significant == "Significant"), size = 5, alpha = 0.8, color = "#FFBB44") +
  geom_errorbar(data = subset(PlotsICB, Significant == "Not significant"),
                aes(xmax = upper, xmin = lower), 
                width = 0, linetype = "longdash", colour = "#2171B5", alpha = 1.2)+
  geom_point(data = subset(PlotsICB, Significant == "Not significant"), size = 5, alpha = 0.8, color = "#2171B5") +
  scale_x_log10()+
  geom_vline(xintercept = 1,linetype ="dotted", color = "grey") +
  xlab("Adjusted hazard ratios sICB")+
  ylab("Covariate")+
  theme_classic()+
  theme(axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(hjust=1, size=12),
        axis.text.y = element_text(size=12),
        axis.line = element_line(colour = "grey30", size = 0.2),
        axis.ticks = element_line(color= "grey30", size=0.2), 
        strip.background =element_rect(fill= "grey95"),
        strip.text = element_text(colour = "grey30"),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16))+
  geom_text(aes(label = paste0("p = ", signif(p_value, digits = 2))),
            vjust = -1.5, hjust = 0.5, size = 5, color = "grey30")

ggsave("/results/Fig2e.pdf", plot = x, width = 8, height = 8)

#CMV is significantly protective p = 0.0089, HR = 0.37 and age is a significant risk factor p = 0.0058, HR = 1.07

################################################################################

#Figure 2f- CMV impact on adjuvant recurrence-free survival 

Adjuvant <- clinicalData %>% 
  filter(!PtName %in% c("193", "197", "300", "408")) %>% # these are stage IV resected 
  filter(cancer == "Melanoma" & grepl("Adj", Rx)) %>% 
  filter(progression == 1 | months_progression > 6) %>% 
  mutate(months_death = as.numeric(months_death)) %>% 
  mutate(months_progression = as.numeric(months_progression)) %>% 
  mutate(months_death = ifelse(months_death > 90, 90, months_death)) %>% 
  mutate(months_progression = ifelse(months_progression > 90, 90, months_progression))
  
fit <- survfit(Surv(months_progression, progression) ~ cmv, data = Adjuvant)
summary(fit)$table
log_rank_test <- survdiff(Surv(months_progression, progression) ~ cmv, data = Adjuvant)
p_value <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)
cox <- coxph(Surv(months_progression, progression) ~ cmv, data = Adjuvant)
sICBsummary <- summary(cox)
coefsICB <- as.data.frame(sICBsummary$coef[, 1])
cisICB <- as.data.frame(sICBsummary$conf.int)
PlotsICB <- data.frame(
  covariate = c("CMV"),
  HR = cisICB[, "exp(coef)"],
  lower = cisICB[,"lower .95"],
  upper = cisICB[,"upper .95"],
  p_value = sICBsummary$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)

time_point_1yr <- 12
time_point_5yr <- 60
surv_prob_1yr <- summary(fit, times = time_point_1yr)$surv
surv_prob_5yr <- summary(fit, times = time_point_5yr)$surv
km <- ggsurvplot(fit,
                 # pval = TRUE,  # Add log rank p-value
                 # pval.method = TRUE,  # Add label for method used,
                 pval.size=5,
                 pval.coord = c(35, 0.70), pval.method.coord = c(35, 0.78),
                 xlab = "Months", ylab = "Probability of progression-free survival", legend.title = "CMV status", legend.labs = c("CMV-", "CMV+"),
                 font.main = 16,font.x =  16, font.y = 16, font.tickslab = 16, font.legend = 16, font.title = c(12, "bold"), break.x.by = 6, censor.shape = 124, panel.border = element_rect(colour = "black", fill=NA, size=4),
                 risk.table = "nrisk_cumevents",risk.table.legend = FALSE, risk.table.y.text.col= FALSE, risk.table.font = 3.5, risk.table.order = c(0,1),
                 ggtheme = theme_bw() + theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(),
                                              plot.title = element_text(hjust = 0.5)),
                 conf.int = F,
                 palette = c("#2171B5","#FFBB44"),
                 xlim = c(0, 72)
)

km$plot<- km$plot +
  geom_segment(x = time_point_1yr, xend = time_point_1yr, y = -Inf, yend = surv_prob_1yr[2] , color = "grey60", lty="dashed")+
  geom_segment(x = time_point_5yr, xend = time_point_5yr, y = -Inf, yend = surv_prob_1yr[2], color = "grey60", lty="dashed")+
  annotate("text", x=time_point_1yr+5, y=surv_prob_1yr[1] -0.18, label=paste0(round(surv_prob_1yr[[1]],3)*100,"%"), size=5, colour="#2171B5")+
  annotate("text", x=time_point_1yr+5, y=surv_prob_1yr[2]+0.05, label=paste0(round(surv_prob_1yr[[2]],3)*100,"%"), size=5, colour="#FFBB44")+
  annotate("text", x=time_point_1yr + 5, y=0, label="1 year", size=5, colour="grey60")+
  annotate("text", x=time_point_5yr+5, y=surv_prob_5yr[1]-0.05, label=paste0(round(surv_prob_5yr[[1]],3)*100,"%"), size=5, colour="#2171B5")+
  annotate("text", x=time_point_5yr+5, y=surv_prob_5yr[2]+0.05, label=paste0(round(surv_prob_5yr[[2]],3)*100,"%"), size=5, colour="#FFBB44")+
  annotate("text", x=time_point_5yr + 5, y=0, label="5 years", size=5, colour="grey60")+
  annotate("text", x=time_point_5yr-20, y = 0.70, label=paste0("HR for progression: ",
                                                               round(PlotsICB$HR, digits = 2),"\n (95% CI: ", round(PlotsICB$lower, digits = 2),
                                                               "-",round(PlotsICB$upper, digits = 2), ") \n Log-rank \n", round(p_value, digits = 3), sep = ""), size=5, fontface="plain")
km
ggsave("/results/Fig2f.pdf", plot = km$plot, width = 8, height = 8)

################################################################################
##FIGURE III- CMV INTERACTION WITH IRAES##

#Mutate ID to match survival data

TimeToTox <- tox %>% 
  mutate(patient = ID) %>% 
  mutate(grade = ifelse(is.na(grade), 0, grade)) %>% 
  dplyr::select(-ID) %>% 
  mutate(Grade_3 = ifelse(grade > 2, 1, 0)) %>% 
  mutate(Grade1_2 = ifelse(grade == 1 | grade == 2, 1, 0)) %>% 
  mutate(time_to_tox = as.numeric(time_to_tox))

#Get survival data for stage 2/3/4 Melanoma

SurvMel <- clinicalData %>% 
  mutate(patient = ID) %>% 
  filter(cancer == "Melanoma") %>% 
  dplyr::select(months_progression, months_death, progression, death_status, Rx, cmv, age, sex, BRAF_status, cancer, patient, subtype) %>% 
  mutate(months_death = as.numeric(months_death)) %>% 
  mutate(months_progression = as.numeric(months_progression)) %>% 
  unique()

#Join toxicity data with survival data 

SurvTox <- left_join(TimeToTox, SurvMel)

#Filter for Melanoma 

SurvToxMel <- SurvTox %>% 
  filter(cancer == "Melanoma")  %>% 
  mutate(Adjuvant = ifelse(Rx == "AdjNivo" | Rx == 'AdjPembro', "Adjuvant intent", "Palliative intent"))

#347 melanoma patients 
#308 metastatic melanoma and 39 stage II/III melanoma
#313 patients have CMV serology
#163 CMV- and 150 CMV+

#Modify time for grade 3+ toxicities  

Grade3 <- SurvToxMel %>% 
  dplyr::select(Grade_3, patient, time_to_tox, cmv, age, sex, Rx, BRAF_status, grade, months_death, death_status, Adjuvant) %>% 
  mutate(time_to_tox = ifelse(Grade_3 == 0, NA, time_to_tox))

#Identify first Grade 3 toxicity 

Grade3 <- Grade3 %>% 
  group_by(patient) %>% 
  slice_min(order_by = time_to_tox, with_ties = F) %>% 
  ungroup()

#Create competing risks column- censor at 24 months as there are no events beyond this point 

Grade3 <- Grade3 %>% 
  mutate(time_to_tox = ifelse(is.na(time_to_tox), months_death, time_to_tox)) %>% 
  mutate(time_to_tox = ifelse(time_to_tox > 24, 24, time_to_tox)) %>% 
  mutate(Status = as.factor(ifelse(Grade_3 == 1, 1, 
                                   ifelse(Grade_3 == 0 & death_status == 0, 0, 
                                          ifelse(death_status == 1 & Grade_3 == 0, 2, NA))))) %>% 
  filter(Rx != "RelatNivo" & Rx != "Ipilimumab") %>% 
  mutate(combination = ifelse(Rx == "IpiNivo", "cICB", "sICB"))

#342 patients excluding Ipi and RelatNivo
#308 with CMV serology 

################################################################################
#FIGURE 3a- Time to grade 3 toxicity according to CMV 

# Grade3 %>% 
#   filter(!is.na(cmv)) %>% 
#   dplyr::mutate(sex = ifelse(sex == 1, "Male", "Female")) %>% 
#   mutate(Treatment = combination) %>% 
#   mutate(across(where(is.character), as.factor)) %>% 
#   mutate(Intent = ifelse(Adjuvant == 0, "Adjuvant", "Palliative")) %>% 
#   dplyr::select(age, sex, BRAF_status, cmv, Treatment, Intent) %>% 
#   gtsummary::tbl_summary(by = cmv) %>% 
#   gtsummary::add_p() %>% 
#   modify_spanning_header(all_stat_cols() ~ "**irAE data**")

CumincCMV <- cuminc(Surv(time_to_tox, Status) ~ cmv , data = Grade3) 

pval <- t(as.data.frame(CumincCMV$cmprsk$Tests[1,2]))
formatted_pval <- paste0("Log-rank: \n ", signif(pval, digits = 2))

#Create the cumulative incidence plot

custom_theme <- theme_survminer()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(fill = F), 
        plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
x <- CumincCMV %>% 
  ggcuminc() +
  scale_x_continuous(breaks = seq(0, 24, 6))+
  labs(x = "Months", y = "Cumulative Incidence of Grade 3+ irAEs \n in Melanoma") +
  add_risktable(risktable_height = 0.3, size = 3.5, theme = custom_theme, 
                risktable_group = "risktable_stats",
                stats_label = c("Number at risk (number of events)"),
                risktable_stats = "{n.risk} ({cum.event})",
                risktable_width = 0.6,
                times = seq(0, 24, by = 6)) +
  add_censor_mark() + 
  theme_pubr(border = T)+
  theme(axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(hjust=1, size=14),
        axis.text.y = element_text(size=14),
        axis.line = element_line(colour = "grey30", size = 0.2),
        axis.ticks = element_line(color= "grey30", size=0.2), 
        strip.background =element_rect(fill= "grey95"),
        strip.text = element_text(colour = "grey30"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        strip.text.x = element_text(size = 18), 
        strip.text.y = element_text(size = 18),
        legend.position = "top",
        legend.text = element_text(size = 18))+
  scale_color_manual(values = c("#EB7926","#134B73"))+
  annotate("text", x = 6, y = 0.2, 
           label = formatted_pval, color = "black", size = 5)

ggsave("/results/Fig3a.pdf", plot = x, width = 8, height = 8)

#p = 2.2e-5
#CI at 6 months: 
#CMV+ = 0.30
#CMV- = 0.52

################################################################################

#Multivariate competing risks, accounting for age, sex, treatment, treatment intent, and BRAF status  

CMVgrade3Tox <- crr(Surv(time_to_tox, Status) ~ age + sex + combination + Adjuvant + BRAF_status + cmv, data = Grade3) %>% 
  tbl_regression(exp = TRUE)

#Extract coefficients 

coef <- CMVgrade3Tox$table_body
coef <- coef %>% 
  dplyr::select(estimate, conf.low, conf.high, p.value, term) %>% 
  na.omit()

################################################################################
#FIGURE 3b- Multivariate grade 3 toxicities HR

PlotCMV <- data.frame(
  covariate = c("Age", "Female vs \n Male","sICB vs \n cICB", "Palliative vs \n Adjuvant intent", "BRAF Wild-Type \n vs Mutant", "CMV+ vs \n CMV-"),
  HR = coef[, "estimate"],
  lower = coef[,"conf.low"],
  upper = coef[,"conf.high"],
  p_value = coef[, "p.value"],
  stringsAsFactors = FALSE
)

PlotCMV <- PlotCMV %>% 
  mutate(Significant = ifelse(p.value < 0.05, "Significant", "Not significant")) %>% 
  mutate(Alpha = ifelse(Significant == "Significant", 1, 0.5))

x <- ggplot(PlotCMV, aes(x = estimate, y = covariate)) +
  geom_errorbar(data = subset(PlotCMV, Significant == "Significant"),
                aes(xmax = conf.high, xmin = conf.low), 
                width = 0, linetype = "longdash", colour = "#EB7926", alpha = 1.2, size = 0.5) +
  geom_point(data = subset(PlotCMV, Significant == "Significant"), size = 5, alpha = 0.8, color = "#EB7926") +
  geom_errorbar(data = subset(PlotCMV, Significant == "Not significant"),
                aes(xmax = conf.high, xmin = conf.low), 
                width = 0, linetype = "longdash", colour = "#134B73", alpha = 1.2) +
  geom_point(data = subset(PlotCMV, Significant == "Not significant"), size = 5, alpha = 0.8, color = "#134B73") +
  scale_x_log10()+
  geom_vline(xintercept = 1,linetype ="dotted", color = "grey") +
  xlab("Adjusted hazard ratios \n Grade 3+ irAEs")+
  ylab("Covariate")+
  theme_classic()+
  theme(axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(hjust=1, size=12),
        axis.text.y = element_text(size=12),
        axis.line = element_line(colour = "grey30", size = 0.2),
        axis.ticks = element_line(color= "grey30", size=0.2), 
        strip.background =element_rect(fill= "grey95"),
        strip.text = element_text(colour = "grey30"),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16))+
  geom_text(aes(label = paste0("p = ", signif(p.value, digits = 2))),
            vjust = -1.5, hjust = 0.5, size = 5, color = "grey30")

ggsave("/results/Fig3b.pdf", plot = x, width = 8, height = 8)

#CMV is significantly protective against developing Grade 3+ irAEs p = 0.0058, HR= 0.60
#sICB treatment which encompasses both metastatic and adjuvant sICB treatment is also protective p = 6.9e-6, HR = 0.30

################################################################################

#Identify the OR of all grade, grade 1-2 toxicity, grade 3+ toxicity, steroids and second-line immunosuppressants 

SurvToxMel <- SurvToxMel %>% 
  filter(Rx != "RelatNivo" & Rx != "Ipilimumab") %>% 
  mutate(combination = ifelse(Rx == "IpiNivo", "cICB", "sICB"))

OverallTox <- SurvToxMel %>% 
  dplyr::select(patient, age, sex, BRAF_status, combination, cmv, Rx,  steroids, AI, grade, second_line, Adjuvant) %>%
  dplyr::select(-grade)

Grade1_2 <- SurvToxMel %>% 
  mutate(Grade1_2 = ifelse(grade == 2 | grade == 1, 1, 0)) %>% 
  dplyr::select(Grade1_2, patient, time_to_tox, cmv, age, sex, Rx, BRAF_status, grade, months_death, death_status, Adjuvant) %>% 
  mutate(time_to_tox = ifelse(Grade1_2 == 0, NA, time_to_tox))

#Identify first Grade 3 toxicity 

Grade1_2 <- Grade1_2 %>% 
  group_by(patient) %>% 
  slice_min(order_by = time_to_tox, with_ties = F) %>% 
  ungroup() %>% 
  dplyr::select(patient, Grade1_2)

Grade3 <- Grade3 %>% 
  dplyr::select(patient, Grade_3)

OverallTox <- left_join(OverallTox, Grade3)
OverallTox <- left_join(OverallTox, Grade1_2)

OverallTox <- OverallTox %>% 
  relocate(Adjuvant)

toxicities <- colnames(OverallTox)[9:13]

#Create an empty list to store results 

results <- list()

#Loop through each grade and management of toxicity 

for (tox in toxicities) {
  
  #Create a dataframe for each grade/management
  
  data <- OverallTox %>% 
    dplyr::select(cmv, age,sex, BRAF_status, Rx, all_of(tox), patient, combination, Adjuvant) %>%
    unique() %>% 
    mutate(!!tox := ifelse(is.na(!!sym(tox)), 0,
                           ifelse(!!sym(tox) == 0, 0, 1))) %>% 
    arrange(-!!sym(tox)) %>% 
    distinct(patient, .keep_all = TRUE)
  
  #Try to fit the logistic regression model
  
  try({
    
    #Fit logistic regression model with increased max iterations
    
    model <- glm(as.formula(paste(tox, "~ age + sex + combination + Adjuvant + BRAF_status + cmv")), 
                 data = data, family = binomial,
                 control = list(maxit = 100))
    
    #Check if model has converged
    
    if (!model$converged) {
      warning(print(paste("Model did not converge for column:", tox)))
      next
    }
    
    #Extract odds ratio and confidence intervals for CMV
    
    odds_ratio <- exp(coef(model)["cmvCMV+"])
    conf_int <- exp(confint(model)["cmvCMV+", ])
    
    #Extract p-value for CMV
    
    p_value <- summary(model)$coefficients["cmvCMV+", "Pr(>|z|)"]
    
    #Store results
    
    results[[tox]] <- c(Odds_Ratio = odds_ratio, 
                        Lower_CI = conf_int[1], 
                        Upper_CI = conf_int[2], 
                        P_Value = p_value)
  }, silent = TRUE)
}

#Combine results into a data frame

ORtoxCMV <- do.call(rbind, results) 

ORtoxCMV <- as.data.frame(ORtoxCMV)

#Filter out non-convergent results 

ORtoxCMV <- ORtoxCMV %>% 
  filter(!is.na(`Lower_CI.2.5 %`) & `Lower_CI.2.5 %` != 0)

ORtoxCMV$Toxicity <- rownames(ORtoxCMV)
results_df <- ORtoxCMV %>% 
  mutate(Significant = ifelse(P_Value < 0.05, "Significant", "Not significant")) %>% 
  mutate(Alpha = ifelse(Significant == "Significant", 1, 0.5)) %>% 
  mutate(Toxicity = ifelse(Toxicity == "steroids", "Steroids",
                           ifelse(Toxicity == "second_line", "Second line \n immunosuppressants",
                                  ifelse(Toxicity == "Grade1_2", "Grade 1-2 irAEs",
                                         ifelse(Toxicity == "Grade_3", "Grade 3+ irAEs", "All grade irAEs"))))) %>% 
  arrange(-`Odds_Ratio.cmvCMV+`) %>%
  mutate(Lower_CI_flag = ifelse(`Lower_CI.2.5 %` < 0.2, TRUE, FALSE)) %>% 
  mutate(`Lower_CI.2.5 %` = ifelse(Toxicity == "Second line \n immunosuppressants", 0.2, `Lower_CI.2.5 %`))

################################################################################
#FIGURE 3c- OR of toxicities according to CMV 

x <- ggplot(results_df, aes(x = `Odds_Ratio.cmvCMV+`, y = Toxicity)) +
  geom_errorbar(data = subset(results_df, Significant == "Significant"),
                aes(xmax = `Upper_CI.97.5 %`, xmin = `Lower_CI.2.5 %`), 
                width = 0, linetype = "longdash", colour = "#EF6548", alpha = 1.2, size = 0.5) +
  geom_point(data = subset(results_df, Significant == "Significant"), size = 5, alpha = 0.8, color = "#EF6548") +
  geom_errorbar(data = subset(results_df, Significant == "Not significant"),
                aes(xmax = `Upper_CI.97.5 %`, xmin = `Lower_CI.2.5 %`), 
                width = 0, linetype = "longdash", colour = "#9E9AC8", alpha = 1.2) +
  geom_segment(data = subset(results_df, Lower_CI_flag == TRUE),
               aes(xend = 0.2, y = Toxicity, yend = Toxicity),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "#EF6548", size = 0.5, linetype = "blank")+
  geom_point(data = subset(results_df, Significant == "Not significant"), size = 5, alpha = 0.8, color = "#9E9AC8") +
  scale_x_log10()+
  geom_hline(yintercept = 3.5, color = "grey30",linetype = "dashed")+
  geom_vline(xintercept = 1,linetype ="dotted", color = "grey") +
  xlab("Adjusted odds ratios in \n CMV seropositivity")+
  ylab("Management of irAEs")+
  theme_classic()+
  theme(axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(hjust=1, size=12),
        axis.text.y = element_text(size=12),
        axis.line = element_line(colour = "grey30", size = 0.2),
        axis.ticks = element_line(color= "grey30", size=0.2), 
        strip.background =element_rect(fill= "grey95"),
        strip.text = element_text(colour = "grey30"),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16))+
  geom_text(aes(label = paste0("p = ", signif(P_Value, digits = 2))),
            vjust = -1, hjust = 0.5, size = 5, color = "grey30")

ggsave("/results/Fig3c.pdf", plot = x, width = 8, height = 8)

#CMV seropositive patients have reduced Grade 3+ irAEs (p = 0.002, OR = 0.45), steroid requirement (p = 0.0032, OR = 0.46), and second line immunosuppression (p = 0.0076, OR = 0.40)

################################################################################

#dplyr::select specific toxicities and look at their incidence in CMV 

ToxTypes <- SurvToxMel %>% 
  pivot_wider(names_from = organ, values_from = AI) %>% 
  mutate(dermatitis = skin) %>% 
  mutate(`adrenal insufficiency` = adrenalitis) %>% 
  dplyr::select(-adrenalitis, -skin) %>% 
  mutate(Adjuvant = ifelse(grepl("Adj", Rx), "Adjuvant intent", "Palliative intent"))

#dplyr::select toxicities where there are a sufficient number of events 

Toxicities <- ToxTypes %>% 
  dplyr::select(hepatitis, arthritis, gastritis, thyroiditis, neuro, hypophysitis, nephritis, `adrenal insufficiency`, dermatitis, pneumonitis, colitis, myalgia)
toxicities <- colnames(Toxicities)

#Create an empty list to store results 

results <- list()

#Loop through each toxicity 

for (tox in toxicities) {
  
  #Prepare data for the current toxicity
  
  data <- ToxTypes %>% 
    dplyr::select(cmv, age,sex, BRAF_status, all_of(tox),Rx, patient, combination, Adjuvant) %>%
    unique() %>% 
    mutate(!!sym(tox) := ifelse(is.na(!!sym(tox)), 0,
                                ifelse(!!sym(tox) == 0, 0, 1))) %>% 
    arrange(-!!sym(tox)) %>% 
    distinct(patient, .keep_all = TRUE)
  
  #Try to fit the logistic regression model
  
  try({
    
    #Fit logistic regression model with increased max iterations
    
    model <- glm(as.formula(paste(tox, "~ age + sex + combination + BRAF_status + Adjuvant + cmv")), 
                 data = data, family = binomial,
                 control = list(maxit = 50))
    
    #Check if model has converged
    
    if (!model$converged) {
      warning(paste("Model did not converge for column:", tox))
      next
    }
    
    #Extract odds ratio and confidence intervals for CMV
    
    odds_ratio <- exp(coef(model)["cmvCMV+"])
    conf_int <- exp(confint(model)["cmvCMV+", ])
    
    #Extract p-value for CMV
    
    p_value <- summary(model)$coefficients["cmvCMV+", "Pr(>|z|)"]
    
    # Store results
    results[[tox]] <- c(Odds_Ratio = odds_ratio, 
                        Lower_CI = conf_int[1], 
                        Upper_CI = conf_int[2], 
                        P_Value = p_value)
  }, silent = TRUE)
}

#Combine results into a data frame

ORtoxCMV <- do.call(rbind, results) 

ORtoxCMV <- as.data.frame(ORtoxCMV)

#Remove toxicities where there was no convergence 

ORtoxCMV <- ORtoxCMV %>% 
  filter(!is.na(`Lower_CI.2.5 %`) & `Lower_CI.2.5 %` != 0)

ORtoxCMV$Toxicity <- rownames(ORtoxCMV)
results_df <- ORtoxCMV %>% 
  mutate(Significant = ifelse(P_Value < 0.05, "Significant", "Not significant")) %>% 
  mutate(Alpha = ifelse(Significant == "Significant", 1, 0.5)) %>% 
  arrange(-`Odds_Ratio.cmvCMV+`)
results_df <- results_df %>%
  mutate(Lower_CI_flag = ifelse(`Lower_CI.2.5 %` < 0.06, TRUE, FALSE)) %>% 
  mutate(`Lower_CI.2.5 %` = ifelse(Toxicity == "myalgia" | Toxicity == "pneumonitis", 0.06, `Lower_CI.2.5 %`))

################################################################################
#FIGURE 3d- OR of toxicities according to CMV 

#Plot 

x <- ggplot(results_df, aes(x = `Odds_Ratio.cmvCMV+`, y = Toxicity)) +
  geom_errorbar(data = subset(results_df, Significant == "Significant"),
                aes(xmax = `Upper_CI.97.5 %`, xmin = `Lower_CI.2.5 %`), 
                width = 0, linetype = "longdash", colour = "#EF6548", alpha = 1.2, size = 0.5) +
  geom_point(data = subset(results_df, Significant == "Significant"), size = 5, alpha = 0.8, color = "#EF6548") +
  geom_errorbar(data = subset(results_df, Significant == "Not significant"),
                aes(xmax = `Upper_CI.97.5 %`, xmin = `Lower_CI.2.5 %`), 
                width = 0, linetype = "longdash", colour = "#9E9AC8", alpha = 1.2) +
  geom_segment(data = subset(results_df, Lower_CI_flag == TRUE),
               aes(x = 0.06, xend = 0.06, y = Toxicity, yend = Toxicity),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "#EF6548", size = 0.5, linetype = "blank") +
  geom_point(data = subset(results_df, Significant == "Not significant"), size = 5, alpha = 0.8, color = "#9E9AC8") +
  scale_x_log10()+
  geom_vline(xintercept = 1,linetype ="dotted", color = "grey") +
  xlab("Adjusted odds ratios in \n CMV seropositivity")+
  ylab("All grade irAEs")+
  theme_classic()+
  theme(axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(hjust=1, size=12),
        axis.text.y = element_text(size=12),
        axis.line = element_line(colour = "grey30", size = 0.2),
        axis.ticks = element_line(color= "grey30", size=0.2), 
        strip.background =element_rect(fill= "grey95"),
        strip.text = element_text(colour = "grey30"),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16))+
  geom_text(aes(label = paste0("p = ", signif(P_Value, digits = 2))),
            vjust = -1, hjust = 0.5, size = 5, color = "grey30")

ggsave("/results/Fig3d.pdf", plot = x, width = 8, height = 8)

#Colitis: p = 0.00078, OR = 0.39
#Pneumonitis: p = 0.028, OR = 0.23
#Myalgia: p = 0.0091, OR = 0.15
#Dermatitis: p = 0.044, OR = 1.66

################################################################################

Colitis <- SurvToxMel %>% 
  filter(Rx == "IpiNivo") %>% 
  dplyr::select(Grade_3, patient, time_to_tox, cmv, age, sex, Rx, BRAF_status, grade, months_death, death_status, organ, cancer) %>%
  mutate(time_to_tox = ifelse(organ != "colitis", NA, time_to_tox))

#Identify first colitis toxicity 

Colitis <- Colitis %>% 
  group_by(patient) %>% 
  slice_min(order_by = time_to_tox, with_ties = F) %>% 
  ungroup() %>% 
  mutate(irAE = ifelse(organ == "colitis", "Colitis", "No/Other irAE")) 

################################################################################
#Figure 3e- OR of colitis according to CMV status in cICB-treated patients

fisherTest <- fisher.test(Colitis$irAE, Colitis$cmv)
OR <- fisherTest$estimate
P <- fisherTest$p.value

Summary <- Colitis %>%
  group_by(cmv, irAE) %>%
  summarise(count = n()) %>%
  ungroup()

p <- ggplot(subset(Summary, !is.na(cmv)), aes(x = irAE, y = count, fill = irAE)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, color = "grey45") +
  facet_wrap(~cmv) +
  labs(
    x = "irAE",
    y = "Number of cICB-treated \n patients") +
  scale_fill_brewer(palette = 3)+
  theme(legend.position = "None")

annotation_data <- data.frame(
  cmv = "CMV+",
  irAE = "Colitis",  
  count = 20,
  label = paste0("OR = ", round(OR, 1), "\nP = ", round(P, 4))
)

p <- p + geom_text(data = annotation_data, aes(label = label), 
                   x = 1, y = 30,
                   size = 5, hjust = 0.5)
p

ggsave("/results/Fig3e.pdf", plot = p, width = 8, height = 8)

################################################################################

#Internal validation of irAE association in Non-melanoma patients 

SurvNonMel <- clinicalData %>% 
  mutate(patient = ID) %>% 
  filter(cancer != "Melanoma") %>% 
  filter(!Rx %in% c("CAPOX", "FOLFOX", "FOLFIRI", "FOLFOX+panitumumab")) %>% 
  dplyr::select(months_progression, months_death, progression, death_status, Rx, cmv, age, sex, BRAF_status, cancer, patient, subtype) %>% 
  mutate(months_death = as.numeric(months_death)) %>% 
  mutate(months_progression = as.numeric(months_progression)) %>% 
  unique()

#join non-melanoma metadata with time to toxicity data 

SurvTox <- left_join(TimeToTox, SurvNonMel)

SurvToxNonMel <- SurvTox %>% 
  filter(!is.na(cancer))

#Pivot data 

ToxTypes <- SurvToxNonMel %>% 
  mutate(Adjuvant = ifelse(grepl("Adj", Rx), "Adjuvant intent", "Palliative intent")) %>% 
  dplyr::select(age, sex, cmv, cancer, patient, organ, Rx) %>% 
  unique() %>% 
  mutate(CMVprotected = ifelse(organ %in% c("myalgia", "colitis", "pneumonitis"), 1, 0)) %>% 
  mutate(OtherTox = ifelse(!organ %in% c("myalgia", "colitis", "pneumonitis", "none"), 1, 0)) %>% 
  mutate(skin = ifelse(organ == "skin", 1, 0)) %>% 
  dplyr::select(-organ) %>% 
  unique()

CMVprotected <- ToxTypes %>% 
  dplyr::select(-OtherTox, -skin) %>% 
  mutate(CMVprotected = ifelse(patient == 409 | patient == 298, 1, CMVprotected)) %>% 
  unique()

fisher.test(CMVprotected$cmv, CMVprotected$CMVprotected, alternative = "less")

# data:  CMVprotected$cmv and CMVprotected$CMVprotected
# p-value = 0.04395
# alternative hypothesis: true odds ratio is less than 1
# 95 percent confidence interval:
#   0.0000000 0.9432429
# sample estimates:
#   odds ratio 
# 0 

OtherTox <- ToxTypes %>% 
  dplyr::select(-CMVprotected, -skin) %>% 
  mutate(OtherTox = ifelse(patient == 409 | patient == 298, 1, OtherTox)) %>% 
  unique()

fisher.test(OtherTox$cmv, OtherTox$OtherTox, alternative = "less")

# data:  OtherTox$cmv and OtherTox$OtherTox
# p-value = 0.9162
# alternative hypothesis: true odds ratio is less than 1
# 95 percent confidence interval:
#   0.000000 4.987983
# sample estimates:
#   odds ratio 
# 1.794632 

################################################################################

##FIGURE IV- TRANSCRIPTION FACTOR DRIVERS OF THE CMV EFFECT##

#Transcription factor analysis

#Filter metastatic melanoma data for C1 expression and create a new sample variable 

metadata <- mmdata %>% 
  filter(patient != 25 & Rx != "Ipilimumab" & sample != "260_C2") %>% 
  filter(cycle == "C1") %>% 
  filter(!is.na(cmv)) %>% 
  mutate(sample = paste0(sample, "_CD8_None")) %>% 
  mutate(sample = as.factor(sample)) 

#For the batch_data, dplyr::select CD8+ T cell samples taken at baseline (format C1_CD8_None)

batch <- batch_data %>% 
  mutate(patient = as.numeric(patient)) %>% 
  mutate(sample = as.factor(sample)) %>% 
  filter(grepl(paste("C1_CD8_None", sep = "|"), sample))

#Merge batch data with metadata which can be used in DESeq2 for batch correction 
#merge using patient, cycle, and sample and discard anything that does not match (unmatched = "drop")

MetaBatch <- left_join(metadata, batch, by = c("patient", "cycle", "sample") , unmatched = "drop", keep = F)

NewSamples <- c(262,264,271, 274, 275, 276, 279, 281, 286, 291, 297, 302,265, 269, 273, 284, 288, 292)
MissingBatches2022.Sep <- c(233, 239, 243, 254, 261)
MissingBatchesMay2021 <- c(224)
MissingBatchesDec2020 <- c(161)
MissingBeforeSep2019 <- c(2, 26, 41, 44, 45, 81, 86)
MissingSep2019 <- c(153, 155)

FilteredMeta <- MetaBatch %>%
  mutate_if(is.character, as.factor) %>% 
  mutate(cmv = ifelse(cmv == "CMV+", "CMVpositive",
                      ifelse(cmv == "CMV-", "CMVnegative", NA))) %>%
  mutate(combination = ifelse(Rx == "IpiNivo", "cICB", "sICB")) %>% 
  mutate(origBatch = case_when(
    patient %in% NewSamples ~ "Jan2024",
    TRUE ~ origBatch)) %>% 
  mutate(origBatch = case_when(patient %in% MissingBatches2022.Sep ~ "2022.Sep", TRUE ~ origBatch)) %>% 
  mutate(origBatch = case_when(patient %in% MissingBatchesMay2021 ~ "May_2021", TRUE ~ origBatch)) %>% 
  mutate(origBatch = case_when(patient %in% MissingBatchesDec2020 ~ "Dec_2020", TRUE ~ origBatch)) %>% 
  mutate(origBatch = case_when(patient %in% MissingBeforeSep2019 ~ "Before_SEP2019", TRUE ~ origBatch)) %>% 
  mutate(origBatch = case_when(patient %in% MissingSep2019 ~ "Sep_19", TRUE ~ origBatch)) %>% 
  mutate(origBatch = ifelse(patient > 302, "June_2024", origBatch)) 

#Ensure the rownames still correspond to the sample ids

rownames(FilteredMeta) <- FilteredMeta$sample

#Subset counts data for baseline CD8+ T cell samples (C1_CD8_None)
#Add Ensembl ids as rownames in this new dataset

#rownames(counts) <- counts$Gene
SubsetCounts <- counts %>%
  dplyr::select(matches(c("C1_CD8_None")))
rownames(SubsetCounts) <- rownames(counts)

#Metadata and counts data do not have matching samples. Match them and subset the matches 
#Ensure that the counts column names are in the rownames of the metadata
#Find where the two dataframes intersect and filter according to the matches

all(colnames(SubsetCounts) %in% rownames(FilteredMeta))
matching <- intersect(colnames(SubsetCounts), rownames(FilteredMeta))
FilteredMeta <- subset(FilteredMeta, (sample %in% matching))
rownames(FilteredMeta) <- FilteredMeta$sample
SubsetCounts <- SubsetCounts %>% 
  dplyr::select(contains(matching))
rownames(SubsetCounts) <- rownames(counts)

#Not only do all the samples have to match in the metadata and the counts, but they have to be in the same order. Reorder the counts data in accordance with the row order in the metadata 
#dplyr::select metadata columns of interest: sample, sex, batch, age as continuous variable and cmv status- remove columns that are not required 

CountsReordered <- SubsetCounts[, match(rownames(FilteredMeta), colnames(SubsetCounts))]
all(colnames(CountsReordered) == rownames(FilteredMeta))
rownames(CountsReordered) = rownames(counts) 
CountsReordered$id <- rownames(CountsReordered)
FinalMetaC1 <- FilteredMeta %>% 
  dplyr::select(sample, patient, sex, age, cmv, origBatch, Rx) %>% 
  mutate(sex = as.factor(sex)) %>% 
  mutate(age = scale(age, center = T))
rownames(FinalMetaC1) <- FinalMetaC1$sample

#Check that columns and rows match again 

all(colnames(CountsReordered) %in% rownames(FinalMetaC1))
matching <- intersect(rownames(FinalMetaC1), colnames(CountsReordered))
CountsReordered2 <- CountsReordered %>% 
  dplyr::select(all_of(matching))
CountsReordered2$id <- CountsReordered$id

#Filter out genes which have low reads as these could easily be classed as differentially expressed. Vignette proposes using a count of 10 as a cutoff 
#Do a final check that the samples in the counts data are matching and in the same order as the samples in the metadata 

keep <- rowMeans(CountsReordered2[,1:(ncol(CountsReordered2)-1)]) >= 10
CountsReordered2 <- CountsReordered2[keep, ] 
CountsFinalC1 <- dplyr::select(CountsReordered2, -id)
rownames(CountsFinalC1) <- CountsReordered2$id
all(colnames(CountsFinalC1) == rownames(FinalMetaC1))

#CollecTRI is a comprehensie resource containing a curated collection of TFs and their transcriptional targets compiled from 12 different resources 
#Retrieve the human version from OmniPath 

Network <- get_collectri(organism = "human", split_complexes = F)

#Filter for transcription factors expressed in CD8+ T cells 

CountsFinalC1_2 <- CountsFinalC1
CountsFinalC1_2$symbol <- mapIds(org.Hs.eg.db, keys = rownames(CountsFinalC1_2), keytype = "ENSEMBL", column = "SYMBOL")
Network <- Network %>% 
  filter(source %in% CountsFinalC1_2$symbol)

metaTF <- FinalMetaC1 %>% 
  dplyr::select(sample, cmv)

ddsTF <- DESeqDataSetFromMatrix(countData = CountsFinalC1, 
                                colData = FinalMetaC1, 
                                design = ~ 1) #No model as we will perform vst transformation 

#Perform variance stabilisation to get normalised data

ddsNorm <- vst(ddsTF)

#Extract normalised counts and map ENSEMBL ids 
#Rownames must be gene names as opposed to ENSEMBL ids

NormCounts <- as.data.frame(assay(ddsNorm))
NormCounts$symbol <- mapIds(org.Hs.eg.db, keys = rownames(NormCounts), keytype = "ENSEMBL", column = "SYMBOL")
NormCounts <- NormCounts %>% 
  filter(!is.na(symbol)) 
NormCounts <- NormCounts %>%
  distinct(symbol, .keep_all = TRUE)
NormCountsFinalC1 <- NormCounts %>% 
  dplyr::select(-symbol)
rownames(NormCountsFinalC1) <- NormCounts$symbol

#Extract differentially expressed genes 

DEGs <- ResultsOrdered %>% 
  dplyr::select(symbol, log2FoldChange, padj, pvalue, stat) %>% 
  filter(!is.na(symbol)) %>% 
  distinct(symbol,.keep_all = T) %>% 
  filter(padj < 0.05)

DEGsFinal <- DEGs %>% 
  dplyr::select(-symbol)
rownames(DEGsFinal) <- DEGs$symbol

#Run ulm to find fold change of transcription factors according to CMV 

TF <- run_ulm(mat = DEGsFinal[, "stat", drop = F],
              network = Network,
              .source = 'source',
              .target = 'target',
              .mor = 'mor',
              minsize = 5)

#FDR for multiple testing 

pvals <- TF$p_value
padj <- p.adjust(pvals, method = "fdr")

TF$padj <- padj

#Filter toppadj#Filter top TFs in both signs

TF <- TF %>% 
  filter(padj < 0.05) %>% 
  mutate(Score = score)

#9 significantly associated transcription factors
#All 9 have induced activity in CMV seropositivity 

################################################################################

#Identify TF that are associated with CMV PC1

TFsInduced <- TF %>% 
  filter(score > 0) %>% 
  dplyr::select(source)

metadata <- mmdata %>% 
  filter(patient != 25 & Rx != "Ipilimumab" & sample != "260_C2") %>% 
  filter(cycle == "C1" | cycle == "C2") %>% 
  mutate(sample = paste0(sample, "_CD8_None")) %>% 
  mutate(sample = as.factor(sample)) 

#For the batch_data 

batch <- batch_data %>% 
  mutate(patient = as.numeric(patient)) %>% 
  mutate(sample = as.factor(sample)) %>% 
  filter(grepl(paste("C1_CD8_None|C2_CD8_None", collapse = "|"), sample))

#Merge batch data with metadata which can be used in DESeq2 for batch correction 
#merge using patient, cycle, and sample and discard anything that does not match

MetaBatch <- left_join(metadata, batch, by = c("patient", "cycle", "sample") , unmatched = "drop", keep = F)
FilteredMeta <- MetaBatch %>% 
  mutate(cmv = ifelse(cmv == "CMV+", "CMVpositive",
                      ifelse(cmv == "CMV-", "CMVnegative", NA)))#this is because DESEQ2 does not like non-conventional characters 

#Ensure that sample names are rownames 

rownames(FilteredMeta) <- FilteredMeta$sample

#Subset counts data for baseline and C2 CD8+ T cells

SubsetCounts <- counts %>% 
  dplyr::select(matches(c("C2_CD8_None", "C1_CD8_None")))
#rownames(SubsetCounts) <- counts$Gene

#Metadata and counts data do not have matching samples. Match them and subset the matches 
#see if the counts column names are in the rownames of the metadata
#Find where the two dataframes intersect and filter according to the matches

all(colnames(SubsetCounts) %in% rownames(FilteredMeta))
Matching <- intersect(colnames(SubsetCounts), rownames(FilteredMeta))
FilteredMeta <- subset(FilteredMeta, (sample %in% Matching))
rownames(FilteredMeta) <- FilteredMeta$sample
SubsetCounts <- SubsetCounts %>% 
  dplyr::select(contains(Matching))
#rownames(SubsetCounts) <- counts$Gene

#Not only do all the samples have to match in the metadata and the counts, but they have to be in the same order. Reorder the counts data 
#dplyr::select metadata columns of interest

SubsetReordered <- SubsetCounts[, match(rownames(FilteredMeta), colnames(SubsetCounts))]
all(colnames(SubsetReordered) == rownames(FilteredMeta))
rownames(SubsetReordered) = rownames(SubsetCounts) 
SubsetReordered$id <- rownames(SubsetReordered)
FinalMeta <- FilteredMeta %>% 
  dplyr::select(sample, sex, cycle, age, cmv, Rx, sample, origBatch, patient, BRAF_status, cancer)

#Filter out genes which have low reads as these could easily be classed as differentially expressed. 10 is the count threshold suggested by the vignette.
#Do a final check that the samples in the counts data are matching and in the same order as the samples in the metadata 

keep <- rowMeans(SubsetReordered[,1:(ncol(SubsetReordered)-1)]) >= 10
SubsetReordered <- SubsetReordered[keep, ] 
rownames(SubsetReordered) <- SubsetReordered$id
CountsFinal <- dplyr::select(SubsetReordered, -id)
rownames(CountsFinal) <- rownames(SubsetReordered)
rownames(FinalMeta) <- FinalMeta$sample
all(colnames(CountsFinal) == rownames(FinalMeta))

NormCounts <- baselineNormaliseTransform(files = CountsFinal, samples = FinalMeta)
NormCounts$symbol <- mapIds(org.Hs.eg.db, keys = rownames(NormCounts), keytype = "ENSEMBL", column = "SYMBOL")
TFinduced <- NormCounts %>% 
  dplyr::select(-X) %>% 
  filter(NormCounts$symbol %in% TFsInduced$source)
rownames(TFinduced) <- TFinduced$symbol
TFinduced <- TFinduced %>% 
  dplyr::select(-symbol)
TFinducedWide <- t(TFinduced)
rownames(TFinducedWide) <- FinalMeta$sample
TFinducedWide <- as.data.frame(TFinducedWide)
TFinducedWide$sample <- rownames(TFinducedWide)
CountsMeta <- left_join(FinalMeta, TFinducedWide)
CountsMeta <- left_join(CountsMeta, PCAmeta)
CountsMetaLong <- CountsMeta %>% 
  pivot_longer(names_to = "TF", values_to = "NormExp", c(11:19))

################################################################################
#FIGURE 4a- TBX21 correlation with CMV score 

CountsMetaLong <- CountsMetaLong %>% 
  mutate(combination = ifelse(Rx == "IpiNivo", "cICB", "sICB"))

x <- ggplot(subset(CountsMetaLong, cycle == "C2" & TF == "TBX21"), aes(x = PC1, y = NormExp))+
  xlab("CMV score")+
  ylab("Post cycle 1 TBX21 expression")+
  geom_point(color = "#005A32", alpha = 0.6, size = 3)+
  geom_smooth(method = "lm", linetype = "dashed", color = "#EB7926", alpha = 0.4, fill = "#C7E9C0")+
  stat_cor(method = "spearman", size = 6)+
  theme_classic()+
  theme(axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(hjust=1, size=14),
        axis.text.y = element_text(size=14),
        axis.line = element_line(colour = "grey30", size = 0.2),
        axis.ticks = element_line(color= "grey30", size=0.2), 
        strip.background =element_rect(fill= "grey95"),
        strip.text = element_text(colour = "grey30"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        strip.text.x = element_text(size = 18), 
        strip.text.y = element_text(size = 18))

ggsave("/results/Fig4a.pdf", plot = x, width = 8, height = 8)

#Spearman Rho of 0.75 p < 2.2e-16

################################################################################

#Explore T-bet dynamics in CMV and treatment type

metadata <- mmdata %>% 
  filter(patient != 25 & Rx != "Ipilimumab" & sample != "260_C2") %>% 
  filter(cycle == "C1"|cycle == "C2") %>% 
  mutate(sample = paste0(sample, "_CD8_None")) %>% 
  mutate(sample = as.factor(sample)) %>% 
  filter(!is.na(Rx)) %>% 
  mutate(combination = ifelse(Rx == "IpiNivo", "cICB", "sICB")) 

#For the batch_data, dplyr::select CD8 T cell samples taken at baseline (format C1_CD8_)

batch <- batch_data %>% 
  mutate(patient = as.numeric(patient)) %>% 
  mutate(sample = as.factor(sample)) %>% 
  filter(grepl(paste("C1_CD8_None|C2_CD8_None", sep = "|"), sample))

#Merge batch data with metadata which can be used in DESeq2 for batch correction 
#merge using patient, cycle, and sample and discard anything that does not match (unmatched = "drop")

MetaBatch <- left_join(metadata, batch, by = c("patient", "cycle", "sample") , unmatched = "drop", keep = F)

FilteredMeta <- MetaBatch %>% 
  mutate(over68 = ifelse(age >= 68, "68+", "<68")) %>% 
  mutate(cmv = ifelse(cmv == "CMV+", "CMVpositive",
                      ifelse(cmv == "CMV-", "CMVnegative", NA))) #this is because DESEQ2 does not like non-conventional characters 

#Ensure the rownames still correspond to the sample ids

rownames(FilteredMeta) <- FilteredMeta$sample

#Subset counts data for baseline CD8 T cell samples (C1_CD8_None format)
#Add Ensembl ids as rownames in this new dataset- still unsure why R always removes rownames when you modify a dataframe

#rownames(counts) <- counts$Gene
SubsetCounts <- counts %>%
  dplyr::select(matches(c("C1_CD8_None","C2_CD8_None")))
rownames(SubsetCounts) <- rownames(counts)

#Metadata and counts data do not have matching samples. Match them and subset the matches 
#Ensure that the counts column names are in the rownames of the metadata
#Find where the two dataframes intersect and filter according to the matches

all(colnames(SubsetCounts) %in% rownames(FilteredMeta))
Matching <- intersect(colnames(SubsetCounts), rownames(FilteredMeta))
FilteredMeta <- subset(FilteredMeta, (sample %in% Matching))
rownames(FilteredMeta) <- FilteredMeta$sample
SubsetCounts <- SubsetCounts %>% 
  dplyr::select(contains(Matching))
rownames(SubsetCounts) <- rownames(counts)

#Not only do all the samples have to match in the metadata and the counts, but they have to be in the same order. Reorder the counts data in accordance with the row order in the metadata 

SubsetReordered <- SubsetCounts[, match(rownames(FilteredMeta), colnames(SubsetCounts))]
all(colnames(SubsetReordered) == rownames(FilteredMeta))
rownames(SubsetReordered) <- rownames(SubsetCounts) 
SubsetReordered$id <- rownames(SubsetReordered)
FinalMeta <- FilteredMeta %>% 
  mutate(`Responder status` = ifelse(months_progression > 60, "Durable Responder",
                                     ifelse(progression == 1 & months_progression < 60, "Progressor", NA))) %>% 
  dplyr::select(sample, patient, sex, age, cmv, Rx, origBatch, combination, `Responder status`, cycle, progression, months_progression, BRAF_status) %>% 
  mutate(sex = as.factor(sex))
rownames(FinalMeta) <- FinalMeta$sample

all(colnames(SubsetReordered) %in% rownames(FinalMeta))
Matching <- intersect(rownames(FinalMeta), colnames(SubsetReordered))
SubsetReordered2 <- SubsetReordered %>% 
  dplyr::select(all_of(Matching))
SubsetReordered2$id <- SubsetReordered$id

#Filter out genes which have low reads as these could easily be classed as differentially expressed. Vignette proposes using a count of 10 as a cutoff 
#Do a final check that the samples in the counts data are matching and in the same order as the samples in the metadata 

keep <- rowMeans(SubsetReordered2[,1:(ncol(SubsetReordered2)-1)]) >= 10
SubsetReordered2 <- SubsetReordered2[keep, ] 
CountsFinal <- dplyr::select(SubsetReordered2, -id)
rownames(CountsFinal) <- SubsetReordered2$id
all(colnames(CountsFinal) == rownames(FinalMeta))

#Extract the normalised counts of all genes using the baselineNormaliseTransform function provided by Ben 

NormCounts <- baselineNormaliseTransform(files = CountsFinal, samples = FinalMeta)

#Ensure the rownames match the ENSEMBL IDs

NormCounts2 <- as.data.frame(NormCounts)

#Map the ENSEMBL ids to symbols using mapIds()

NormCounts2$symbol <- mapIds(org.Hs.eg.db, keys = rownames(NormCounts2), keytype = "ENSEMBL", column = "SYMBOL")

#Subset the counts to include gene of interest
#Merge normalised counts for genes of interest with metadata 

Subset <- NormCounts2 %>% 
  filter(symbol %in% TF$source) %>% 
  as.data.frame()
rownames(Subset) <- Subset$symbol
Subset <- Subset %>% 
  dplyr::select(-c(symbol, X))
SubsetWide <- t(Subset)
rownames(SubsetWide) <- FinalMeta$sample
SubsetWide <- as.data.frame(SubsetWide)
SubsetWide$sample <- rownames(SubsetWide)
FinalMeta <- left_join(FinalMeta, SubsetWide, by = "sample")

#Correlating T-bet with CD8+ T cell flow subsets 

TbetClones <- FinalMeta %>% 
  filter(grepl("_C1_", sample)) %>% 
  dplyr::select(TBX21, patient) %>% 
  mutate(patient = as.numeric(patient))

#Join flow data with metadata 

FlowMeta <- left_join(flow, mmdata)

#Identify major subsets of interest 

SubsetsCD8 <- c("TEM_CD8pos_CD3pos", "TNaive_CD8pos_CD3pos", "TCM_CD8pos_CD3pos", "TEMRA_CD8pos_CD3pos")
FlowMetaSubsetCD8 <- FlowMeta %>% 
  filter(cycle == "C1") %>% 
  filter(!is.na(cancer)) %>% 
  filter(metric %in% SubsetsCD8) %>% 
  mutate(metric = ifelse(metric == "TNaive_CD8pos_CD3pos", "T Naive",
                         ifelse(metric == "TCM_CD8pos_CD3pos", "TCM",
                                ifelse(metric == "TEMRA_CD8pos_CD3pos", "TEMRA", "TEM"))))

FlowMetaSubsetCD8$metric <- factor(FlowMetaSubsetCD8$metric, levels = c("T Naive", "TCM", "TEM", "TEMRA"))

#Join flow data with TBX21 counts 

TbetFlow <- left_join(FlowMetaSubsetCD8, TbetClones)

################################################################################
#FIGURE 4b- TBX21 correlation with CD8+ T cell subsets

# TbetFlow %>%  
#   filter(!is.na(TBX21)) %>% 
#   dplyr::mutate(sex = ifelse(sex == 1, "Male", "Female")) %>% 
#   mutate(across(where(is.character), as.factor)) %>% 
#   dplyr::select(age, sex, BRAF_status, cmv, patient) %>% 
#   unique() %>% 
#   dplyr::select(-patient) %>% 
#   gtsummary::tbl_summary(by = cmv) %>% 
#   gtsummary::add_p() %>% 
#   modify_spanning_header(all_stat_cols() ~ "**Baseline TBX21-Flow correlation**")

x <- ggplot(TbetFlow, aes(x = TBX21, y = value))+
  xlab("Normalised TBX21 expression")+
  ylab("Baseline proportion of CD8+ T cells (%)")+
  geom_point(color = "#005A32", alpha = 0.6, size = 2)+
  geom_smooth(method = "lm", linetype = "dashed", color = "#EB7926", alpha = 0.4, fill = "#C7E9C0")+
  theme_classic()+
  theme(axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(hjust=1, size=14),
        axis.text.y = element_text(size=14),
        axis.line = element_line(colour = "grey30", size = 0.2),
        axis.ticks = element_line(color= "grey30", size=0.2), 
        strip.background =element_rect(fill= "grey95"),
        strip.text = element_text(colour = "grey30"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        strip.text.x = element_text(size = 18), 
        strip.text.y = element_text(size = 18))+
  facet_wrap(~metric, scales = "free", nrow = 1)+
  stat_cor(data = TbetFlow, method = "spearman", size = 6)

ggsave("/results/Fig4b.pdf", plot = x, width = 8, height = 8)

#89 baseline samples with Flow and bulk data 

################################################################################
#FIGURE 4c- TBX21 dyanamics with treatment according to CMV 

PlotData <- FinalMeta %>% 
  group_by(patient) %>% 
  filter(n_distinct(cycle) > 1) %>% 
  ungroup() %>% 
  filter(!is.na(cmv)) %>% 
  mutate(cmv = ifelse(cmv == "CMVpositive", "CMV+", "CMV-")) %>% 
  mutate(cycle = ifelse(cycle == "C1", "Baseline", "Post cycle 1"))

x <- ggplot(PlotData, aes(x = cycle, y = TBX21))+
  ggbeeswarm::geom_quasirandom(aes(fill = cycle),
                               colour='grey30',
                               size=1.5,
                               pch=21,
                               alpha=0.6,
                               width=0.075)+   
  geom_line(aes(group = patient), color = "grey", alpha = 0.5, linetype = "dashed") +
  geom_boxplot(colour = "grey30", aes(fill = cycle), alpha=0.5, width = 0.5)+
  scale_fill_brewer(palette = 11)+
  facet_grid(vars(combination), vars(cmv))+
  xlab("Cycle")+
  ylab("Normalised TBX21 expression")+
  theme(legend.position = "none")+
  labs(fill = "Cycle", color = "Cycle")+
  # stat_compare_means(method = "wilcox", label.x.npc = "center", paired = T, comparisons = list(c("Baseline", "Post cycle 1")), size = 5)
  stat_compare_means(comparisons = list(c("Baseline", "Post cycle 1")), paired = T, label = "p.signif", size = 6)

ggsave("/results/Fig4c.pdf", plot = x, width = 8, height = 8)

  #CMV- cICB: p = 0.0033
  #CMV- sICB: p = 0.0026
  #CMV+ cICB: p = 1 
  #CMV+ sICB: p = 0.96

################################################################################

SubsetGZMB <- FlowMeta %>% 
  filter(cycle == "C2"|cycle == "C1") %>% 
  mutate(cycle = ifelse(cycle == "C1", "Baseline", "Post cycle 1")) %>% 
  filter(!is.na(cancer)) %>%
  filter(!is.na(cmv)) %>% 
  mutate(combination = ifelse(Rx == "IpiNivo", "cICB", "sICB")) %>% 
  filter(metric == "GZMBpos_CD8pos_CD3pos") %>% 
  group_by(patient) %>% 
  filter(n_distinct(cycle) >1) %>% 
  ungroup()

################################################################################
#FIGURE 4d- Convergence of GZMB+ CD8+ T cells in cICB-treated, which is not observed in sICB-treated 

# SubsetGZMB %>%  
#   filter(cycle == "Baseline") %>% 
#   mutate(Treatment = combination) %>% 
#   dplyr::mutate(sex = ifelse(sex == 1, "Male", "Female")) %>% 
#   mutate(across(where(is.character), as.factor)) %>% 
#   dplyr::select(age, sex, BRAF_status, cmv, Treatment) %>% 
#   gtsummary::tbl_summary(by = cmv) %>% 
#   gtsummary::add_p() %>% 
#   modify_spanning_header(all_stat_cols() ~ "**Paired GZMB+ Flow data**")

x <- ggplot(SubsetGZMB, aes(x = cycle, y = value))+
  ggbeeswarm::geom_quasirandom(aes(fill= cycle),
                               colour='grey30',
                               size=1.5,
                               pch=21,
                               alpha=0.6,
                               width=0.075)+   
  geom_boxplot(colour = "grey30", aes(fill = cycle), alpha=0.5, width = 0.5)+
  geom_line(aes(group = patient), color = "grey", alpha = 0.5, linetype = "dashed") +
  scale_fill_brewer(palette = 2)+
  theme(legend.position = "none")+
  facet_grid(combination~cmv)+
  xlab("Cycle")+
  ylab("Proportion of GZMB+ CD8+ T cells (%)")+
  #stat_compare_means(method = "wilcox", label.x.npc = "center", comparisons = list(c("Baseline", "Post cycle 1")), paired =T, size = 5)
  stat_compare_means(comparisons = list(c("Baseline", "Post cycle 1")), paired = T, label = "p.signif", size = 6)

ggsave("/results/Fig4d.pdf", plot = x, width = 8, height = 8)

  #CMV- cICB: p = 0.0024
  #CMV- sICB: p = 0.84 
  #CMV+ cICB: p = 0.12
  #CMV+ sICB: p = 0.44
  
################################################################################

#Determine which patients had a durable response to treatment 

#samples from 23 patients with metastatic melanoma

cd8@meta.data <- cd8@meta.data %>% 
  mutate(Response = ifelse(progression == 1 & months_progression < 36, "Progressive disease",
                           ifelse(months_progression >= 36, "Durable Response", NA)))

cd8@meta.data <- cd8@meta.data %>% 
  mutate(Subset = ifelse(celltype.l2 == "CD8.TEMRA" | celltype.l2 == "CD8.TEMRA.CMC1", "Effector",
                         ifelse(celltype.l2 == "CD8.TEM", "EM",
                                ifelse(celltype.l2 == "CD8.TCM.CCL5" | celltype.l2 == "CD8.TCM", "CM",
                                       ifelse(celltype.l2 == "MAIT", "MAIT",
                                              ifelse(celltype.l2 == "CD8.Proliferating", "Proliferating",
                                                     ifelse(celltype.l2 == "CD8.Naive", "Naive",
                                                            ifelse(celltype.l2 == "dnT", "DN", 
                                                                   ifelse(celltype.l2 == "gdT", "GD", NA)))))))))

################################################################################
#FIGURE 4e- TBX21 by CD8+ T cell subset 

medianTBX21 <- cd8@meta.data %>%
  filter(cycle == "C1") %>% 
  filter(combination == "sICB") %>% 
  group_by(Subset) %>%
  summarise(median_expression = median(TBX21imputed, na.rm = TRUE))

cd8@meta.data <- cd8@meta.data %>%
  mutate(Subset = factor(Subset, levels = medianTBX21$Subset[order(medianTBX21$median_expression)]))

x <- ggplot(subset(cd8@meta.data, combination == "sICB" & !is.na(Subset) & cycle == "C1"), aes(x = TBX21imputed, y = Subset,  fill = Subset), alpha = 0.5)+
  geom_density_ridges()+
  ylab("")+
  xlab("Baseline imputed TBX21 expression")+
  theme_classic()+
  scale_fill_brewer(palette = 5)+
  theme(axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(hjust=1, size=14),
        axis.text.y = element_text(size=14),
        axis.line = element_line(colour = "grey30", size = 0.2),
        axis.ticks = element_line(color= "grey30", size=0.2), 
        strip.background =element_rect(fill= "grey95"),
        strip.text = element_text(colour = "grey30"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        strip.text.x = element_text(size = 18), 
        strip.text.y = element_text(size = 18),
        legend.text = element_text(size =14),
        legend.title = element_text(size = 18))

ggsave("/results/Fig4e.pdf", plot = x, width = 8, height = 8)

################################################################################

#Pseudobulked TBX21 by clone in patients treated with sICB

#Determine the clone count of each effector clone 

CloneSizes <- cd8@meta.data %>%
  filter(!is.na(TBX21imputed)) %>% 
  mutate(combination = ifelse(Rx == "IpiNivo", 'cICB', "sICB")) %>% 
  filter(Rx == "Nivo" | Rx == "Pembro") %>% 
  mutate(sample = paste(patient, cycle, sep = "_")) %>% 
  filter(celltype.l2 %in% c('CD8.TEMRA', 'CD8.TEMRA.CMC1')) %>% 
  filter(cycle == "C1" | cycle == "C2") %>% 
  mutate(cycle = ifelse(cycle == "C1", "Baseline",
                        ifelse(cycle == "C2", "Post cycle 1", NA))) %>% 
  filter(!is.na(clone_id) & clone_id != "NA" & nTRB < 2) %>% 
  mutate(Clone = paste(TCR_TRA1_v_gene,TCR_TRA1_j_gene, TCR_TRA1_cdr3, TCR_TRB1_chain, TCR_TRB1_v_gene,TCR_TRB1_j_gene, TCR_TRB1_cdr3, sep = "_")) %>% 
  group_by(sample, Clone) %>% 
  summarise(cloneCount = n()) %>% 
  ungroup() 

#Determine the total count per sample 

TotalCount <- cd8@meta.data %>%
  filter(!is.na(TBX21imputed)) %>% 
  mutate(combination = ifelse(Rx == "IpiNivo", 'cICB', "sICB")) %>% 
  filter(Rx == "Nivo" | Rx == "Pembro") %>% 
  mutate(sample = paste(patient, cycle, sep = "_")) %>% 
  filter(cycle == "C1" | cycle == "C2") %>% 
  mutate(cycle = ifelse(cycle == "C1", "Baseline",
                        ifelse(cycle == "C2", "Post cycle 1", NA))) %>% 
  filter(!is.na(clone_id) & clone_id != "NA" & nTRB < 2) %>% 
  group_by(sample) %>% 
  summarise(TotalCount = n())

CloneSizes <- left_join(CloneSizes, TotalCount)

#Determine the repertoire occupancy of each clone 

CloneSizes <- CloneSizes %>% 
  mutate(Proportion = cloneCount/TotalCount * 100) %>% 
  mutate(Size = ifelse(Proportion > 0.5, "> 0.5%",
                       ifelse(Proportion <= 0.5, "< 0.5%", NA)))

#Determine the median imputed TBX21 expression of each effector clone 

TBX21clones <- cd8@meta.data %>%
  filter(!is.na(TBX21imputed)) %>% 
  mutate(combination = ifelse(Rx == "IpiNivo", 'cICB', "sICB")) %>% 
  filter(Rx == "Nivo" | Rx == "Pembro") %>% 
  mutate(sample = paste(patient, cycle, sep = "_")) %>% 
  filter(celltype.l2 %in% c('CD8.TEMRA', 'CD8.TEMRA.CMC1')) %>% 
  filter(cycle == "C1" | cycle == "C2") %>% 
  mutate(cycle = ifelse(cycle == "C1", "Baseline",
                        ifelse(cycle == "C2", "Post cycle 1", NA))) %>% 
  filter(!is.na(clone_id) & clone_id != "NA" & nTRB < 2) %>% 
  mutate(Clone = paste(TCR_TRA1_v_gene,TCR_TRA1_j_gene, TCR_TRA1_cdr3, TCR_TRB1_chain, TCR_TRB1_v_gene,TCR_TRB1_j_gene, TCR_TRB1_cdr3, sep = "_")) %>% 
  group_by(sample, Clone) %>% 
  mutate(tbet = median(TBX21imputed)) %>% 
  distinct(sample, Clone, .keep_all=T) %>%
  ungroup()

#Median 1572 cells at baseline per sample 
#Median 1683 cells at post cycle 1 per sample 
#2988 effector clones from 18 sICB-treated patients 
#1431 clones from baseline samples and #1557 from post cycle 1 samples 

################################################################################
#FIGURE 4f- Imputed TBX21 expression according to effector clone size at day 21

TBX21clones <- left_join(TBX21clones, CloneSizes)
TBX21clones$Size <- factor(TBX21clones$Size, levels = c("< 0.5%", "> 0.5%"))

x <- ggplot(subset(TBX21clones, cycle == "Post cycle 1"), aes(x = Size, y = tbet))+
  ggbeeswarm::geom_quasirandom(aes(fill = Size),
                               colour='grey30',
                               size=1.5,
                               pch=21,
                               alpha=0.6,
                               width=0.075)+
  geom_boxplot(colour = "grey30", aes(fill = Size), alpha=0.5, width = 0.5)+
  scale_fill_brewer(palette = 3)+
  theme(legend.position = "none")+
  xlab("Clonal repertoire occupancy")+
  ylab("Median imputed TBX21 expression \n by effector clone")+
  stat_compare_means(comparisons = list(c("> 0.5%", "< 0.5%")), label = "p.signif", size = 7)
  #stat_compare_means(method = 'wilcox', comparisons = list(c("> 0.5%", "< 0.5%")), size = 5)

ggsave("/results/Fig4f.pdf", plot = x, width = 8, height = 8)

#p = 8.5e-7

################################################################################
#Figure 4g- CMV-reactive TRB and TBX21

CMVreactive <- cd8@meta.data %>% 
  filter(!is.na(TBX21imputed)) %>% 
  filter(!is.na(clone_id) & clone_id != "NA" & nTRB < 2) %>% 
  mutate(combination = ifelse(Rx == "IpiNivo", 'cICB', "sICB")) %>% 
  filter(Rx == "Nivo" | Rx == "Pembro") %>% 
  mutate(sample = paste(patient, cycle, sep = "_")) %>% 
  filter(celltype.l2 %in% c('CD8.TEMRA', 'CD8.TEMRA.CMC1')) %>% 
  filter(cycle == "C1" | cycle == "C2")

x <- ggplot(subset(CMVreactive, cycle == "C1"), aes(x = cmv_trb, y = TBX21imputed))+
  ggbeeswarm::geom_quasirandom(aes(fill = cmv_trb),
                               colour='grey30',
                               size=1.5,
                               pch=21,
                               alpha=0.6,
                               width=0.075)+
  geom_boxplot(colour = "grey30", aes(fill = cmv_trb), alpha=0.5, width = 0.5)+
  scale_fill_brewer(palette = 4)+
  theme(legend.position = "none")+
  xlab("CMV-reactive TRB")+
  ylab("Imputed TBX21 expression \n by effector cell")+
  stat_compare_means(comparisons = list(c("Putative CMV clone", "Not CMV-reactive")), label = "p.signif", size = 7)
#stat_compare_means(method = 'wilcox', comparisons = list(c("Putative CMV clone", "Not CMV-reactive")), size = 5)

ggsave("/results/Fig4g.pdf", plot = x, width = 8, height = 8)

#CMV beta chains have higher imputed TBX21 in the effector subset- 1.6e-14

################################################################################

#Figure 4i- Explore the bystander effect of CMV on effector clones 

CMVreactive2 <- cd8@meta.data %>%
  filter(!is.na(TBX21imputed)) %>% 
  mutate(combination = ifelse(Rx == "IpiNivo", 'cICB', "sICB")) %>% 
  filter(Rx == "Nivo" | Rx == "Pembro") %>% 
  mutate(sample = paste(patient, cycle, sep = "_")) %>% 
  filter(celltype.l2 %in% c('CD8.TEMRA', 'CD8.TEMRA.CMC1')) %>% 
  filter(cycle == "C1" | cycle == "C2") %>% 
  mutate(cycle = ifelse(cycle == "C1", "Baseline",
                        ifelse(cycle == "C2", "Post cycle 1", NA))) %>% 
  filter(!is.na(clone_id) & clone_id != "NA" & nTRB < 2) %>% 
  mutate(Clone = paste(TCR_TRA1_v_gene,TCR_TRA1_j_gene, TCR_TRA1_cdr3, TCR_TRB1_chain, TCR_TRB1_v_gene,TCR_TRB1_j_gene, TCR_TRB1_cdr3, sep = "_")) %>% 
  group_by(sample, Clone) %>% 
  mutate(tbet = median(TBX21imputed)) %>% 
  distinct(sample, Clone, .keep_all=T) %>%
  ungroup() %>% 
  group_by(patient, Clone) %>% 
  filter(n_distinct(cycle) >1 ) %>% 
  ungroup() %>% 
  mutate(id = paste(patient, Clone, sep = "_"))

x <- ggplot(subset(CMVreactive2, cmv_trb == "Not CMV-reactive"), aes(x = cycle, y = tbet))+
  geom_line(aes(group = id), color = "grey", alpha = 0.5, linetype = "dashed") +
  ggbeeswarm::geom_quasirandom(aes(fill = cycle),
                               colour='grey30',
                               size=1.5,
                               pch=21,
                               alpha=0.6,
                               width=0.075)+   
  geom_boxplot(colour = "grey30", aes(fill = cycle), alpha=0.5, width = 0.5)+
  scale_fill_brewer(palette = 4)+
  facet_wrap(vars(CMV))+
  theme(legend.position = "none")+
  theme(axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(hjust=1, size=14),
        axis.text.y = element_text(size=14),
        axis.line = element_line(colour = "grey30", size = 0.2),
        axis.ticks = element_line(color= "grey30", size=0.2), 
        strip.background =element_rect(fill= "grey95"),
        strip.text = element_text(colour = "grey30"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        strip.text.x = element_text(size = 18), 
        strip.text.y = element_text(size = 18),
        legend.text = element_text(size =18),
        legend.title = element_text(size = 18))+
  xlab("Cycle")+
  ylab("Median imputed TBX21 expression \n by effector clone")+
  labs(fill = "Response", color = "Response")+
  #stat_compare_means(method = "wilcox", label.x.npc = "center", comparisons = list(c("Baseline", "Post cycle 1")), size = 5, paired =T)
  stat_compare_means(comparisons = list(c("Baseline", "Post cycle 1")), label = "p.signif", size = 6, paired = T)

ggsave("/results/Fig4i.pdf", plot = x, width = 8, height = 8)

################################################################################

##FIGURE V- TBX21 DETERMINES SURVIVAL##

################################################################################
#FIGURE 5a- Imputed TBX21 expression in durable response and progressive disease- disease control at 3 years

TBX21response <- TBX21clones %>% 
  group_by(patient, Clone) %>% 
  filter(n_distinct(cycle) >1 ) %>% 
  ungroup() %>% 
  mutate(id = paste(patient, Clone, sep = "_"))

x <- ggplot(subset(TBX21response, !is.na(Response)), aes(x = cycle, y = tbet))+
  geom_line(aes(group = id), color = "grey", alpha = 0.5, linetype = "dashed") +
  ggbeeswarm::geom_quasirandom(aes(fill = cycle),
                               colour='grey30',
                               size=1.5,
                               pch=21,
                               alpha=0.6,
                               width=0.075)+   
  geom_boxplot(colour = "grey30", aes(fill = cycle), alpha=0.5, width = 0.5)+
  scale_fill_brewer(palette = 13)+
  facet_wrap(vars(Response))+
  theme(legend.position = "none")+
  theme(axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(hjust=1, size=14),
        axis.text.y = element_text(size=14),
        axis.line = element_line(colour = "grey30", size = 0.2),
        axis.ticks = element_line(color= "grey30", size=0.2), 
        strip.background =element_rect(fill= "grey95"),
        strip.text = element_text(colour = "grey30"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        strip.text.x = element_text(size = 14), 
        strip.text.y = element_text(size = 18),
        legend.text = element_text(size =14),
        legend.title = element_text(size = 18))+
  xlab("Cycle")+
  ylab("Median imputed TBX21 expression \n by effector clone")+
  labs(fill = "Response", color = "Response")+
  #stat_compare_means(method = "wilcox", label.x.npc = "center", comparisons = list(c("Baseline", "Post cycle 1")), size = 5, paired =T)
  stat_compare_means(comparisons = list(c("Baseline", "Post cycle 1")), label = "p.signif", size = 7, paired = T)

ggsave("/results/Fig5a.pdf", plot = x, width = 8, height = 8)

  #Durable response: p = 6.2e-11
  #Progressive disease: p < 2.22e-16

################################################################################
#FIGURE 5b- TBX21 expression according to response at 3y

PlotData <- FinalMeta %>% 
  mutate(progression = ifelse(months_progression > 36, "Durable Response",
                              ifelse(months_progression < 36 & progression == 1, "Progressive disease", NA))) %>% 
  mutate(cycle = ifelse(cycle == "C1", "Baseline", 'Post cycle 1')) 

# PlotData %>%  
#   filter(cycle == "Post cycle 1") %>% 
#   mutate(Treatment = combination) %>% 
#   mutate(CMV = ifelse(cmv == "CMVpositive", "CMV+", "CMV-")) %>% 
#   dplyr::mutate(sex = ifelse(sex == 1, "Male", "Female")) %>% 
#   mutate(across(where(is.character), as.factor)) %>% 
#   dplyr::select(age, sex, BRAF_status, CMV, Treatment) %>% 
#   gtsummary::tbl_summary(by = CMV) %>% 
#   gtsummary::add_p() %>% 
#   modify_spanning_header(all_stat_cols() ~ "**Post cycle 1 CD8+ Bulk RNA-seq data**")

x <- ggplot(subset(PlotData, !is.na(progression)), aes(x = progression, y = TBX21))+
  ggbeeswarm::geom_quasirandom(aes(fill = progression),
                               colour='grey30',
                               size=1.5,
                               pch=21,
                               alpha=0.6,
                               width=0.075)+   
  geom_boxplot(colour = "grey30", aes(fill = progression), alpha=0.5, width = 0.5)+
  scale_fill_brewer(palette = 13)+
  facet_wrap(vars(cycle))+
  theme(legend.position = "none")+
  xlab("Disease control at 3 years")+
  ylab("Normalised TBX21 expression")+
  labs(fill = "Disease control at 3 years", color = "Disease control at 3 years")+
  #stat_compare_means(method = "wilcox", label.x.npc = "center", comparisons = list(c("Durable Response", "Progressive disease")), size = 5)
  stat_compare_means(comparisons = list(c("Durable Response", "Progressive disease")), label = "p.signif", size = 7)

ggsave("/results/Fig5b.pdf", plot = x, width = 8, height = 8)

#Baseline: p =. 0.049 
#Post cycle 1: p = 0.008

################################################################################

#Check to see if module expression is associated with TF count at C2

FinalMetaC2 <- FinalMeta %>% 
  filter(cycle == "C2") %>% 
  mutate(patient = as.character(patient)) %>% 
  mutate(ID = as.numeric(patient))

#Survival according to TBX21 expression at C2 

survivalTBX21 <- FinalMetaC2 %>% 
  dplyr::select(ID, TBX21)

survivalTBX21 <- left_join(surv_last_update, survivalTBX21)

################################################################################
#FIGURE 5c- PFS Kaplan Meier of TBX21 expression post-cycle 1

MedianTBX21median <- median(subset(survivalTBX21, !is.na(TBX21))$TBX21)
survivalTBX21 <- survivalTBX21 %>% 
  mutate(DichotTBX21 = ifelse(TBX21 >= MedianTBX21median, "> Median TBX21", "< Median TBX21"))

fit <- survfit(Surv(months_progression, progression) ~ DichotTBX21, data = survivalTBX21)
summary(fit)$table
log_rank_test <- survdiff(Surv(months_progression, progression) ~ DichotTBX21, data = survivalTBX21)
p_value <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)
cox <- coxph(Surv(months_progression, progression) ~ DichotTBX21, data = survivalTBX21)
sICBsummary <- summary(cox)
coefsICB <- as.data.frame(sICBsummary$coef[, 1])
cisICB <- as.data.frame(sICBsummary$conf.int)
PlotsICB <- data.frame(
  covariate = c("TBX21 expression"),
  HR = cisICB[, "exp(coef)"],
  lower = cisICB[,"lower .95"],
  upper = cisICB[,"upper .95"],
  p_value = sICBsummary$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)

time_point_1yr <- 12
time_point_5yr <- 60
surv_prob_1yr <- summary(fit, times = time_point_1yr)$surv
surv_prob_5yr <- summary(fit, times = time_point_5yr)$surv
km <- ggsurvplot(fit,
                 # pval = TRUE,  # Add log rank p-value
                 # pval.method = TRUE,  # Add label for method used,
                 pval.size=5,
                 pval.coord = c(35, 0.70), pval.method.coord = c(35, 0.78),
                 xlab = "Months", ylab = "Probability of progression-free survival", legend.title = "TBX21 expression", legend.labs = c("< Median", "> Median"),
                 font.main = 16,font.x =  16, font.y = 16, font.tickslab = 16, font.legend = 16, font.title = c(12, "bold"), break.x.by = 6, censor.shape = 124, panel.border = element_rect(colour = "black", fill=NA, size=4),
                 risk.table = "nrisk_cumevents",risk.table.legend = FALSE, risk.table.y.text.col= FALSE, risk.table.font = 3.5, risk.table.order = c(0,1),
                 ggtheme = theme_bw() + theme(panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(),
                                              plot.title = element_text(hjust = 0.5)),
                 conf.int = F,
                 palette = c("#2171B5","#FFBB44"),
                 xlim = c(0, 90)
)

km$plot<- km$plot +
  geom_segment(x = time_point_1yr, xend = time_point_1yr, y = -Inf, yend = 0.75 , color = "grey60", lty="dashed")+
  geom_segment(x = time_point_5yr, xend = time_point_5yr, y = -Inf, yend = 0.60, color = "grey60", lty="dashed")+
  annotate("text", x=time_point_1yr+4, y=surv_prob_1yr[1] -0.1, label=paste0(round(surv_prob_1yr[[1]],3)*100,"%"), size=5, colour="#2171B5")+
  annotate("text", x=time_point_1yr+4, y=surv_prob_1yr[2]+0.02, label=paste0(round(surv_prob_1yr[[2]],3)*100,"%"), size=5, colour="#FFBB44")+
  annotate("text", x=time_point_1yr, y=surv_prob_1yr[1]+0.35, label="1 year", size=5, colour="grey60")+
  annotate("text", x=time_point_5yr+4, y=surv_prob_5yr[1]-0.05, label=paste0(round(surv_prob_5yr[[1]],3)*100,"%"), size=5, colour="#2171B5")+
  annotate("text", x=time_point_5yr+4, y=surv_prob_5yr[2]+0.05, label=paste0(round(surv_prob_5yr[[2]],3)*100,"%"), size=5, colour="#FFBB44")+
  annotate("text", x=time_point_5yr, y=surv_prob_5yr[1]+0.38, label="5 years", size=5, colour="grey60")+
  annotate("text", x=time_point_5yr-20, y = 0.85, label=paste0("HR for progression: ",
                                                               round(PlotsICB$HR, digits = 2),"\n (95% CI: ", round(PlotsICB$lower, digits = 2),
                                                               "-",round(PlotsICB$upper, digits = 2), ") \n Log-rank \n", round(p_value, digits = 4), sep = ""), size=5, fontface="plain")
km

ggsave("/results/Fig5c.pdf", plot = km$plot, width = 8, height = 8)

################################################################################

survivalTBX21 <- survivalTBX21 %>% 
  filter(!is.na(TBX21)) %>% 
  mutate(Sex = sex) %>% 
  mutate(Age = age) %>% 
  mutate(Treatment = Rx) %>% 
  mutate(`BRAF status` = BRAF_status)

TBX21cox <- coxph(Surv(months_death, death_status) ~ Sex + Age + combination + `BRAF status` + cmv + TBX21,
                  data = subset(survivalTBX21))

TBX21summary <- summary(TBX21cox)
TBX21coef <- TBX21summary$coef[, 1]
ciTBX21 <- as.data.frame(TBX21summary$conf.int)

#Create a data frame for plotting

PlotTBX21 <- data.frame(
  covariate = c("Female vs \n Male", "Age", "sICB vs \n cICB"," BRAF Wild-Type \n vs Mutant", "CMV+ vs \n CMV-", "TBX21 expression"),
  HR = ciTBX21[, "exp(coef)"],
  lower = ciTBX21[,"lower .95"],
  upper = ciTBX21[,"upper .95"],
  p_value = TBX21summary$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)

PlotTBX21 <- PlotTBX21 %>% 
  mutate(Significant = ifelse(p_value < 0.05, "Significant", "Not significant")) %>% 
  mutate(Alpha = ifelse(Significant == "Significant", 1, 0.5))

################################################################################
#FIGURE 5d- Survival according to TBX21 expression post cycle 1

x <- ggplot(PlotTBX21, aes(x = HR, y = covariate)) +
  geom_errorbar(data = subset(PlotTBX21, Significant == "Significant"),
                aes(xmax = upper, xmin = lower), 
                width = 0, linetype = "longdash", colour = "#FFBB44", alpha = 1.2, size = 0.5) +
  geom_point(data = subset(PlotTBX21, Significant == "Significant"), size = 5, alpha = 0.8, color = "#FFBB44") +
  geom_errorbar(data = subset(PlotTBX21, Significant == "Not significant"),
                aes(xmax = upper, xmin = lower), 
                width = 0, linetype = "longdash", colour = "#2171B5", alpha = 1.2) +
  geom_point(data = subset(PlotTBX21, Significant == "Not significant"), size = 5, alpha = 0.8, color = "#2171B5") +
  scale_x_log10()+
  geom_vline(xintercept = 1,linetype ="dotted", color = "grey") +
  xlab("Adjusted hazard ratios post cycle 1")+
  ylab("Covariate")+
  theme_classic()+
  theme(axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(hjust=1, size=12),
        axis.text.y = element_text(size=12),
        axis.line = element_line(colour = "grey30", size = 0.2),
        axis.ticks = element_line(color= "grey30", size=0.2), 
        strip.background =element_rect(fill= "grey95"),
        strip.text = element_text(colour = "grey30"),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16))+
  geom_text(aes(label = paste0("p = ", signif(p_value, digits = 2))),
            vjust = -1.5, hjust = 0.5, size = 5, color = "grey30")

ggsave("/results/Fig5d.pdf", plot = x, width = 8, height = 8)

#TBX21: HR= 0.61, p = 0.035
#Age: HR= 1.03, p = 0.023


################################################################################

#Figure 5e- TBX21 expression in tumour bulk RNA-seq

#Requires Riaz et al. data

#TBX21 replication in Riaz data 

# Riaz <- read.table("/data/RiazCountMatrix.txt", sep = "\t", row.names = 1, header = T)
# RiazMeta <- fread('/data/RiazMeta.csv')
# 
# Riaz <- Riaz %>% 
#   filter(!is.na(HUGO) & HUGO != "") %>% 
#   distinct(HUGO, .keep_all = TRUE)
# 
# OnTreatment <- Riaz %>% dplyr::select(ends_with("_On"), ends_with("_Pre"))
# rownames(OnTreatment) <- Riaz$HUGO
# 
# Meta <- as.data.frame(x = colnames(OnTreatment)) 
# Meta <- Meta %>% 
#   mutate(Samples = `colnames(OnTreatment)`) %>% 
#   dplyr::select(-`colnames(OnTreatment)`)
# 
# NormCounts <- baselineNormaliseTransform(files = OnTreatment, samples = Meta)
# Subset <- NormCounts %>% 
#   filter(rownames(NormCounts) %in% c("TBX21")) %>% 
#   as.data.frame()
# SubsetWide <- as.data.frame(t(Subset))
# SubsetWide$Samples <- rownames(SubsetWide)
# Meta <- left_join(Meta, SubsetWide)
# Meta$Patient <- str_split(Meta$Samples, pattern='_', n=2, simplify = T)[,1]
# Meta$cycle <- str_split(Meta$Samples, pattern='_', n=2, simplify = T)[,2]
# OnTreatment <- Meta %>% filter(cycle=='On')
# 
# RiazMeta <- left_join(RiazMeta, OnTreatment)
# 
# RiazMeta <- RiazMeta %>% 
#   filter(Subtype=='CUTANEOUS') %>% 
#   mutate(TBX21 = as.numeric(TBX21))
# 
# RiazMeta$outcome<-ifelse(RiazMeta$Response %in% c('PD'),'Progressive disease','Not progressive disease')
# 
# RiazMeta <- RiazMeta %>%
#   mutate(outcome=factor(outcome, levels = c('Progressive disease','Not progressive disease'))) 
# 
# x <- ggplot(RiazMeta, aes(x = outcome, y = TBX21))+
#   ggbeeswarm::geom_quasirandom(aes(fill = outcome),
#                                colour='grey30',
#                                size=1.5,
#                                pch=21,
#                                alpha=0.6,
#                                width=0.075)+
#   geom_boxplot(colour = "grey30", aes(fill = outcome), alpha=0.5, width = 0.5)+
#   scale_fill_brewer(palette = 13)+
#   theme(legend.position = "none")+
#   xlab("Clinical outcome")+
#   ylab(paste("Normalised TBX21 expression"))+
#   #stat_compare_means(method = "wilcox", label.x.npc = "center", comparisons = list(c("Progressive disease", "Not progressive disease")), size = 5)
#   stat_compare_means(label = "p.signif", size = 7, comparisons = list(c("Progressive disease", "Not progressive disease")))
# 
# ggsave("/results/Fig5e.pdf", plot = x, width = 8, height = 8)

#p= 0.036

################################################################################

##FIGURE VI- EPIDEMIALOGICAL STUDY##

#Rename Melanoma subtypes and identify individuals who have non-melanoma cancers 

SurvivalAllCancer <- clinicalData %>% 
  filter(!is.na(age)) %>% 
  dplyr::select(months_death, months_progression, progression, death_status, cmv, Rx, PtName, age, sex, BRAF_status, subtype, cancer, age_mel, UnknownPrimary) %>% 
  mutate(Cancer = ifelse(cancer == "Melanoma", "Melanoma", "Non-Melanoma")) %>% 
  mutate(intent = ifelse(grepl("Adj", Rx), "adjuvant", "palliative")) %>% 
  mutate(subtype = ifelse(Cancer == "Melanoma" & subtype == "E" | Cancer == "Melanoma" & subtype == "NOS" | Cancer == "Melanoma" & subtype == "small cell", "Cutaneous",
                          ifelse(Cancer == "Melanoma" & subtype == "M", "Mucosal",
                                 ifelse(Cancer == "Melanoma" & subtype == "U", "Uveal", NA))))

#Cohort summary 
#375 patients with serology 

Melanoma <- SurvivalAllCancer %>% 
  filter(PtName != "25") %>% 
  filter(subtype == "Cutaneous") %>% 
  filter(Rx != "RelatNivo" & Rx != "Ipilimumab") %>% 
  filter(cancer == "Melanoma")

wilcox.test(age~ intent, data = Melanoma, conf.int =T) 

# data:  age by intent
# W = 4639.5, p-value = 0.2604
# alternative hypothesis: true location shift is not equal to 0
# 95 percent confidence interval:
#   -6.999942  2.000047
# sample estimates:
#   difference in location 
# -2.999941 

#Prepare cohort data for merging with Biobank data 

AllCancer <- SurvivalAllCancer %>% 
  filter(cancer == "Melanoma") %>% 
  filter(PtName != "25") %>% 
  filter(Rx != "Ipilimumab" & Rx != "RelatNivo") %>% 
  mutate(CMV = ifelse(cmv == "CMV+", 1,
                      ifelse(cmv == 'CMV-', 0, NA))) %>% 
  mutate(status = ifelse(Cancer == "Melanoma" & intent == "adjuvant", "Mel-Adj",
                         ifelse(Cancer == "Melanoma" & intent == "palliative" & BRAF_status == "Wild-type" & subtype == "Cutaneous", "MM-WT", 
                                ifelse(Cancer == "Melanoma" & intent == "palliative" & BRAF_status == "Mutant" & subtype == "Cutaneous", "MM-Mut", NA)))) %>% 
  mutate(Status = ifelse(Cancer == "Melanoma" & intent == "palliative", "Metastatic Melanoma", "Other")) %>% 
  mutate(assessment_centre = "None") %>% 
  filter(subtype == "Cutaneous") %>% 
  filter(!PtName %in% c("193", "197", "300", "408")) %>% #These patients were stage IV resected 
  mutate(sex = ifelse(sex == 1, "Male", "Female")) %>% 
  dplyr::select(age, CMV, status, assessment_centre, Status, sex) %>% 
  mutate(cancer = "Cancer") %>% 
  filter(!is.na(age))

#Prepare UK Biobank data 
#7885 donors 
#4425 CMV+
#3460 CMV-

#For pre-processing Biobank data 

# Biobank <- allmat %>% 
#   mutate(cmv = ifelse(cmv == 2, "CMV+", 
#                       ifelse(cmv == 1, "CMV-", NA)))
# 
# Biobank <- Biobank %>% 
#   mutate(CMV = ifelse(cmv == "CMV+", 1, 0)) %>% 
#   mutate(cancer = "False") %>% 
#   mutate(status = "Healthy") %>% 
#   mutate(Status = "Healthy") %>% 
#   dplyr::select(-cmv) %>% 
#   filter(!is.na(age))
# 
# #Bind Biobank and cancer cohort data and remove patients outside the recruitment age for the Biobank (40-70 years old)
# 
# HealthyCancerCombined <- rbind(AllCancer, Biobank)
# HealthyCancerCombined <- HealthyCancerCombined %>% 
#   filter(age >= 40 & age <= 70) %>% 
#   mutate(HealthyWT = ifelse(status == "Healthy", "Healthy",
#                             ifelse(status == "Mel-Adj", "Adjuvant Melanoma",
#                                    ifelse(status == "MM-WT", "Wild-type Metastatic",
#                                           ifelse(status == "MM-Mut", "Mutant Metastatic",NA)))))

#Compute the OR of CMV seropositivity in cancer groups compared to healthy donors 

AllCancer2 <- AllCancer %>% 
  filter(age >= 40 & age <= 70) %>% 
    mutate(HealthyWT = ifelse(status == "Healthy", "Healthy",
                              ifelse(status == "Mel-Adj", "Adjuvant Melanoma",
                                     ifelse(status == "MM-WT", "Wild-type Metastatic",
                                            ifelse(status == "MM-Mut", "Mutant Metastatic",NA))))) %>% 
  group_by(CMV, HealthyWT) %>% 
  summarise(Freq = n()) %>%
  ungroup() %>% 
  filter(!is.na(CMV)) %>% 
  mutate(cmv = ifelse(CMV == 1, "CMV+", "CMV-")) %>% 
  dplyr::select(cmv, HealthyWT, Freq, -CMV)

Biobank$HealthyWT <- "Healthy"

Biobank2 <- Biobank %>% 
  dplyr::select(-centre) %>% 
  group_by(cmv, HealthyWT) %>% 
  summarise(Freq = sum(Freq))

HealthyCancerCombined <- rbind(Biobank2, AllCancer2)

#Metastatic BRAF wild-type 

HealthyMMwildtype <- HealthyCancerCombined %>% 
  filter(HealthyWT != "Adjuvant Melanoma" & HealthyWT != "Mutant Metastatic" & !is.na(HealthyWT))

table_matrix <- xtabs(Freq ~ cmv + HealthyWT, data = HealthyMMwildtype)

#fisherMMwildtype <- fisher.test(HealthyMMwildtype$cmv, HealthyMMwildtype$HealthyWT)
fisherMMwildtype <- fisher.test(table_matrix)
oddsMMwildtype <- fisherMMwildtype$estimate
confIntMMwildtype <- fisherMMwildtype$conf.int
pvalMMwildtype <- fisherMMwildtype$p.value

#Adjuvant melanoma 

HealthyAdj <- HealthyCancerCombined %>% 
  filter(HealthyWT != "Wild-type Metastatic" & HealthyWT != "Mutant Metastatic" & !is.na(HealthyWT))

table_matrix <- xtabs(Freq ~ cmv + HealthyWT, data = HealthyAdj)

#fisherAdj <- fisher.test(HealthyAdj$CMV, HealthyAdj$HealthyWT)
fisherAdj <- fisher.test(table_matrix)
oddsAdj <- fisherAdj$estimate
confIntAdj <- fisherAdj$conf.int
pvalAdj <- fisherAdj$p.value

#Metastatic BRAF mutant 

HealthyMMmutant <- HealthyCancerCombined %>% 
  filter(HealthyWT != "Adjuvant Melanoma" & HealthyWT != "Wild-type Metastatic" & !is.na(HealthyWT))

table_matrix <- xtabs(Freq ~ cmv + HealthyWT, data = HealthyMMmutant)

#fisherMMmutant <- fisher.test(HealthyMMmutant$CMV, HealthyMMmutant$HealthyWT)
fisherMMmutant <- fisher.test(table_matrix)
oddsMMmutant <- fisherMMmutant$estimate
confIntMMmutant <- fisherMMmutant$conf.int
pvalMMmutant <- fisherMMmutant$p.value

AllMetaMel <- HealthyCancerCombined
AllMetaMel$Status <- "Metastatic"
AllMetaMel <- AllMetaMel %>% 
  mutate(Status = ifelse(HealthyWT == "Adjuvant Melanoma", "Adj",
                         ifelse(HealthyWT == "Healthy", "Healthy", Status))) %>% 
  mutate(Status = ifelse(is.na(Status), "Metastatic", Status)) %>% 
  filter(Status != "Adj")

table_matrix <- xtabs(Freq ~ cmv + Status, data = AllMetaMel)

#fisherMel <- fisher.test(AllMetaMel$CMV, AllMetaMel$Status)
fisherMel <- fisher.test(table_matrix)
oddsMel <- fisherMel$estimate
confIntMel <- fisherMel$conf.int
pvalMel <- fisherMel$p.value

#Create a dataframe for plotting results 

results <- data.frame(
  Comparison = c("Stage II/III resectable","BRAF Mutant \n Metastatic", "BRAF Wild-type \n Metastatic",  "Stage IV/Unresectable Stage III-Metastatic"),
  Odds_Ratio = c(oddsAdj, oddsMMmutant, oddsMMwildtype, oddsMel),
  Lower_CI = c(confIntAdj[1], confIntMMmutant[1], confIntMMwildtype[1], confIntMel[1]),
  Upper_CI = c(confIntAdj[2], confIntMMmutant[2], confIntMMwildtype[2], confIntMel[2]),
  pval = c(pvalAdj, pvalMMmutant, pvalMMwildtype, pvalMel)
)

results <- results %>% 
  mutate(Significant = ifelse(pval < 0.05, "Significant", "Not significant")) %>% 
  mutate(Alpha = ifelse(Significant == "Significant", 1, 0.5)) %>% 
  mutate(Comparison = as.factor(Comparison))

################################################################################
#FIGURE 6a- OR of CMV in melanoma patients compared to healthy donors

results <- results %>% 
  mutate(Significant = ifelse(pval < 0.05, "Significant", "Not significant")) %>% 
  mutate(Alpha = ifelse(Significant == "Significant", 1, 0.5)) 

results$Comparison <- factor(results$Comparison, levels = c("Stage II/III resectable", "BRAF Mutant \n Metastatic", "BRAF Wild-type \n Metastatic", "Stage IV/Unresectable Stage III-Metastatic"))

x <- ggplot(results, aes(x = `Odds_Ratio`, y = Comparison)) +
  geom_errorbar(data = subset(results, Significant == "Significant"),
                aes(xmax = `Upper_CI`, xmin = `Lower_CI`), 
                width = 0, linetype = "longdash", colour = "#EF6548", alpha = 1.2, size = 0.5) +
  geom_point(data = subset(results, Significant == "Significant"), size = 5, alpha = 0.8, color = "#EF6548") +
  geom_errorbar(data = subset(results, Significant == "Not significant"),
                aes(xmax = `Upper_CI`, xmin = `Lower_CI`), 
                width = 0, linetype = "longdash", colour = "#9E9AC8", alpha = 1.2) +
  geom_point(data = subset(results, Significant == "Not significant"), size = 5, alpha = 0.8, color = "#9E9AC8") +
  scale_x_log10()+
  geom_vline(xintercept = 1,linetype ="dotted", color = "grey") +
  xlab("Odds ratios of CMV seropositivity \n compared to UKB donors")+
  ylab("Cutaneous Melanoma")+
  theme_classic()+
  theme(axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(hjust=1, size=12),
        axis.text.y = element_text(size=12),
        axis.line = element_line(colour = "grey30", size = 0.2),
        axis.ticks = element_line(color= "grey30", size=0.2), 
        strip.background =element_rect(fill= "grey95"),
        strip.text = element_text(colour = "grey30"),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16))+
  geom_text(aes(label = paste0("p = ", signif(pval, digits = 2))),
            vjust = -1, hjust = 0.5, size = 5, color = "grey30")

ggsave("/results/Fig6a.pdf", plot = x, width = 8, height = 8)

#BRAF mutant: 0.44, p = 0.0018
#BRAF WT: 0.58, p = 0.026
#MM: 0.52, p = 0.00018
#Adjuvant Melanoma: 1.2, p = 0.83

################################################################################

#Requires raw UK biobank data 

#Check validity by age and sex matching 

# MetastaticMelanoma <- SurvivalAllCancer %>% 
#   filter(PtName != "25") %>% 
#   filter(Rx != "Ipilimumab" & Rx != "RelatNivo") %>% 
#   mutate(CMV = ifelse(cmv == "CMV+", 1,
#                                           ifelse(cmv == 'CMV-', 0, NA))) %>% 
#   mutate(status = ifelse(Cancer == "Non-Melanoma", "Non-melanoma",
#                          ifelse(Cancer == "Melanoma" & intent == "adjuvant", "Mel-Adj",
#                                 ifelse(Cancer == "Melanoma" & intent == "palliative", "MM-Mel", NA)))) %>% 
#   mutate(sex = ifelse(sex == 1, "Male", "Female")) %>% 
#   mutate(assessment_centre = "None") %>% 
#   dplyr::select(age, CMV, status, assessment_centre, sex) %>% 
#   filter(age >= 40 & age <= 70) %>% 
#   filter(status == "MM-Mel") %>% 
#   mutate(cancer = "Cancer") %>% 
#   filter(!is.na(age) & !is.na(CMV))
# 
# Biobank2 <- Biobank %>% 
#   dplyr::select(age, CMV, status, assessment_centre, sex, cancer)
# AgeMatched <- rbind(Biobank2, MetastaticMelanoma)
# AgeMatched <- AgeMatched %>% 
#   mutate(Status = ifelse(status == "Healthy", 0, 1))
# 
# #Perform age and sex matching with multiple iterations are there will be more than one match for each patient in the biobank
# 
# Iterations <- 1000
# 
# #Store p-values from each iteration
# 
# Pvals <- numeric(Iterations)
# 
# #Create lists to store outputs
# 
# MatchedData <- vector("list", Iterations)
# PermutationResults <- vector("list", Iterations)
# 
# #Perform the repeated matching and permutation test
# 
# set.seed(123)
# 
# for (i in 1:Iterations) {
#   
#   #Use MatchIt to store match Biobank-patient matching according to age and sex 
#   
#   MatchOutput <- matchit(Status ~ age + sex, data = AgeMatched, method = "nearest", ratio = 1)
#   
#   #Extract matched data
#   
#   MatchedData[[i]] <- match.data(MatchOutput)
#   
#   MatchedData[[i]]$CMV <- as.factor(MatchedData[[i]]$CMV)
#   
#   #Perform Fisher-Pitman permutation test
#   
#   PermutationResults[[i]] <- independence_test(CMV ~ Status, data = MatchedData[[i]], distribution = approximate(nresample = 10000))
# }
# 
# #Aggregate the results by taking the median p-value over 1000 iterations
# 
# Pvals <- sapply(PermutationResults, function(x) pvalue(x))
# 
# FinalResult <- median(Pvals)
# 
# #p = 0.0053
# #Range 0.0036 0.0078

################################################################################
#FIGURE 6b- Age according to CMV status in Melanoma and Non-Melanoma 

#Median CMV+ 71 years 
#Median CMV- 64 years 
# data:  age by cmv
# W = 6054.5, p-value = 0.001404
# alternative hypothesis: true location shift is not equal to 0
# 95 percent confidence interval:
#   -9.000034 -2.000055
# sample estimates:
#   difference in location 
# -5.999971 

# SurvivalAllCancer %>% 
#   filter(!is.na(cmv) & intent == "palliative" &  !Rx %in% c("CAPOX", "FOLFOX", "FOLFIRI", "FOLFOX+panitumumab")) %>% 
#   dplyr::mutate(sex = ifelse(sex == 1, "Male", "Female")) %>% 
#   mutate(across(where(is.character), as.factor)) %>% 
#   dplyr::select(age, sex, cmv, Cancer) %>% 
#   gtsummary::tbl_summary(by = cmv) %>% 
#   gtsummary::add_p() %>% 
#   modify_spanning_header(all_stat_cols() ~ "**Melanoma/Non-Melanoma Epidemiology**")

SurvivalAllCancer2 <- SurvivalAllCancer %>% 
  filter(cancer == "Melanoma" & subtype == "Cutaneous" | Cancer == "Non-Melanoma") %>% 
  filter(PtName != "25") %>% 
  filter(Rx != "Ipilimumab" & Rx != "RelatNivo") 

x <- ggplot(subset(SurvivalAllCancer2, !is.na(cmv) & intent == "palliative" &  !Rx %in% c("CAPOX", "FOLFOX", "FOLFIRI", "FOLFOX+panitumumab")), aes(x = cmv, y = age))+
  ggbeeswarm::geom_quasirandom(aes(fill = cmv),
                               colour='grey30',
                               size=1.5,
                               pch=21,
                               alpha=0.6,
                               width=0.075)+   
  geom_boxplot(colour = "grey30", aes(fill = cmv), alpha=0.5, width = 0.5)+
  scale_fill_brewer(palette = 7)+
  facet_wrap(~Cancer)+
  theme(legend.position = "none")+
  xlab("CMV status")+
  ylab(paste("age at initiation of \n systemic treatment (years)"))+
  labs(fill = "CMV status", color = "CMV status")+
  #stat_compare_means(method = "wilcox", label.x.npc = "center", comparisons = list(c("CMV+", "CMV-")), size = 5)
  stat_compare_means(comparisons = list(c("CMV+", "CMV-")), label = "p.signif", size = 7)

ggsave("/results/Fig6b.pdf", plot = x, width = 8, height = 8)

################################################################################

#Focus on BRAF mutational status 

BRAFMMdata <- SurvivalAllCancer %>% 
  filter(PtName != "25") %>% 
  filter(Rx != "Ipilimumab" & Rx != "RelatNivo") %>% 
  filter(cancer == "Melanoma") %>% 
  filter(!is.na(cmv) & !is.na(BRAF_status)) %>% 
  filter(!grepl("Adj", Rx)) %>% 
  filter(subtype == "Cutaneous")

# BRAFMMdata %>% 
#   dplyr::mutate(sex = ifelse(sex == 1, "Male", "Female")) %>% 
#   mutate(across(where(is.character), as.factor)) %>% 
#   dplyr::select(age, sex, BRAF_status, cmv) %>% 
#   gtsummary::tbl_summary(by = cmv) %>% 
#   gtsummary::add_p() %>% 
#   modify_spanning_header(all_stat_cols() ~ "**BRAF epidemiology data**")s

#No significant age difference between CMV negative BRAF Mutant and Wild-type

CMVnegBRAF <- BRAFMMdata %>% 
  filter(cmv == "CMV-")

wilcox.test(age ~ BRAF_status, conf.int = T, data = CMVnegBRAF)

# data:  age by BRAF_status
# W = 2008, p-value = 0.4735
# alternative hypothesis: true location shift is not equal to 0
# 95 percent confidence interval:
#   -7.00000  2.99995
# sample estimates:
#   difference in location 
# -1.999961 

wilcox.test(age ~ BRAF_status, conf.int = T, data = BRAFMMdata)

# data:  age by BRAF_status
# W = 4926, p-value = 6.91e-05
# alternative hypothesis: true location shift is not equal to 0
# 95 percent confidence interval:
#   -11.000013  -4.000014
# sample estimates:
#   difference in location 
# -7.999958 


################################################################################
#FIGURE 6c- BRAF status age according to CMV serostatus

#Median CMV+ = 74
#Median CMV- = 65
# data:  age by cmv
# W = 1839, p-value = 0.0001294
# alternative hypothesis: true location shift is not equal to 0
# 95 percent confidence interval:
#   -12.999952  -4.000002
# sample estimates:
#   difference in location 
# -8.999977 

x <- ggplot(subset(BRAFMMdata), aes(x = cmv, y = age))+
  ggbeeswarm::geom_quasirandom(aes(fill = cmv),
                               colour='grey30',
                               size=1.5,
                               pch=21,
                               alpha=0.6,
                               width=0.075)+   
  geom_boxplot(colour = "grey30", aes(fill = cmv), alpha=0.5, width = 0.5)+
  scale_fill_brewer(palette = 7)+
  facet_wrap(~BRAF_status, nrow = 1)+
  theme(legend.position = "none")+
  xlab("CMV status")+
  ylab(paste("age at initiation of \n systemic treatment (years)"))+
  labs(fill = "CMV status", color = "CMV status")+
  #stat_compare_means(method = "wilcox", label.x.npc = "center", comparisons = list(c("CMV+", "CMV-")), size = 5)
  stat_compare_means(comparisons = list(c("CMV+", "CMV-")), label = "p.signif", size = 7)

ggsave("/results/Fig6c.pdf", plot = x, width = 8, height = 8)

#Interaction test for BRAF mutational status and CMV status 

model <- lm(age ~ sex + cmv + BRAF_status + cmv:BRAF_status, data = BRAFMMdata)
summary(model)

#p = 0.0029

################################################################################
#FIGURE 6d- Depletion of CMV in BRAF mutants: Cutaneous Metatstatic melanoma only 

PlotBRAF <- as.data.frame(table(cmv = BRAFMMdata$cmv, BRAF = BRAFMMdata$BRAF_status))
PlotBRAF <- PlotBRAF %>% 
  mutate(`CMV status` = cmv)
fisherTest <- fisher.test(BRAFMMdata$BRAF_status, BRAFMMdata$cmv)
OR <- fisherTest$estimate
P <- fisherTest$p.value

Summary <- BRAFMMdata %>%
  group_by(cmv, BRAF_status) %>%
  summarise(count = n()) %>%
  ungroup()

p <- ggplot(Summary, aes(x = BRAF_status, y = count, fill = BRAF_status)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, color = "grey45") +
  facet_wrap(~cmv) +
  labs(
    x = "BRAF status",
    y = "Number of patients") +
  scale_fill_brewer(palette = 5)+
  theme(legend.position = "None")

annotation_data <- data.frame(
  cmv = "CMV+",
  BRAF_status = "Mutant",  
  count = 20,
  label = paste0("OR = ", round(OR, 1), "\nP = ", round(P, 4))
)

p <- p + geom_text(data = annotation_data, aes(label = label), 
                   x = 1, y = 40,
                   size = 5, hjust = 0.5)
p

ggsave("/results/Fig6d.pdf", plot = p, width = 8, height = 8)

# Mutant Wild-type
# CMV-     61        71
# CMV+     31        82

# data:  BRAFMMdata$cmv and BRAFMMdata$BRAF_status
# p-value = 0.005431
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.235316 3.851759
# sample estimates:
#   odds ratio 
# 2.167779 

################################################################################

# ##END##-------------------------------------------------------------------------
