##SCRIPT FINAL-CMV SUPPLEMENTARY##
#Authors: Gusztav Milotay, Martin Little, Robert Watson & Benjamin Fairfax
#Last update: 15/01/2025
#Fully annotated: Yes
#Reviewed: Yes

##START-------------------------------------------------------------------------
##LOAD PACKAGES-----------------------------------------------------------------

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
library(org.Hs.eg.db) #package version is 3.16.0 due to R version on code ocean being 4.2- this impacts TF analysis outcomes 
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

cmv <- read.delim("/data/CMVdata.txt",sep = "\t")

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

tox <- tox[1:772,]

#Haematological data 

bloods <- read.csv("/data/HaemData.csv")

#UK Biobank data 

Biobank <- load("/data/UKB.Rd")

#Single cell data for CD8+ T cells  

cd8 <- readRDS("/data/scData.rds")

#Gene list used for calculation of cytotoxicity score- same as used in Watson et al. 2021 

cytotoxicity <- read.csv("/data/BulkScores.csv")

clinicalData <- clinicalData %>%
  mutate(PtName = as.character(ID))%>%
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

#Up to date clinical data- last patient with sufficient followup for survival analysis- defined as having at least 6 months of followup if the patient has not yet progressed 
#Censor at 90 months 

surv_last_update <- survival_mm %>% 
  filter(progression == 1 | months_progression > 6) %>% 
  mutate(months_progression = ifelse(months_progression > 90, 90, months_progression)) %>% 
  mutate(months_death = ifelse(months_death > 90, 90, months_death))

#Pembrolizumab and Nivolumab will be treated as anti-PD1 sICB

surv_last_update <- surv_last_update %>% 
  mutate(combination = ifelse(Rx == "IpiNivo", "cICB", "sICB"))

#Set plotting theme 

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

#Extracting normalised counts from DESEQ2

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

#Geometric mean calculator for Cytotoxicity and Mitotic scores

GeometricMean <- function(x) {
  x <- x[x > 0]  
  if (length(x) == 0) {
    return(NA) 
  }
  exp(mean(log(x)))  
  
}

################################################################################
#EXTENDED FIGURE I- NLR IS DETRIMENTAL TO SURVIVAL ON IMMUNOTHERAPY##
################################################################################

#Join bloods data with metadata containing CMV status 

bloods <- bloods %>% 
  mutate(patient = ID) 
NLR <- left_join(surv_last_update, bloods)

#Perform a Cox regression for overall survival, correcting for age, sex, treatment, and BRAF status

NLRcox <- coxph(Surv(months_death, death_status) ~ age + sex + combination + `BRAF_status` + `NLR`,
                  data = subset(NLR, Rx != "Ipilimumab" & Rx != "RelatNivo"))
NLRsummary <- summary(NLRcox)
coefsNLR <- NLRsummary$coef[, 1]
ciNLR <- as.data.frame(NLRsummary$conf.int)

#Create a data frame for plotting

PlotNLR <- data.frame(
  covariate = c("Age", "Female vs \n Male", "sICB vs \n cICB", "BRAF Wild-Type \n vs Mutant", "NLR"),
  HR = ciNLR[, "exp(coef)"],
  lower = ciNLR[,"lower .95"],
  upper = ciNLR[,"upper .95"],
  p_value = NLRsummary$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)

PlotNLR <- PlotNLR %>% 
  mutate(Significant = ifelse(p_value < 0.05, "Significant", "Not significant")) %>% 
  mutate(Alpha = ifelse(Significant == "Significant", 1, 0.5))

################################################################################
#FIGURE S1- Multivariate survival analysis incorporating NLR  

x <- ggplot(PlotNLR, aes(x = HR, y = covariate)) +
  geom_errorbar(data = subset(PlotNLR, Significant == "Significant"),
                aes(xmax = upper, xmin = lower), 
                width = 0, linetype = "longdash", colour = "#FFBB44", alpha = 1.2, size = 0.5) +
  geom_point(data = subset(PlotNLR, Significant == "Significant"), size = 5, alpha = 0.8, color = "#FFBB44") +
  geom_errorbar(data = subset(PlotNLR, Significant == "Not significant"),
                aes(xmax = upper, xmin = lower), 
                width = 0, linetype = "longdash", colour = "#2171B5", alpha = 1.2)+
  geom_point(data = subset(PlotNLR, Significant == "Not significant"), size = 5, alpha = 0.8, color = "#2171B5") +
  scale_x_log10()+
  geom_vline(xintercept = 1,linetype ="dotted", color = "grey") +
  xlab("Adjusted hazard ratios")+
  ylab("Covariate")+
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
  geom_text(aes(label = paste0("p = ", signif(p_value, digits = 2))),
            vjust = -1.5, hjust = 0.5, size = 5, color = "grey30")

  ggsave("/results/FigS1.pdf", plot = x, width = 8, height = 8)

################################################################################
##EXTENDED FIGURE II- BULK SCORES##
################################################################################

#Determine a CMV score- this will be a PCA analysis of top up and downregulated genes identified at baseline and subsequently applied to baseline and post cycle 1 samples 

#Filter metastatic melanoma data for C1 expression and create a new sample variable 

metadata <- mmdata %>% 
  filter(patient != 25 & Rx != "Ipilimumab" & sample != "260_C2") %>% 
  filter(cycle == "C1") %>% 
  mutate(sample = paste0(sample, "_CD8_None")) %>% 
  mutate(sample = as.factor(sample)) 

#For the batch_data, select CD8+ T cell samples taken at baseline (format C1_CD8_None)

batch <- batch_data %>% 
  mutate(patient = as.numeric(patient)) %>% 
  mutate(sample = as.factor(sample)) %>% 
  filter(grepl(paste("C1_CD8_None", sep = "|"), sample))

#Merge batch data with metadata which can be used in DESeq2 for batch correction 
#merge using patient, cycle, and sample and discard anything that does not match (unmatched = "drop")

MetadataBatch <- left_join(metadata, batch, by = c("patient", "cycle", "sample") , unmatched = "drop", keep = F)

#Latest samples have not yet been assigned a batch 

NewSamples <- c(262,264,271, 274, 275, 276, 279, 281, 286, 291, 297, 302,265, 269, 273, 284, 288, 292)
MissingBatches2022.Sep <- c(233, 239, 243, 254, 261)
MissingBatchesMay2021 <- c(224)
MissingBatchesDec2020 <- c(161)
MissingBeforeSep2019 <- c(2, 26, 41, 44, 45, 81, 86)
MissingSep2019 <- c(153, 155)

FilteredMeta <- MetadataBatch %>%
  mutate_if(is.character, as.factor) %>% 
  mutate(cmv = ifelse(cmv == "CMV+", "CMVpositive",
                      ifelse(cmv == "CMV-", "CMVnegative", NA))) %>%
  mutate(origBatch = case_when(
    patient %in% NewSamples ~ "Jan2024",
    TRUE ~ origBatch)) %>% 
  mutate(origBatch = case_when(patient %in% MissingBatches2022.Sep ~ "2022.Sep", TRUE ~ origBatch)) %>% 
  mutate(origBatch = case_when(patient %in% MissingBatchesMay2021 ~ "May_2021", TRUE ~ origBatch)) %>% 
  mutate(origBatch = case_when(patient %in% MissingBatchesDec2020 ~ "Dec_2020", TRUE ~ origBatch)) %>% 
  mutate(origBatch = case_when(patient %in% MissingBeforeSep2019 ~ "Before_SEP2019", TRUE ~ origBatch)) %>% 
  mutate(origBatch = case_when(patient %in% MissingSep2019 ~ "Sep_19", TRUE ~ origBatch)) %>% 
  mutate(origBatch = ifelse(patient > 302, "June_2024", origBatch)) %>% 
  filter(!is.na(cmv))

#Ensure the rownames still correspond to the sample ids

rownames(FilteredMeta) <- FilteredMeta$sample

#Subset counts data for baseline CD8+ T cell samples (C1_CD8_None)
#Add Ensembl ids as rownames in this new dataset

rownames(counts) <- counts$Gene
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
#Select metadata columns of interest: sample, sex, cycle, batch, age as continuous variable and cmv status- remove columns that are not required 

CountsReordered <- SubsetCounts[, match(rownames(FilteredMeta), colnames(SubsetCounts))]
all(colnames(CountsReordered) == rownames(FilteredMeta))
rownames(CountsReordered) = rownames(counts) 
CountsReordered$id <- rownames(CountsReordered)
FinalMeta <- FilteredMeta %>% 
  dplyr::select(sample, patient, sex, age, cmv, origBatch, Rx) %>% 
  mutate(sex = as.factor(sex)) %>% 
  mutate(age = scale(age, center = T))
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

#206 C1 CD8 T cell bulk samples with CMV serology 
#111 CMV- and 95 CMV+

#Create a DESEQ object with age, sex and batch as covariates  

dds <- DESeqDataSetFromMatrix(countData = CountsFinal, 
                              colData = FinalMeta, 
                              design = ~ age  + sex + origBatch + cmv)

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

#7576 differentially expressed genes 

#Filter metadata for baseline and C2 samples 

metadata <- mmdata %>% 
  filter(patient != 25 & Rx != "Ipilimumab" & sample != "260_C2") %>% 
  filter(cycle == "C1" | cycle == "C2") %>% 
  mutate(sample = paste0(sample, "_CD8_None")) %>% 
  mutate(sample = as.factor(sample)) 

#For the batch_data, select cells of interest in the format 'CX_CD8' referring to cycle and cell marker

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
#Select metadata columns of interest: sample, sex, cycle, batch, age as continuous variable, cmv sereostatus, etc...

SubsetReordered <- SubsetCounts[, match(rownames(FilteredMeta), colnames(SubsetCounts))]
all(colnames(SubsetReordered) == rownames(FilteredMeta))
rownames(SubsetReordered) <- rownames(counts) 
SubsetReordered$id <- rownames(SubsetReordered)
FinalMeta <- FilteredMeta %>% 
  dplyr::select(sample, sex, cycle, age, cmv, Rx, patient)

#Filter out genes which have low reads as these could easily be classed as differentially expressed. 10 is the count threshold suggested by the vignette.
#Do a final check that the samples in the counts data are matching and in the same order as the samples in the metadata 

keep <- rowMeans(SubsetReordered[,1:(ncol(SubsetReordered)-1)]) >= 10
SubsetReordered <- SubsetReordered[keep, ] 
rownames(SubsetReordered) <- SubsetReordered$id
CountsFinal <- dplyr::select(SubsetReordered, -id)
rownames(CountsFinal) <- rownames(SubsetReordered)
rownames(FinalMeta) <- FinalMeta$sample
all(colnames(CountsFinal) == rownames(FinalMeta))

#417 CD8 samples from C1 and C2 
#236 C1 and 181 C2

#Identify a threshold of significance for the CMV score based on a ROC analysis

genesPval <- ResultsOrdered$padj

#Thresholds where selected from all significant genes (7175) down to 44 significant genes 

thresholds <- 10^seq(log10(0.05), -15, by = -0.5) 

#Create an empty list for ROC results 

rocResults <- list()

#Store AUC for each threshold 

AUCresults <- data.frame(threshold = numeric(0), auc = numeric(0))

#Extract normalised counts of genes from principle component analysis- PC1 of top up and downregulated genes will be used as the CMV score 

NormCounts <- baselineNormaliseTransform(files = CountsFinal, samples = FinalMeta)
NormCounts$symbol <- mapIds(org.Hs.eg.db, keys = rownames(NormCounts), keytype = "ENSEMBL", column = "SYMBOL")
test <- FinalMeta
rownames(test) <- FinalMeta$sample

#Loop through each threshold 

for (threshold in thresholds) {
  
  #Select genes that meet the threshold for significance 
  
  selectedGenes <- rownames(ResultsOrdered)[which(genesPval <= threshold & !is.na(genesPval))]
  length = length(selectedGenes)
 
  if (length(selectedGenes) > 1) {
    CMVinduced <- NormCounts %>% 
      dplyr::select(-c(X, symbol)) %>% 
      filter(rownames(NormCounts) %in% selectedGenes) 
    CMVinducedWide <- t(CMVinduced)
    rownames(CMVinducedWide) <- rownames(test)
    CMVinducedWide <- as.data.frame(CMVinducedWide)
    CMVinducedWide$sample <- rownames(CMVinducedWide)
    CountsMeta <- left_join(test, CMVinducedWide)
    
    #Perform PCA analysis to create a score for each individual on the basis of CMV-induced genes- this will be PC1
    
    GeneCounts <- CMVinducedWide %>% 
      dplyr::select(-sample)
    PCA <- prcomp(x = GeneCounts)
    PCAdata <- as.data.frame(PCA$x)
    CMVscore <- PCAdata %>% 
      dplyr::select(PC1)
    CMVscore$sample <- rownames(CMVscore)
    PCAmeta <- left_join(test, CMVscore)
    ROCcurve <- roc(PCAmeta$cmv, PCAmeta$PC1, levels = c("CMVnegative", "CMVpositive"))
    rocResults[[as.character(threshold)]] <- ROCcurve
    AUCvalue <- auc(ROCcurve)
    AUCresults <- rbind(AUCresults, data.frame(threshold = threshold, auc = AUCvalue))
  }
}

MaxAUC <- max(AUCresults$auc)

#Determine the optimal threshold 

BestThreshold <- AUCresults$threshold[which.max(AUCresults$auc)]
AUCresults <- AUCresults %>% 
  mutate(ChosenThreshold = ifelse(threshold == BestThreshold, "Max", "Other"))

################################################################################
#FIGURE S2a- ROC analysis to determine optimal threshold for CMV score 

xlab <- expression(-log[10]~ (adjP))
ylab <- expression(AUC)

x <- ggplot(AUCresults, aes(x = -log10(threshold), y = auc))+
  geom_point(aes(color = ChosenThreshold), alpha = 0.6, size = 2)+
  scale_color_manual(values=c("#EB7926", "#005A32"))+
  xlab(xlab)+
  ylab(ylab)+
  geom_vline(aes(xintercept = -log10(BestThreshold)), linetype = "dotted", color = "grey")+
  geom_hline(aes(yintercept = MaxAUC), linetype = "dotted", color = "grey")+
  theme_classic()+
  theme(legend.position = "None")+
  annotate("text", x = -log10(5e-13),
           y = MaxAUC,
           label = paste0("Max AUC: ", round(MaxAUC, 2), "\n Threshold: ",round(BestThreshold, 15)),
           vjust = 3) +
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

ggsave("/results/FigS2a.pdf", plot = x, width = 8, height = 8)
################################################################################

#Filter metadata for baseline and C2 samples 

metadata <- mmdata %>% 
  filter(patient != 25 & Rx != "Ipilimumab" & sample != "260_C2") %>%  
  filter(cycle == "C1" | cycle == "C2") %>% 
  mutate(sample = paste0(sample, "_CD8_None")) %>% 
  mutate(sample = as.factor(sample)) 

#For the batch_data, select cells of interest in the format 'CX_CD8' referring to cycle and cell marker

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
#Select metadata columns of interest: sample, sex, cycle, batch, age as continuous variable, cmv sereostatus, etc...

SubsetReordered <- SubsetCounts[, match(rownames(FilteredMeta), colnames(SubsetCounts))]
all(colnames(SubsetReordered) == rownames(FilteredMeta))
rownames(SubsetReordered) <- rownames(counts) 
SubsetReordered$id <- rownames(SubsetReordered)
FinalMeta <- FilteredMeta %>% 
  dplyr::select(sample, sex, cycle, age, cmv, Rx, patient)

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
NormCounts <- baselineNormaliseTransform(files = CountsFinal, samples = FinalMeta)
NormCounts$symbol <- mapIds(org.Hs.eg.db, keys = rownames(NormCounts), keytype = "ENSEMBL", column = "SYMBOL")
CMVinduced <- NormCounts %>% 
  dplyr::select(-X) %>% 
  filter(rownames(NormCounts) %in% rownames(CMVgenes)) %>% 
  dplyr::select(-c(symbol))
CMVinducedWide <- t(CMVinduced)

#Add genes to metadata 

rownames(CMVinducedWide) <- rownames(FinalMeta)
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
CMVscore$sample <- FinalMeta$sample
PCAmeta <- left_join(FinalMeta, CMVscore)

################################################################################
#FIGURE S2b- PC1 of CMV genes according to CMV status at baseline

PlotData <- PCAmeta %>% 
  filter(cycle == "C1") %>% 
  mutate(cmv = ifelse(cmv == "CMVpositive", "CMV+", "CMV-")) %>% 
  filter(!is.na(cmv)) %>% 
  mutate(over68 = ifelse(age >= 68, "Over", "Under"))

x <- ggplot(PlotData, aes(x = cmv, y = PC1))+
  ggbeeswarm::geom_quasirandom(aes(fill = cmv),
                               colour='grey30',
                               size=1.5,
                               pch=21,
                               alpha=0.6,
                               width=0.075)+   
  geom_boxplot(colour = "grey30", aes(fill = cmv), alpha=0.5, width = 0.5)+
  scale_fill_brewer(palette = 7)+
  xlab("CMV status")+
  ylab("Baseline CMV score")+
  theme(legend.position = "none")+
  labs(fill = "CMV status", color = "CMV status")+
  #stat_compare_means(method = "wilcox", label.x.npc = "center", comparisons = list(c("CMV+", "CMV-")), size = 5)
  stat_compare_means(comparisons = list(c("CMV+", "CMV-")), label = "p.signif", size = 7)

#p < 2.2e-16

ggsave("/results/FigS2b.pdf", plot = x, width = 8, height = 8)

################################################################################

#Determine total TRB read counts for each individual 

#Large clones were previously associated with 

cloneCounts <- tcr_processed %>%
  filter(grepl(paste("TRB", sep = "|"), bestVHit)) %>% 
  group_by(ID) %>%
  summarise(total_cell_count = sum(cloneCount))

#Filter for C2 beta chains

BetaSizeC2 <- tcr_processed %>% 
  filter(grepl(paste("TRB", sep = "|"), bestVHit)) %>% 
  filter(grepl("_C2_", ID))

#Join total counts with filtered counts 

BetaSizeC2 <- left_join(BetaSizeC2, cloneCounts)

#Determine the number and proportion of clones by size per individual at day 21

BetaSizeC2 <- BetaSizeC2 %>% 
  mutate(repertoire_occupancy = cloneCount/total_cell_count * 100) %>% 
  mutate(size = ifelse(repertoire_occupancy > 0.5, "large", 
                       ifelse(repertoire_occupancy <= 0.5, "intermediate", NA))) %>% 
  group_by(sample) %>% 
  summarise(`> 0.5%` = sum(grepl("large", size)),
            `< 0.5%` = sum(grepl("intermediate", size)))

#join with metadata 

BetaSizeC2 <- left_join(BetaSizeC2, mmdata)

#Pivot each size bin into a single variable 

BetaSizeC2long <- BetaSizeC2 %>% 
  filter(sample != "250_C2") %>% 
  filter(!grepl("Adj", Rx)) %>% 
  filter(cancer == "Melanoma") %>% 
  mutate(combination = ifelse(Rx == "IpiNivo", "cICB", "sICB")) %>% 
  pivot_longer(names_to = "size", cols = c(`> 0.5%`, `< 0.5%`), values_to = "size_counts")

#Relevel factor 

BetaSizeC2long$size <- factor(BetaSizeC2long$size, levels = c("< 0.5%", "> 0.5%"))

#Join with CMV score 

BetaSizeC2long <- BetaSizeC2long %>% 
  dplyr::select(size, size_counts, cmv, patient, sample) %>% 
  mutate(sample = paste0(sample, "_CD8_None"))

CMVscoreClones <- PCAmeta %>% 
  filter(cycle == "C2") %>% 
  dplyr::select(patient, PC1, sample, age, sex, Rx)

BetaSize <- left_join(BetaSizeC2long, CMVscoreClones)

################################################################################
#FIGURE S2c- Number of clones by size according to CMV score post cycle 1

# BetaSize %>% 
#   filter(!is.na(cmv)) %>% 
#   filter(!is.na(PC1)) %>% 
#   mutate(Treatment = ifelse(Rx == "IpiNivo", "cICB", "sICB")) %>% 
#   dplyr::mutate(sex = ifelse(sex == 1, "Male", "Female")) %>% 
#   mutate(across(where(is.character), as.factor)) %>% 
#   select(age, sex, BRAF_status, cmv, patient, Treatment) %>% 
#   unique() %>% 
#   select(-patient) %>% 
#   gtsummary::tbl_summary(by = cmv) %>% 
#   gtsummary::add_p() %>% 
#   modify_spanning_header(all_stat_cols() ~ "**Day 21 MiXCR data**")

x <- ggplot(BetaSize, aes(x = PC1, y = size_counts))+
  xlab("CMV score")+
  ylab("Post cycle 1 number of clones")+
  facet_wrap(~size, nrow = 1, scales = 'free')+
  geom_point(color = "#134B73", alpha = 0.6, size = 2)+
  geom_smooth(method = "lm", linetype = "dashed", color = "#EB7926", alpha = 0.4, fill ="#ABC9C8" )+
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
        strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16))

ggsave("/results/FigS2c.pdf", plot = x, width = 8, height = 8)

#rho = 0.42 with p = 4.4e-9 for > 0.5%
#rho = 0.51 with p = 2.4e-13 for < 0.5%

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
rownames(SubsetCounts) <- counts$Gene

#Metadata and counts data do not have matching samples. Match them and subset the matches 
#see if the counts column names are in the rownames of the metadata
#Find where the two dataframes intersect and filter according to the matches

all(colnames(SubsetCounts) %in% rownames(FilteredMeta))
Matching <- intersect(colnames(SubsetCounts), rownames(FilteredMeta))
FilteredMeta <- subset(FilteredMeta, (sample %in% Matching))
rownames(FilteredMeta) <- FilteredMeta$sample
SubsetCounts <- SubsetCounts %>% 
  dplyr::select(contains(Matching))
rownames(SubsetCounts) <- counts$Gene

#Not only do all the samples have to match in the metadata and the counts, but they have to be in the same order. Reorder the counts data 
#Select metadata columns of interest

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
#FIGURE S2d- Cytotoxicity score pre- and post-treatment according to treatment type and CMV status

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

x <- ggplot(subset(PlotData, combination == "sICB"), aes(x = cycle, y = Cytotoxicity))+
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

  ggsave("/results/FigS2d.pdf", plot = x, width = 8, height = 8)

################################################################################
##EXTENDED FIGURE III- IRAES##
################################################################################

#creating time-dependent covariates for toxicities
# temp df with just patient ID and OS time/censoring

tox2 <- tox %>% 
  mutate(Grade3 = ifelse(grade > 2, "Yes", NA)) %>% 
  mutate(Grade1_2 = ifelse(grade == 1 | grade == 2, "Yes", NA)) %>% 
  pivot_wider(names_from = organ, values_from = AI) %>% 
  mutate(across(18:45, ~ ifelse(. == "1", "Yes", NA))) %>% 
  mutate(dermatitis = skin) %>% 
  mutate(patient = as.character(ID)) %>% 
  dplyr::select(-notes)

temp <- surv_last_update %>% 
  filter(Rx != "Ipilimumab" & Rx != "RelatNivo") %>% 
  mutate(patient = as.character(ID)) %>% 
  dplyr::select(patient, months_death, death_status)

surv_last_update <- surv_last_update %>% 
  mutate(patient = as.character(ID))

TimeToTox <- left_join(tox2, surv_last_update)

TimeToTox <- TimeToTox %>% 
  filter(Rx != "Ipilimumab" & Rx != "RelatNivo") %>% 
  filter(cancer == "Melanoma" & !grepl("Adj", Rx))

#OS sets range for tstart and tstop

OS <- tmerge(temp, temp, id = patient, prog = event(months_death, death_status))

# TimeCovariateIrAEos merges OS with surv_last_update  and creates new long format for each patient
# with one line for each measured toxicity and the associated time to toxicity
# from either baseline or the previous toxicity.
# tdc() highlights which variables are time-dependent and states their time measurement (time_to_tox).
# surv_last_update is set so that once toxicity has occurred it remains so at later timepoints.
# TimeCovariateIrAEos is then merged with df to allow inclusion of other time-independent covariates in cox model

TimeCovariateIrAEos <- tmerge(OS, TimeToTox, id = patient, thyroiditis = tdc(time_to_tox, thyroiditis), pneumonitis = tdc(time_to_tox, pneumonitis),
                              neuro = tdc(time_to_tox, neuro), nephritis = tdc(time_to_tox, nephritis), myalgia = tdc(time_to_tox, myalgia),
                              hypophysitis = tdc(time_to_tox, hypophysitis), hepatitis = tdc(time_to_tox, hepatitis), gastritis = tdc(time_to_tox, gastritis),
                              dermatitis = tdc(time_to_tox, dermatitis), colitis = tdc(time_to_tox, colitis), arthritis = tdc(time_to_tox, arthritis),
                              grade_2 = tdc(time_to_tox, Grade1_2), grade_3 = tdc(time_to_tox, Grade3)) %>%
  mutate_all(~replace_na(., "No")) %>% 
  merge(surv_last_update, by = c("patient"), all = T)

#Final cox regression using time-independent and time-dependent covariates.

SpecificToxOS <- coxph(Surv(tstart, tstop, prog) ~ age + sex + BRAF_status + combination + thyroiditis + pneumonitis + neuro + nephritis + myalgia + hypophysitis + hepatitis + gastritis + dermatitis + colitis + arthritis, TimeCovariateIrAEos)

AllToxOS <- coxph(Surv(tstart, tstop, prog) ~ age + sex + BRAF_status + combination + grade_2, TimeCovariateIrAEos)

Grade3ToxOS <- coxph(Surv(tstart, tstop, prog) ~ age + sex + BRAF_status + combination + grade_3, TimeCovariateIrAEos)

#Extract results for specific toxicities 

SpecificToxsummary <- summary(SpecificToxOS)
coefSpecificTox <- SpecificToxsummary$coef[, 1]
coefSpecificTox <- coefSpecificTox[5:15]
ciSpecificTox <- as.data.frame(SpecificToxsummary$conf.int)
ciSpecificTox <- ciSpecificTox[5:15,]

#Create a data frame for plotting

PlotSpecificTox <- data.frame(
  covariate = c("thyroiditis", "pneumonitis", "neuro", "nephritis", "myalgia", "hypophysitis", "hepatitis", "gastritis", "dermatitis", "colitis", "arthritis"),
  HR = ciSpecificTox[, "exp(coef)"],
  lower = ciSpecificTox[,"lower .95"],
  upper = ciSpecificTox[,"upper .95"],
  p_value = SpecificToxsummary$coefficients[5:15, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)

PlotSpecificTox <- PlotSpecificTox %>% 
  mutate(Significant = ifelse(p_value < 0.05, "Significant", "Not significant")) %>% 
  mutate(Alpha = ifelse(Significant == "Significant", 1, 0.5))

#Extract results for any grade 3 irAE 

Grade3Toxsummary <- summary(Grade3ToxOS)
coefGrade3Tox <- Grade3Toxsummary$coef[, 1]
coefGrade3Tox <- coefGrade3Tox[5]
ciGrade3Tox <- as.data.frame(Grade3Toxsummary$conf.int)
ciGrade3Tox <- ciGrade3Tox[5,]

#Create a data frame for plotting

PlotGrade3Tox <- data.frame(
  covariate = c("Grade 3+ irAEs"),
  HR = ciGrade3Tox[, "exp(coef)"],
  lower = ciGrade3Tox[,"lower .95"],
  upper = ciGrade3Tox[,"upper .95"],
  p_value = Grade3Toxsummary$coefficients[5, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)

PlotGrade3Tox <- PlotGrade3Tox %>% 
  mutate(Significant = ifelse(p_value < 0.05, "Significant", "Not significant")) %>% 
  mutate(Alpha = ifelse(Significant == "Significant", 1, 0.5))

#Extract results for irAE regardless of grade 

AllToxsummary <- summary(AllToxOS)
coefAllTox <- AllToxsummary$coef[, 1]
coefAllTox <- coefAllTox[5]
ciAllTox <- as.data.frame(AllToxsummary$conf.int)
ciAllTox <- ciAllTox[5,]

#Create a data frame for plotting

PlotAllTox <- data.frame(
  covariate = c("Grade 1-2 irAEs"),
  HR = ciAllTox[, "exp(coef)"],
  lower = ciAllTox[,"lower .95"],
  upper = ciAllTox[,"upper .95"],
  p_value = AllToxsummary$coefficients[5, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)

PlotAllTox <- PlotAllTox %>% 
  mutate(Significant = ifelse(p_value < 0.05, "Significant", "Not significant")) %>% 
  mutate(Alpha = ifelse(Significant == "Significant", 1, 0.5))

Plot <- rbind(PlotSpecificTox, PlotGrade3Tox, PlotAllTox)
Plot <- Plot %>% 
  mutate(Grade = ifelse(covariate == "Grade 1-2 irAEs" | covariate == "Grade 3+ irAEs", "Grade", "Specific"))

Plot <- Plot %>%
  mutate(covariate = fct_relevel(covariate, "Grade 3+ irAEs", "Grade 1-2 irAEs"))


################################################################################
#FIGURE S3a- Survival impact of different irAEs

x <- ggplot(Plot, aes(x = HR, y = covariate)) +
  theme_classic()+
  theme(axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(hjust=1, size=12),
        axis.text.y = element_text(size = ifelse(Plot$covariate %in% c("Grade 1-2 irAEs", "Grade 3+ irAEs"), 14, 12)),
        axis.line = element_line(colour = "grey30", size = 0.2),
        axis.ticks = element_line(color= "grey30", size=0.2), 
        strip.background =element_rect(fill= "grey95"),
        strip.text = element_text(colour = "grey30"),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16))+
  geom_errorbar(data = subset(Plot, Significant == "Not significant" & Grade == "Specific"),
                aes(xmax = upper, xmin = lower), 
                width = 0, linetype = "longdash", colour ="#2171B5", alpha = 0.5) +
  geom_point(data = subset(Plot, Significant == "Not significant" & Grade == "Specific"), size = 5, alpha = 0.5, color = "#2171B5") +
  geom_errorbar(data = subset(Plot, Significant == "Significant" & Grade == "Specific"),
                aes(xmax = upper, xmin = lower), 
                width = 0, linetype = "longdash", colour = "#FFBB44", alpha = 0.5, size = 0.5) +
  geom_point(data = subset(Plot, Significant == "Significant" & Grade == "Specific"), size = 5, alpha = 0.5, color = "#FFBB44") +
  geom_errorbar(data = subset(Plot, Significant == "Not significant" & Grade == "Grade"),
                aes(xmax = upper, xmin = lower), 
                width = 0, linetype = "longdash", colour = "#2171B5", alpha = 1) +
  geom_point(data = subset(Plot, Significant == "Not significant" & Grade == "Grade"), size = 7, alpha = 1, color = "#2171B5") +
  geom_errorbar(data = subset(Plot, Significant == "Significant" & Grade == "Grade"),
                aes(xmax = upper, xmin = lower), 
                width = 0, linetype = "longdash", colour = "#FFBB44", alpha = 1, size = 0.5) +
  geom_point(data = subset(Plot, Significant == "Significant" & Grade == "Grade"), size = 7, alpha = 1, color = "#FFBB44") +
  scale_x_log10()+
  geom_vline(xintercept = 1,linetype ="dotted", color = "grey30") +
  geom_hline(yintercept = 11.75, color = "grey30",linetype = "dashed")+
  xlab("Adjusted hazard ratios")+
  ylab("Covariate")+
  geom_text(aes(label = paste0("p = ", signif(p_value, digits = 2))),
            vjust = -1, hjust = 0.5, size = 5, color = "grey30")

ggsave("/results/FigS3a.pdf", plot = x, width = 8, height = 8)

#Grade 1-2 irAEs: p = 0.00022, HR = 0.50
#Arthritis: p = 0.039, HR= 0.46
#Dermatitis: p = 0.0088, HR = 0.60


################################################################################

#Mutate ID to match survival data

TimeToTox <- tox %>% 
  mutate(patient = ID) %>% 
  mutate(grade = ifelse(is.na(grade), 0, grade)) %>% 
  dplyr::select(-ID) %>% 
  mutate(Grade_3 = ifelse(grade > 2, 1, 0)) %>% 
  mutate(Grade1_2 = ifelse(grade == 1 | grade == 2, 1, 0)) %>% 
  mutate(time_to_tox = as.numeric(time_to_tox))

#Get survival data for grade 3/4 Melanoma

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
  mutate(Adjuvant = ifelse(Rx == "AdjNivo" | Rx == 'AdjPembro', 0, 1))

#347 melanoma patients 
#308 metastatic melanoma and 39 stage II/III melanoma
#313 patients have CMV serology
#163 CMV- and 150 CMV+

#Modify time for grade 1-2 toxicities  

Grade12 <- SurvToxMel %>% 
  dplyr::select(Grade1_2, patient, time_to_tox, cmv, age, sex, Rx, BRAF_status, grade, months_death, death_status, Adjuvant) %>% 
  mutate(time_to_tox = ifelse(Grade1_2 == 0, NA, time_to_tox))

#Identify first grade 1-2 toxicity 

Grade12 <- Grade12 %>% 
  group_by(patient) %>% 
  slice_min(order_by = time_to_tox, with_ties = F) %>% 
  ungroup()

#Create competing risks column- censor at 24 months as there are no events beyond this point 

Grade12 <- Grade12 %>% 
  mutate(time_to_tox = ifelse(is.na(time_to_tox), months_death, time_to_tox)) %>% 
  mutate(time_to_tox = ifelse(time_to_tox > 24, 24, time_to_tox)) %>% 
  mutate(Status = as.factor(ifelse(Grade1_2 == 1, 1, 
                                   ifelse(Grade1_2 == 0 & death_status == 0, 0, 
                                          ifelse(death_status == 1 & Grade1_2 == 0, 2, NA))))) %>% 
  filter(Rx != "RelatNivo" & Rx != "Ipilimumab") %>% 
  mutate(combination = ifelse(Rx == "IpiNivo", "cICB", "sICB"))

#342 patients excluding Ipi and RelatNivo
#308 with CMV serology 

################################################################################
#FIGURE 3Sb- Multivariate all grade toxicities HR

CumincCMV <- cuminc(Surv(time_to_tox, Status) ~ cmv , data = Grade12) 

#Extract p-value

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
  labs(x = "Months", y = "Cumulative Incidence of Grade 1-2 irAEs") +
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

ggsave("/results/FigS3b.pdf", plot = x, width = 8, height = 8)

#p = 0.99
  
################################################################################

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
#FIGURE S3c- Time to grade 3 toxicity according to CMV 

# #Analysis restricted to metastatic melanoma patients receiving IpiNivo

MMgrade3 <- Grade3 %>%
  filter(Rx == "IpiNivo")

CumincCMV <- cuminc(Surv(time_to_tox, Status) ~ cmv , data = MMgrade3)

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
  labs(x = "Months", y = "Cumulative Incidence of Grade 3+ irAEs \n in cICB treated") +
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
  
  ggsave("/results/FigS3c.pdf", plot = x, width = 8, height = 8)

#p = 0.0061

################################################################################
##EXTENDED FIGURE IV- TRANSCRIPTION FACTOR ANALYSIS##
################################################################################

#Filter metastatic melanoma data for C1 expression and create a new sample variable 

metadata <- mmdata %>% 
  filter(Rx != "Ipilimumab" & patient != 25 & sample != "360_C2") %>% 
  filter(cycle == "C1") %>% 
  mutate(sample = paste0(sample, "_CD8_None")) %>% 
  mutate(sample = as.factor(sample)) 

#For the batch_data, select CD8+ T cell samples taken at baseline (format C1_CD8_None)

batch <- batch_data %>%
  mutate(patient = as.numeric(patient)) %>% 
  mutate(sample = as.factor(sample)) %>% 
  filter(grepl(paste("C1_CD8_None", sep = "|"), sample))

#Merge batch data with metadata which can be used in DESeq2 for batch correction 
#merge using patient, cycle, and sample and discard anything that does not match (unmatched = "drop")

MetadataBatch <- left_join(metadata, batch, by = c("patient", "cycle", "sample") , unmatched = "drop", keep = F)

NewSamples <- c(262,264,271, 274, 275, 276, 279, 281, 286, 291, 297, 302,265, 269, 273, 284, 288, 292)
MissingBatches2022.Sep <- c(233, 239, 243, 254, 261)
MissingBatchesMay2021 <- c(224)
MissingBatchesDec2020 <- c(161)
MissingBeforeSep2019 <- c(2, 26, 41, 44, 45, 81, 86)
MissingSep2019 <- c(153, 155)

FilteredMeta <- MetadataBatch %>%
  mutate_if(is.character, as.factor) %>% 
  mutate(cmv = ifelse(cmv == "CMV+", "CMVpositive",
                      ifelse(cmv == "CMV-", "CMVnegative", NA))) %>%
  mutate(origBatch = case_when(
    patient %in% NewSamples ~ "Jan2024",
    TRUE ~ origBatch)) %>% 
  mutate(origBatch = case_when(patient %in% MissingBatches2022.Sep ~ "2022.Sep", TRUE ~ origBatch)) %>% 
  mutate(origBatch = case_when(patient %in% MissingBatchesMay2021 ~ "May_2021", TRUE ~ origBatch)) %>% 
  mutate(origBatch = case_when(patient %in% MissingBatchesDec2020 ~ "Dec_2020", TRUE ~ origBatch)) %>% 
  mutate(origBatch = case_when(patient %in% MissingBeforeSep2019 ~ "Before_SEP2019", TRUE ~ origBatch)) %>% 
  mutate(origBatch = case_when(patient %in% MissingSep2019 ~ "Sep_19", TRUE ~ origBatch)) %>% 
  mutate(origBatch = ifelse(patient > 302, "June_2024", origBatch)) %>% 
  filter(!is.na(cmv))

#Ensure the rownames still correspond to the sample ids

rownames(FilteredMeta) <- FilteredMeta$sample

#Subset counts data for baseline CD8+ T cell samples (C1_CD8_None)
#Add Ensembl ids as rownames in this new dataset

rownames(counts) <- counts$Gene
SubsetCounts <- counts %>%
  dplyr::select(matches(c("C1_CD8_None")))
rownames(SubsetCounts) <- counts$Gene

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
#Select metadata columns of interest: sample, sex, batch, age as continuous variable and cmv status- remove columns that are not required 

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

#206 C1 samples with CMV serology

#CollecTRI is a comprehensie resource containing a curated collection of TFs and their transcriptional targets compiled from 12 different resources 
#Retrieve the human version from OmniPath 

Network <- readRDS("/data/TFnetwork.rds")

#Filter for transcription factors expressed in CD8+ T cells 

CountsFinalC12 <- CountsFinalC1
CountsFinalC12$symbol <- mapIds(org.Hs.eg.db, keys = rownames(CountsFinalC12), keytype = "ENSEMBL", column = "SYMBOL")
Network <- Network %>% 
  filter(source %in% CountsFinalC12$symbol)

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
  filter(padj < 0.05)
DEGs <- DEGs %>%
  distinct(symbol, .keep_all = TRUE)
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

#Filter top TFs in both signs

pvals <- TF$p_value
padj <- p.adjust(pvals, method = "fdr")

TF$padj <- padj

TF <- TF %>% 
  filter(padj < 0.05) %>% 
  mutate(Score = score)

#Identify TFs in DEGs according to CMV at baseline 

ResultsOrdered <- ResultsOrdered %>% 
  mutate(TF = ifelse(symbol %in% TF$source, "TF", "Other")) %>% 
  mutate(TFdirection = ifelse(TF == "TF" & padj < 0.05 & log2FoldChange > 0, "TFhigh",
                              ifelse(TF == "TF" & padj < 0.05 & log2FoldChange < 0, "TFlow", 
                                     ifelse(TF == "TF" & padj > 0.05, "TFNone", 'Other')))) 

Labels <- ResultsOrdered %>% 
  filter(TF == "TF")

################################################################################
#FIGURE S4a- TFs on CMV volcano plot 

ylab <- expression(-log[10]~ (adjP))
xlab <- expression(log[2]~FC )
x <- ggplot(ResultsOrdered,
       aes(x = log2FoldChange,
           y = -log10(padj),
           col = TFdirection,
           label = symbol,
           alpha = ifelse(TF == "TF", 3, 0.1)))+
  geom_point(aes(color = TFdirection), size = 2)+
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-0.5, 17))+
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", alpha = 0.5)+
  geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.5)+
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
  theme(legend.position = "None")+
  geom_text_repel(data = Labels, colour = "grey30", size = 3.75, max.overlaps = 10, alpha = 0.6)+
  scale_color_manual(values=c("grey95", "#FFBB44", "#2171B5", "#ABC9C8"))

  ggsave("/results/FigS4a.pdf", plot = x, width = 8, height = 8)

################################################################################

#Identify Transcription factor correlations with the CMV score- focus on TFs that are induced by CMV

metadata <- mmdata %>% 
  filter(Rx != "Ipilimumab" & patient != 25 & sample != "260_C2") %>% 
  filter(cycle == "C1" | cycle == "C2") %>% 
  mutate(sample = paste0(sample, "_CD8_None")) %>% 
  mutate(sample = as.factor(sample)) 

#For the batch_data, select cells of interest in the format 'CX_CDY' referring to cycle and cell marker (eg CD8, CD19, etc...)

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

#Subset counts data for cells of interest in the format '_CX_CD8_Stim'
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
#Select metadata columns of interest: sample, sex, cycle, batch, age as continuous variable, cmv sereostatus, etc...

SubsetReordered <- SubsetCounts[, match(rownames(FilteredMeta), colnames(SubsetCounts))]
all(colnames(SubsetReordered) == rownames(FilteredMeta))
rownames(SubsetReordered) <- rownames(counts) 
SubsetReordered$id <- rownames(SubsetReordered)
FinalMeta <- FilteredMeta %>% 
  dplyr::select(sample, sex, cycle, age, cmv, Rx, patient)

#Filter out genes which have low reads as these could easily be classed as differentially expressed. 10 is the count threshold suggested by the vignette.
#Do a final check that the samples in the counts data are matching and in the same order as the samples in the metadata 

keep <- rowMeans(SubsetReordered[,1:(ncol(SubsetReordered)-1)]) >= 10
SubsetReordered <- SubsetReordered[keep, ] 
rownames(SubsetReordered) <- SubsetReordered$id
CountsFinal <- dplyr::select(SubsetReordered, -id)
rownames(CountsFinal) <- rownames(SubsetReordered)
rownames(FinalMeta) <- FinalMeta$sample
all(colnames(CountsFinal) == rownames(FinalMeta))

#417 samples

#Filter for genes that reach threshold of significance to be included in the score 

CMVgenes <- ResultsOrdered %>% 
  filter(padj < 1.581139e-14)
NormCounts <- baselineNormaliseTransform(files = CountsFinal, samples = FinalMeta)
NormCounts$symbol <- mapIds(org.Hs.eg.db, keys = rownames(NormCounts), keytype = "ENSEMBL", column = "SYMBOL")
CMVinduced <- NormCounts %>% 
  dplyr::select(-X) %>% 
  filter(rownames(NormCounts) %in% rownames(CMVgenes)) %>% 
  dplyr::select(-c(symbol))
CMVinducedWide <- t(CMVinduced)
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

#Filter for transcription factors whose activity is induced in CMV seropositivity 

TFsInduced <- TF %>% 
  filter(score > 0) %>% 
  dplyr::select(source)
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

#Join TF data with CMV score 

CountsMeta <- left_join(CountsMeta, PCAmeta)
CountsMetaLong <- CountsMeta %>% 
  pivot_longer(names_to = "TF", values_to = "NormExp", c(8:16))

################################################################################
#FIGURE S4b- Active TFs correlations with CMV score at day 21

x <- ggplot(subset(CountsMetaLong, cycle == "C2"), aes(x = PC1, y = NormExp))+
  xlab("CMV score")+
  ylab("Normalised TF expression")+
  geom_point(color = "#ABC9C8", alpha = 0.6)+
  geom_smooth(method = "lm", linetype = "dashed", color = "#EB7926", alpha = 0.4, fill = "grey60")+
  stat_cor(method = "spearman", size = 5)+
  facet_wrap(~TF, scales = "free")+
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
        strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16))

ggsave("/results/FigS4b.pdf", plot = x, width = 8, height = 8)

#TBX21 is the most correlated TF 

################################################################################

#Explore T-bet dynamics in CMV and treatment type

metadata <- mmdata %>% 
  filter(patient != 25 & Rx != "Ipilimumab" & sample != "260_C2") %>% 
  filter(cycle == "C1"|cycle == "C2") %>% 
  mutate(sample = paste0(sample, "_CD8_None")) %>% 
  mutate(sample = as.factor(sample)) %>% 
  mutate(combination = ifelse(Rx == "IpiNivo", "cICB", "sICB")) 

#For the batch_data, select CD8 T cell samples taken at baseline (format C1_CD8_)

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

rownames(counts) <- counts$Gene
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
rownames(SubsetReordered) = rownames(counts) 
SubsetReordered$id <- rownames(SubsetReordered)
FinalMeta <- FilteredMeta %>% 
  mutate(`Responder status` = ifelse(months_progression > 60, "Durable Responder",
                                     ifelse(progression == 1 & months_progression < 60, "Progressor", NA))) %>% 
  dplyr::select(sample, patient, sex, age, cmv, Rx, origBatch, combination, `Responder status`, cycle, progression, months_progression) %>% 
  mutate(sex = as.factor(sex)) %>% 
  mutate(age = scale(age, center = T))
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
  filter(symbol %in% c("TBX21")) %>% 
  as.data.frame()
rownames(Subset) <- Subset$symbol
Subset <- Subset %>% 
  dplyr::select(-c(symbol, X))
SubsetWide <- t(Subset)
rownames(SubsetWide) <- rownames(FinalMeta)
SubsetWide <- as.data.frame(SubsetWide)
SubsetWide$sample <- rownames(SubsetWide)
FinalMeta <- left_join(FinalMeta, SubsetWide, by = "sample")

#Correlating T-bet with clone size 

TbetClones <- FinalMeta %>% 
  filter(grepl("_C2_", sample)) %>% 
  dplyr::select(TBX21, patient, sample) %>% 
  mutate(patient = as.numeric(patient))

#Determine the number of clones by size according to CMV 

#Determine total read counts for each individual 

cloneCounts <- tcr_processed %>%
  filter(grepl(paste("TRB", sep = "|"), bestVHit)) %>% 
  group_by(ID) %>%
  filter(sample != "250_C2") %>% 
  summarise(total_cell_count = sum(cloneCount))

#Filter for cycle 2 beta chains

BetaSizeC2 <- tcr_processed %>% 
  filter(grepl(paste("TRB", sep = "|"), bestVHit)) %>% 
  filter(sample != "250_C2") %>% 
  filter(grepl("_C2_", ID))

#Join total counts with filtered counts 

BetaSizeC2 <- left_join(BetaSizeC2, cloneCounts)

#Determine the number and proportion of clones by size per individual at day 21

BetaSizeC2 <- BetaSizeC2 %>% 
  mutate(repertoire_occupancy = cloneCount/total_cell_count * 100) %>% 
  mutate(size = ifelse(repertoire_occupancy > 0.5, "large", 
                       ifelse(repertoire_occupancy <= 0.5, "intermediate",  NA))) %>% 
  group_by(sample) %>% 
  summarise(`> 0.5%` = sum(grepl("large", size)),
            `< 0.5%` = sum(grepl("intermediate", size))) %>% 
  unique()

#join with metadata 

BetaSizeC2 <- left_join(BetaSizeC2, mmdata)

BetaSizeC2long <- BetaSizeC2 %>% 
  filter(!grepl("Adj", Rx)) %>% 
  filter(cancer == "Melanoma") %>% 
  mutate(combination = ifelse(Rx == "IpiNivo", "cICB", "sICB")) %>% 
  pivot_longer(names_to = "size", cols = c(`> 0.5%`, `< 0.5%`), values_to = "size_counts") %>% 
  dplyr::select(size, size_counts, cancer, Rx, combination, patient, sample) %>% 
  mutate(sample = paste0(sample, "_CD8_None"))

#Relevel factor 

BetaSizeC2long$size <- factor(BetaSizeC2long$size, levels = c("< 0.5%", "> 0.5%"))

#Join with TBX21 data 

TbetClones <- left_join(BetaSizeC2long, TbetClones)

################################################################################
#FIGURE S4c- TBX21 correlation with clone size 

x <- ggplot(TbetClones, aes(x = TBX21, y = size_counts))+
  xlab("Normalised TBX21 expression")+
  ylab("Number of clones")+
  geom_point(color = "#134B73", alpha = 0.6, size = 2)+
  geom_smooth(method = "lm", linetype = "dashed", color = "#EB7926", alpha = 0.4, fill ="#ABC9C8" )+
  stat_cor(method = "spearman", size = 5)+
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
        strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 18))+
  facet_wrap(~size, scales = "free", nrow = 1)

#Large clones: rho = 0.5, p = 5.1e-13
#Small clones: rho = 0.43, p = 2.7e-9

ggsave("/results/FigS4c.pdf", plot = x, width = 8, height = 8)

################################################################################

#CMV-reactive clones and TBX21 

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

#Plot a UMAP of cell subsets in the single cell data- Nature Medicine Review 

SeuratMeta <- read.csv("/ceph/project/fairfaxlab/gmilotay/scRNAseq/data/scvi_umap_embeddings.csv")

#Subset the barcodes and UMAP coordinates 

SeuratMeta <- SeuratMeta %>% 
  select(barcode_pool, X_umap_1, X_umap_2)

cd8@meta.data <- left_join(cd8@meta.data, SeuratMeta)

#Subset sICB patients used in analysis 

rownames(cd8@meta.data) <- cd8$barcode_pool

sICB <- subset(cd8, subset = Rx != "IpiNivo")

sICB <- subset(sICB, subset = cycle %in% c("C1", "C2"))

#Plot UMAP 

UMAP <- as.matrix(sICB@meta.data[, c("X_umap_1", "X_umap_2")])

#Add UMAP coordinates to the embeddings slot

sICB[["umap"]] <- CreateDimReducObject(embeddings = UMAP, key = "UMAP_", assay = DefaultAssay(sICB))

################################################################################
#FIGURE S4d- Single Cell UMAP

DimPlot(sICB,
        label = T,
        group.by = "Subset", 
        label.color = "black",
        label.box = T,
        label.size = 3, 
        repel = T)+
  ggtitle("")+
  theme(legend.position = "None")+
  coord_cartesian(xlim = c(-7, 7), ylim = c(-10, 7))

################################################################################

#Plot a map of differential gene expression of each cluster based on marker genes 

Idents(sICB) <- "Subset"

#Normalise and Scale the data 

sICB <- NormalizeData(sICB)
sICB <- ScaleData(sICB)

#Perform differential expression analysis 

DEGs <- FindAllMarkers(sICB)

#Identify the top 5 markers for each cluster

TopMarkers <- DEGs %>%
  group_by(cluster) %>%
  top_n(5, wt = avg_log2FC) %>%  
  arrange(cluster, desc(avg_log2FC)) %>% 
  select(gene) 

#Extract scaled expression of marker genes 

TopMarkers <- TopMarkers$gene
TopMarkers <- unique(TopMarkers)
TopMarkers <- TopMarkers[which(TopMarkers %in% rownames(sICB@assays$RNA@scale.data))]
Data <- GetAssayData(object = sICB[["RNA"]], slot = "scale.data")[TopMarkers, ]
Data <- as.data.frame(t(as.matrix(Data)))
Data$ID <- rownames(Data)

sICB$ID <- rownames(sICB@meta.data)

Subsets <- sICB@meta.data %>% 
  select(ID, Subset)

Meta <- left_join(Subsets, Data)

Meta <- Meta %>%
  pivot_longer(names_to = "Gene", values_to = "Expression", cols = c(3:39)) %>% 
  group_by(Subset, Gene) %>% 
  summarise(Expression = mean(Expression)) %>% 
  pivot_wider(names_from = Gene, values_from = Expression)

Plot <- as.matrix(Meta[, -1])
rownames(Plot) <- Meta$Subset

Plot <- t(Plot)
col1 <- rgb(228, 28, 30, max = 255)  
col2 <- rgb(16, 98, 188, max = 255) 

min<- round(min(Plot, na.rm = TRUE), 1)
max <- round(max(Plot, na.rm = TRUE), 1)

col_fun <- colorRamp2(c(min, 0, max), c(col2, "white", col1))

################################################################################
#FIGURE S4e- Differential expression of marker genes in single cell data 

Heatmap(Plot, name = "Scaled\nExpression", 
        border = TRUE, 
        col = col_fun, 
        column_names_side = "top",
        column_names_rot = 45, 
        column_names_centered = TRUE,
        row_names_side = "left",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        column_names_gp = gpar(fontsize = 11),
        row_names_gp = gpar(fontsize = 9),
        rect_gp = gpar(col = "white", lwd = 0.1),
        heatmap_legend_param = list(title = "Scaled\nExpression", at = c(min, 0, max), 
                                    legend_height = unit(2, "cm"))
) 

################################################################################
#FIGURE S4f-Feature plot for TBX21 

FeaturePlot(sICB, features = "TBX21")+
  ggtitle("")+
  coord_cartesian(xlim = c(-7, 7), ylim = c(-10, 7))

################################################################################

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

################################################################################
#Figure S4g- CMV-reactive clones expression of TBX21 

x <- ggplot(subset(CMVreactive2, cmv_trb == "Putative CMV clone"), aes(x = cycle, y = tbet))+
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

ggsave("/results/FigS4d.pdf", plot = x, width = 8, height = 8)

################################################################################
##EXTENDED FIGURE V- TBX21 SURVIVAL##
################################################################################

FinalMetaC2 <- FinalMeta %>% 
  filter(cycle == "C2") %>% 
  mutate(patient = as.character(patient))

#Survival according to TBX21 expression at C2 

survivalTBX21 <- FinalMetaC2 %>% 
  dplyr::select(TBX21, patient)

survivalTBX21 <- left_join(surv_last_update, survivalTBX21)

################################################################################
#FIGURE S5a- OS Kaplan Meier of TBX21 expression at day 21

MedianTBX21median <- median(subset(survivalTBX21, !is.na(TBX21))$TBX21)
survivalTBX21 <- survivalTBX21 %>% 
  mutate(DichotTBX21 = ifelse(TBX21 >= MedianTBX21median, "> Median TBX21", "< Median TBX21"))

fit <- survfit(Surv(months_death, death_status) ~ DichotTBX21, data = survivalTBX21)
summary(fit)$table
log_rank_test <- survdiff(Surv(months_death, death_status) ~ DichotTBX21, data = survivalTBX21)
p_value <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)
cox <- coxph(Surv(months_death, death_status) ~ DichotTBX21, data = survivalTBX21)
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
                 xlab = "Months", ylab = "Probability of overall survival", legend.title = "TBX21 expression", legend.labs = c("< Median", "> Median"),
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
  geom_segment(x = time_point_5yr, xend = time_point_5yr, y = -Inf, yend = 0.7, color = "grey60", lty="dashed")+
  annotate("text", x=time_point_1yr+4.5, y=surv_prob_1yr[1] -0.15, label=paste0(round(surv_prob_1yr[[1]],3)*100,"%"), size=5, colour="#2171B5")+
  annotate("text", x=time_point_1yr+4, y=surv_prob_1yr[2]+0.05, label=paste0(round(surv_prob_1yr[[2]],3)*100,"%"), size=5, colour="#FFBB44")+
  annotate("text", x=time_point_1yr, y=surv_prob_1yr[1]+0.2, label="1 year", size=5, colour="grey60")+
  annotate("text", x=time_point_5yr+4, y=surv_prob_5yr[1]-0.1, label=paste0(round(surv_prob_5yr[[1]],3)*100,"%"), size=5, colour="#2171B5")+
  annotate("text", x=time_point_5yr+4.5, y=surv_prob_5yr[2]+0.05, label=paste0(round(surv_prob_5yr[[2]],3)*100,"%"), size=5, colour="#FFBB44")+
  annotate("text", x=time_point_5yr, y=surv_prob_5yr[1]+0.38, label="5 years", size=5, colour="grey60")+
  annotate("text", x=time_point_5yr-20, y = 0.85, label=paste0("HR for death: ",
                                                               round(PlotsICB$HR, digits = 2),"\n (95% CI: ", round(PlotsICB$lower, digits = 2),
                                                               "-",round(PlotsICB$upper, digits = 2), ") \n Log-rank \n", round(p_value, digits = 3), sep = ""), size=5, fontface="plain")
km

ggsave("/results/FigS5a.pdf", plot = km$plot, width = 8, height = 8)

#p = 0.026
#Median 34.4 vs NA
#HR = 0.62

################################################################################

survivalTBX21 <- survivalTBX21 %>% 
  filter(!is.na(TBX21)) %>% 
  mutate(Sex = sex) %>% 
  mutate(Age = age) %>% 
  mutate(Treatment = Rx) %>% 
  mutate(`BRAF status` = BRAF_status)

TBX21cox <- coxph(Surv(months_progression, progression) ~ Sex + Age + combination + `BRAF status`+ cmv + `TBX21`,
                  data = survivalTBX21)

TBX21summary <- summary(TBX21cox)
TBX21coef <- TBX21summary$coef[, 1]
ciTBX21 <- as.data.frame(TBX21summary$conf.int)

#Create a data frame for plotting

PlotTBX21 <- data.frame(
  covariate = c("Female vs \n Male", "Age", "sICB vs \n cICB"," BRAF Wild-Type vs \n Mutant", "CMV+ vs \n CMV-", "TBX21 expression"),
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
#FIGURE S5b- PFS according to TBX21 expression at day 21 

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
  xlab("Adjusted hazard ratios day 21")+
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

ggsave("/results/FigS5b.pdf", plot = x, width = 8, height = 8)

################################################################################

##EXTENDED FIGURE VI- EPIDEMIALOGICAL STUDY##

#Determine if patients are metastatic, melanoma and their melanoma subtype 

SurvivalAllCancer <- clinicalData %>% 
  filter(!is.na(cmv)) %>% 
  dplyr::select(months_death, months_progression, progression, death_status, cmv, Rx, PtName, age, sex, BRAF_status, subtype, cancer, age_mel, UnknownPrimary) %>% 
  mutate(Cancer = ifelse(cancer == "Melanoma", "Melanoma", "Non-Melanoma")) %>% 
  mutate(intent = ifelse(grepl("Adj", Rx), "adjuvant", "palliative")) %>% 
  mutate(subtype = ifelse(Cancer == "Melanoma" & subtype == "E" | Cancer == "Melanoma" & subtype == "NOS" | Cancer == "Melanoma" & subtype == "small cell", "Cutaneous",
                          ifelse(Cancer == "Melanoma" & subtype == "M", "Mucosal",
                                 ifelse(Cancer == "Melanoma" & subtype == "U", "Uveal", NA))))

#Prepare data for joining with Biobank data 

AllCancer <- SurvivalAllCancer %>% 
  filter(subtype == "Cutaneous") %>% 
  filter(Rx != "Ipilimumab" & Rx != "RelatNivo") %>% 
  mutate(CMV = ifelse(cmv == "CMV+", 1,
                      ifelse(cmv == 'CMV-', 0, NA))) %>% 
  mutate(status = ifelse(Cancer == "Melanoma" & intent == "adjuvant", "Mel-Adj",
                                ifelse(Cancer == "Melanoma" & intent == "palliative" & BRAF_status == "Wild-type" & subtype == "Cutaneous", "MM-WT", 
                                       ifelse(Cancer == "Melanoma" & intent == "palliative" & BRAF_status == "Mutant" & subtype == "Cutaneous", "MM-Mut", NA)))) %>% 
  mutate(Status = ifelse(Cancer == "Melanoma" & intent == "palliative", "Metastatic Melanoma", "Other")) %>% 
  mutate(assessment_centre = "None") %>% 
  filter(!PtName %in% c("193", "197", "300", "408")) %>% #These patients were stage IV resected 
  dplyr::select(age, CMV, status, assessment_centre, Status) %>% 
  mutate(cancer = "Cancer") %>% 
  filter(!is.na(age))

#Prepare Biobank data 

Biobank <- allmat %>% 
  mutate(cmv = ifelse(cmv == 2, "CMV+", 
                      ifelse(cmv == 1, "CMV-", NA)))

Biobank <- Biobank %>% 
  dplyr::select(-sex) %>% 
  mutate(CMV = ifelse(cmv == "CMV+", 1, 0)) %>% 
  mutate(cancer = "False") %>% 
  mutate(status = "Healthy") %>% 
  mutate(Status = "Healthy") %>% 
  dplyr::select(-cmv) %>% 
  filter(!is.na(age))

#Bind healthy donors to cancer cohort 

HealthyCancerCombined <- rbind(AllCancer, Biobank)

#Remove anyone above the age of 70 or below the age of 40 

MetastaticHealthy <- HealthyCancerCombined %>% 
  filter(age >= 40 & age <= 70) %>% 
  filter(Status %in% c("Healthy", "Metastatic Melanoma")) 

assessmentCentres <- MetastaticHealthy %>% 
  filter(assessment_centre != "None")
assessmentCentres <- unique(assessmentCentres$assessment_centre)

#Initialize results dataframe

results <- data.frame(
  Center = character(0),
  Odds_Ratio = numeric(0),
  Lower_CI = numeric(0),
  Upper_CI = numeric(0),
  pval = numeric(0)
)

#Loop through assessment centers and perform Fisher's exact test

for (center in assessmentCentres) {
  center_data <- MetastaticHealthy %>% 
    filter(assessment_centre == center | assessment_centre == "None")
  
  #Fisher's test for CMV status in metastatic melanoma vs healthy donors
  
  fisher_test <- fisher.test(center_data$CMV, center_data$assessment_centre)
  odds_ratio <- fisher_test$estimate
  conf_int <- fisher_test$conf.int
  p_value <- fisher_test$p.value
  number <- center_data %>% 
    filter(assessment_centre == center) 
  number <- nrow(number)
  
  #Store the results
  results <- rbind(results, data.frame(
    Center = center,
    Odds_Ratio = odds_ratio,
    Lower_CI = conf_int[1],
    Upper_CI = conf_int[2],
    pval = p_value,
    Participants = number
  ))
}

################################################################################
#FIGURE S6a- Incidence of CMV in metastatic melanoma cohort compared to individual recruitment centres for the UK Biobank

results <- results %>% 
  mutate(Significant = ifelse(pval < 0.05, "Significant", "Not significant")) %>% 
  mutate(Alpha = ifelse(Significant == "Significant", 1, 0.5)) 
results <- results %>%
  mutate(Lower_CI_flag = ifelse(`Lower_CI` < 0.1, TRUE, FALSE)) %>% 
  mutate(`Lower_CI` = ifelse(Center == "11023", 0.1, `Lower_CI`))

x <- ggplot(results, aes(x = `Odds_Ratio`, y = Center)) +
  geom_errorbar(data = subset(results, Significant == "Significant"),
                aes(xmax = `Upper_CI`, xmin = `Lower_CI`), 
                width = 0, linetype = "longdash", colour = "#EF6548", alpha = 1.2, size = 0.5) +
  geom_point(data = subset(results, Significant == "Significant"), aes(size = Participants), alpha = 0.8, color = "#EF6548") +
  geom_segment(data = subset(results, Lower_CI_flag == TRUE),
               aes(x = 0.15, xend = 0.1, y = Center, yend = Center),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "#9E9AC8", size = 0.5, linetype = "blank") +
  geom_errorbar(data = subset(results, Significant == "Not significant"),
                aes(xmax = `Upper_CI`, xmin = `Lower_CI`), 
                width = 0, linetype = "longdash", colour = "#9E9AC8", alpha = 1.2) +
  geom_point(data = subset(results, Significant == "Not significant"), aes(size = Participants), alpha = 0.8, color = "#9E9AC8") +
  scale_x_log10()+
  geom_vline(xintercept = 1,linetype ="dotted", color = "grey") +
  xlab("Odds ratios of CMV seropositivity in metastatic melanoma \n patients compared to UKB donors")+
  ylab("Centre")+
  theme_classic()+
  theme(axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(hjust=1, size=14),
        axis.text.y = element_text(size=14),
        axis.line = element_line(colour = "grey30", size = 0.2),
        axis.ticks = element_line(color= "grey30", size=0.2), 
        strip.background =element_rect(fill= "grey95"),
        strip.text = element_text(colour = "grey30"),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        strip.text.x = element_text(size = 18), 
        strip.text.y = element_text(size = 18),
        legend.text = element_text(size =14),
        legend.title = element_text(size = 18))+
  geom_text(aes(label = paste0("p = ", signif(pval, digits = 2))),
            vjust = -0.6, hjust = 1.2, size = 4, color = "grey30")

ggsave("/results/FigS6a.pdf", plot = x, width = 8, height = 8)

################################################################################
#FIGURE S6b- Age of primary diagnosis according to BRAF status and CMV serostatus 

x <- ggplot(subset(SurvivalAllCancer, !is.na(cmv) & Cancer == "Melanoma" &
                !is.na(BRAF_status) & BRAF_status != "?" & !Rx %in% c("Ipilimumab", "RelatNivo") &
                intent == "palliative" & subtype == "Cutaneous" & is.na(UnknownPrimary)), aes(x = cmv, y = age_mel))+
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
  ylab(paste("age at primary diagnosis (years)"))+
  labs(fill = "CMV status", color = "CMV status")+
  #stat_compare_means(method = "wilcox", label.x.npc = "center", comparisons = list(c("CMV+", "CMV-")), size = 5)
  stat_compare_means(comparisons = list(c("CMV+", "CMV-")), label = "p.signif", size = 7)

ggsave("/results/FigS6b.pdf", plot = x, width = 8, height = 8)

#Wild-type: p = 0.0019
#Mutant: p = 0.25

model <- lm(age_mel ~ BRAF_status*cmv, data = subset(SurvivalAllCancer, !is.na(cmv) & Cancer == "Melanoma" &
                                                       !is.na(BRAF_status) & BRAF_status != "?" & !Rx %in% c("Ipilimumab", "RelatNivo") &
                                                       intent == "palliative" & subtype == "Cutaneous" & is.na(UnknownPrimary)))

summary(model)

################################################################################

##END---------------------------------------------------------------------------