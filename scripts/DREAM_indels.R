#New file to check INDELS metrics for DREAM datasets results



#Set your own directory
setwd('~/Descargas/INDELinsilicos')



#LIBRARIES

load.lib<-c("stringr", "viridis", "data.table","ggplot2","gridExtra","UpSetR",
            "caret","catspec","devtools","tidyr","plyr","dplyr",
            "rpart","rpart.plot","rattle","e1071","openssl","sqldf")

install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) {install.packages(lib,dependencies=TRUE)}
sapply(load.lib,library,character=TRUE)



metrics_calculator_indel <- function(insilico_vcf, truth_insilico_vcf){
  
  truth_tab <- truth_insilico_vcf
  truth_tab_formated <- read.table(truth_tab)
  colnames(truth_tab_formated) <- c("#CHROM",  "POS" ,    "ID"  ,    "REF" ,    "ALT"  ,   "QUAL" ,   "FILTER" , "INFO" ,   "FORMAT",  "SPIKEIN")
  truth_tab_formated$mut <- paste0(truth_tab_formated$"#CHROM","_",truth_tab_formated$POS,"_",truth_tab_formated$REF,"_",truth_tab_formated$ALT)
  
  insilico_tab <- insilico_vcf
  insilico_tab_formated <- read.table(insilico_tab, fill= TRUE, row.names=NULL)
  insilico_tab_formated <- insilico_tab_formated[,1:7] 
  colnames(insilico_tab_formated) <- c("#CHROM",  "POS" ,    "ID"  ,    "REF" ,    "ALT"  ,   "QUAL" ,   "FILTER")
  insilico_tab_formated <- insilico_tab_formated[insilico_tab_formated$FILTER == "PASS",]
  insilico_tab_formated$mut <- paste0(insilico_tab_formated$"#CHROM","_",insilico_tab_formated$POS,"_",insilico_tab_formated$REF,"_",insilico_tab_formated$ALT)
  
  options(scipen = 999)
  TP <- insilico_tab_formated[insilico_tab_formated$mut %in% truth_tab_formated$mut,]
  FP <- insilico_tab_formated[!insilico_tab_formated$mut %in% truth_tab_formated$mut,]
  FN <- truth_tab_formated[!truth_tab_formated$mut %in% insilico_tab_formated$mut,]
  Recall <- nrow(TP)/nrow(truth_tab_formated) 
  Precision <-nrow(TP)/(nrow(TP)+nrow(FP)) 
  F1_Score <- 2 * Precision * Recall / (Precision +Recall)
  
  metrics_X <- data.frame(round(Recall,3), round(Precision,3), round(F1_Score,3), nrow(TP), nrow(FP), nrow(FN))
  colnames(metrics_X) <- c("Recall","Precision","F1_Score","FP","TP","FN")
  metrics_X
  
  return(metrics_X)}



#Set your own directory
setwd('~/Downloads')


# truth_insilico_1_indel_vcf <-"./truth.indel.synthetic.challenge.set1.vcf"
# truth_insilico_2_indel_vcf <-"./truth.indel.synthetic.challenge.set2.vcf"
truth_insilico_3_indel_vcf <-"./truth.indel.synthetic.challenge.set3.vcf"


# insilico_1_indels_bsc <- "./insilico_1_indels_bsc.vcf"
# insilico_1_result_indels_bsc <- metrics_calculator_indel(insilico_1_indels_bsc,truth_insilico_1_indel_vcf)
# insilico_1_result_indels_bsc
# 
# insilico_2_indels_bsc <- "./insilico_2_indels_bsc.vcf"
# insilico_2_result_indels_bsc <- metrics_calculator_indel(insilico_2_indels_bsc,truth_insilico_2_indel_vcf)
# insilico_2_result_indels_bsc

insilico_3_indels_bsc <- "./insilico_3_indels_bsc.vcf"
insilico_3_result_indels_bsc <- metrics_calculator_indel(insilico_3_indels_bsc,truth_insilico_3_indel_vcf)
# insilico_3_result_indels_bsc




# insilico_1_indels_hmf <- "./insilico_1_indels_hmf.vcf"
# insilico_1_result_indels_hmf <- metrics_calculator_indel(insilico_1_indels_hmf,truth_insilico_1_indel_vcf)
# insilico_1_result_indels_hmf
# 
# insilico_2_indels_hmf <- "./insilico_2_indels_hmf.vcf"
# insilico_2_result_indels_hmf <- metrics_calculator_indel(insilico_2_indels_hmf,truth_insilico_2_indel_vcf)
# insilico_2_result_indels_hmf

insilico_3_indels_hmf <- "./insilico_3_indels_hmf.vcf"
insilico_3_result_indels_hmf <- metrics_calculator_indel(insilico_3_indels_hmf,truth_insilico_3_indel_vcf)
# insilico_3_result_indels_hmf



# insilico_1_indels_oicr <- "./insilico_1_indels_oicr.vcf"
# insilico_1_result_indels_oicr <- metrics_calculator_indel(insilico_1_indels_oicr,truth_insilico_1_indel_vcf)
# insilico_1_result_indels_oicr
# 
# insilico_2_indels_oicr <- "./insilico_2_indels_oicr.vcf"
# insilico_2_result_indels_oicr <- metrics_calculator_indel(insilico_2_indels_oicr,truth_insilico_2_indel_vcf)
# insilico_2_result_indels_oicr

insilico_3_indels_oicr <- "./insilico_3_indels_oicr.vcf"
insilico_3_result_indels_oicr <- metrics_calculator_indel(insilico_3_indels_oicr,truth_insilico_3_indel_vcf)
# insilico_3_result_indels_oicr



# insilico_1_indels_charite <- "./insilico_1_indels_charite.vcf"
# insilico_1_result_indels_charite <- metrics_calculator_indel(insilico_1_indels_charite,truth_insilico_1_indel_vcf)
# insilico_1_result_indels_charite
# 
# insilico_2_indels_charite <- "./insilico_2_indels_charite.vcf"
# insilico_2_result_indels_charite <- metrics_calculator_indel(insilico_2_indels_charite,truth_insilico_2_indel_vcf)
# insilico_2_result_indels_charite

insilico_3_indels_charite <- "./insilico_3_indels_charite.vcf"
insilico_3_result_indels_charite <- metrics_calculator_indel(insilico_3_indels_charite,truth_insilico_3_indel_vcf)
# insilico_3_result_indels_charite


insilico_3_result_indels_bsc
insilico_3_result_indels_charite
insilico_3_result_indels_oicr
insilico_3_result_indels_hmf



insilico_3_indels_curie <- "./insilico_3_indels_curie.vcf"

insilico_3_result_indels_curie <- metrics_calculator_indel(insilico_3_indels_curie,truth_insilico_3_indel_vcf) 

insilico_3_result_indels_curie








