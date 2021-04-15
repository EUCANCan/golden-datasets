#!/usr/bin/env Rscript

args<-commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: compute_polym.r <inputList> <output> <sample name><bed polym> <minDP (optional, default=10)>", call.=FALSE)
}

inputTestFile <- args[1]
inputTruthFile<-args[2]
outputDir<-args[3]
sampleName<-args[4]

# Handle path & Output Names
outputFile=paste(outputDir,"/",sampleName,sep = "")

# Load data
test_df = read.table(inputTestFile, header = T, stringsAsFactors = F, sep = ",")
#head(test_df)
colnames(test_df)=c("chrom","pos","ref","alt")
#head(test_df)

truth_df = read.table(inputTruthFile, header = T, stringsAsFactors = F,sep = ",")
#head(truth_df)
colnames(test_df)=c("chrom","pos","ref","alt")

# Benchmark function
metrics_calculator_snv <- function(test_df,truth_df){
  test_df$mut <- paste0(test_df$chrom,"_",test_df$pos,"_",test_df$ref,"_",test_df$alt)
  truth_df$mut <- paste0(truth_df$chrom,"_",truth_df$pos,"_",truth_df$ref,"_",truth_df$alt)
  
  options(scipen=999)
  TP <- test_df[test_df$mut %in% truth_df$mut,]
  FP <- test_df[!test_df$mut %in% truth_df$mut,]
  FN <- truth_df[!truth_df$mut %in% test_df$mut,]
  Recall <- nrow(TP)/nrow(truth_df) 
  Precision <-nrow(TP)/(nrow(TP)+nrow(FP)) 
  F1_Score <- 2 * Precision * Recall / (Precision +Recall)
  
  metrics_X <- data.frame(round(Recall,3), round(Precision,3), round(F1_Score,3), nrow(TP), nrow(FP), nrow(FN))
  colnames(metrics_X) <- c("Recall","Precision","F1_Score","FP","TP","FN")
  metrics_X
  
  return(metrics_X)
}

results = metrics_calculator_snv(test_df, truth_df)
print(sampleName[1])
results

 #truth_insilico_1_snv_vcf <-"./truth.snvs.synthetic.challenge.set1.vcf"
 #truth_insilico_2_snv_vcf <-"./truth.snvs.synthetic.challenge.set2.vcf"
 #truth_insilico_3_snv_vcf <-"./truth.snvs.synthetic.challenge.set3.vcf"
 
 
 #insilico_1_snv_bsc <- "./insilico_1_snvs_bsc.vcf"
 #insilico_1_result_snv_bsc <- metrics_calculator_snv(insilico_1_snv_bsc,truth_insilico_1_snv_vcf)
 # insilico_1_result_snv_bsc
 
 #insilico_2_snv_bsc <- "./insilico_2_snvs_bsc.vcf"
 #insilico_2_result_snv_bsc <- metrics_calculator_snv(insilico_2_snv_bsc,truth_insilico_2_snv_vcf)
 # insilico_2_result_snv_bsc

 #insilico_3_snv_bsc <- "./insilico_3_snvs_bsc.vcf"
 #insilico_3_result_snv_bsc <- metrics_calculator_snv(insilico_3_snv_bsc,truth_insilico_3_snv_vcf)
 # insilico_3_result_snv_bsc




 #insilico_1_snv_hmf <- "./insilico_1_snvs_hmf.vcf"
 #insilico_1_result_snv_hmf <- metrics_calculator_snv(insilico_1_snv_hmf,truth_insilico_1_snv_vcf)
 # insilico_1_result_snv_hmf
 
 #insilico_2_snv_hmf <- "./insilico_2_snvs_hmf.vcf"
 #insilico_2_result_snv_hmf <- metrics_calculator_snv(insilico_2_snv_hmf,truth_insilico_2_snv_vcf)
 # insilico_2_result_snv_hmf
 
 #insilico_3_snv_hmf <- "./insilico_3_snvs_hmf.vcf"
 #insilico_3_result_snv_hmf <- metrics_calculator_snv(insilico_3_snv_hmf,truth_insilico_3_snv_vcf)
 # insilico_3_result_snv_hmf
 
 
 
 #insilico_1_snv_oicr <- "./insilico_1_snvs_oicr.vcf"
 #insilico_1_result_snv_oicr <- metrics_calculator_snv(insilico_1_snv_oicr,truth_insilico_1_snv_vcf)
 # insilico_1_result_snv_oicr
 
 #insilico_2_snv_oicr <- "./insilico_2_snvs_oicr.vcf"
 #insilico_2_result_snv_oicr <- metrics_calculator_snv(insilico_2_snv_oicr,truth_insilico_2_snv_vcf)
 # insilico_2_result_snv_oicr
 
 #insilico_3_snv_oicr <- "./insilico_3_snvs_oicr.vcf"
 #insilico_3_result_snv_oicr <- metrics_calculator_snv(insilico_3_snv_oicr,truth_insilico_3_snv_vcf)
 # insilico_3_result_snv_oicr
 
 #insilico_1_snv_charite <- "./insilico_1_snvs_charite.vcf"
 #insilico_1_result_snv_charite <- metrics_calculator_snv(insilico_1_snv_charite,truth_insilico_1_snv_vcf)
 # insilico_1_result_snv_charite
 
 #insilico_2_snv_charite <- "./insilico_2_snvs_charite.vcf"
 #insilico_2_result_snv_charite <- metrics_calculator_snv(insilico_2_snv_charite,truth_insilico_2_snv_vcf)
 # insilico_2_result_snv_charite
 
 #insilico_3_snv_charite <- "./insilico_3_snvs_charite.vcf"
 #insilico_3_result_snv_charite <- metrics_calculator_snv(insilico_3_snv_charite,truth_insilico_3_snv_vcf)
 # insilico_3_result_snv_charite
 
 #insilico_1_result_snv_bsc
 #insilico_1_result_snv_charite
 #insilico_1_result_snv_oicr
 #insilico_1_result_snv_hmf
 
 #insilico_2_result_snv_bsc
 #insilico_2_result_snv_charite
 #insilico_2_result_snv_oicr
 #insilico_2_result_snv_hmf
 
 #insilico_3_result_snv_bsc
 #insilico_3_result_snv_charite
 #insilico_3_result_snv_oicr
 #insilico_3_result_snv_hmf
 

 
 
 #DEVEL
 
 #insilico_1_snv_curie <- "./insilico_1_snvs_curie.vcf"
 #insilico_1_result_snv_curie <- metrics_calculator_snv(insilico_1_snv_curie,truth_insilico_1_snv_vcf)
 # insilico_1_result_snv_curie
 
 #insilico_2_snv_curie <- "./insilico_2_snvs_curie.vcf"
 #insilico_2_result_snv_curie <- metrics_calculator_snv(insilico_2_snv_curie,truth_insilico_2_snv_vcf)
 # insilico_2_result_snv_curie
 
 #insilico_3_snv_curie <- "./insilico_3_snvs_curie.vcf"
 #insilico_3_result_snv_curie <- metrics_calculator_snv(insilico_3_snv_curie,truth_insilico_3_snv_vcf)
 # insilico_3_result_snv_curie

 
 #insilico_1_result_snv_curie
 #insilico_2_result_snv_curie
 #insilico_3_result_snv_curie

