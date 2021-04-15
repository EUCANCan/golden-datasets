#New file to check SV metrics for DREAM datasets results


#Set your own directory
setwd('~/Descargas/SVinsilicos')

#LIBRARIES

load.lib<-c("stringr", "viridis", "data.table","ggplot2","gridExtra","UpSetR",
            "caret","catspec","devtools","tidyr","plyr","dplyr",
            "rpart","rpart.plot","rattle","e1071","openssl","sqldf")

install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) {install.packages(lib,dependencies=TRUE)}
sapply(load.lib,library,character=TRUE)



metrics_calculator_tab <- function(insilico_sv_tab, truth_insilico_sv_tab, SV_types_in_sample){
  
  
  #FORMATING TRUTH FILE
  truth_tab <- truth_insilico_sv_tab
  truth_tab_formated <- read.table(truth_tab)
  colnames(truth_tab_formated) <- c("#CHROM",  "POS" ,    "ID"  ,    "REF" ,    "ALT"  ,   "QUAL" ,   "FILTER" , "INFO" ,   "FORMAT",  "SPIKEIN")
  truth_tab_formated <- truth_tab_formated[!grepl("MSK", truth_tab_formated$ALT),]
  
  if("./truth.SV.synthetic.challenge.set3.vcf" %in% truth_insilico_sv_tab ){
    INS <- truth_tab_formated[grepl("INS", truth_tab_formated$ALT),]
    columns <- str_split_fixed(INS$INFO, ";", Inf)
    columns[,2] <- gsub("SVTYPE=", "", columns[,2])
    columns[,3] <- gsub("END=", "", columns[,3])
    columns[,4] <- gsub("SVLEN=", "", columns[,4])
    colnames(columns) <- c("no", "SV_TYPE", "END", "SVLEN")
    t <- as.data.frame(INS)
    t1 <- as.data.frame(columns)
    final_truth <- cbind(t, t1)
    # final_truth = final_truth[-1,]
    final_truth2 <- as.data.table(final_truth)
    final_truth2$SVLEN <- as.integer(as.character(final_truth2$SVLEN))
    final_truth2$SV_TYPE <- as.character(final_truth2$SV_TYPE)
    
    rest_variants<- truth_tab_formated[!grepl("INS", truth_tab_formated$ALT),]
    columns_rest <- str_split_fixed(rest_variants$INFO, ";", Inf)
    columns_rest_a <- columns_rest
    columns_rest_a[,1] <- columns_rest[,2]
    columns_rest_a[,3] <- gsub("END=", "",   columns_rest[,1])
    columns_rest_a[,4] <- gsub("SVLEN=", "", columns_rest[,3])
    columns_rest_a[,2] <- gsub("SVTYPE=", "",columns_rest[,4])
    colnames(columns_rest_a) <- c("no", "SV_TYPE", "END", "SVLEN")
    t_rest <- as.data.frame(rest_variants)
    t1_rest <- as.data.frame(columns_rest_a)
    final_truth_rest <- cbind(t_rest, t1_rest)
    # final_truth = final_truth[-1,]
    final_truth2_rest <- as.data.table(final_truth_rest)
    final_truth2_rest$SVLEN <- as.integer(as.character(final_truth2_rest$SVLEN))
    final_truth2_rest$SV_TYPE <- as.character(final_truth2_rest$SV_TYPE)
    
    truth_file <- rbind(final_truth2, final_truth2_rest)
    table(truth_file$SV_TYPE)
    
    truth_file$"#CHROM" <- as.character(truth_file$"#CHROM")
    
  } else {
    columns <- str_split_fixed(truth_tab_formated$INFO, ";", Inf)
    columns[,2] <- gsub("SVTYPE=", "", columns[,2])
    columns[,3] <- gsub("END=", "", columns[,3])
    columns[,4] <- gsub("SVLEN=", "", columns[,4])
    colnames(columns) <- c("no", "SV_TYPE", "END", "SVLEN")
    t <- as.data.frame(truth_tab_formated)
    t1 <- as.data.frame(columns)
    final_truth <- cbind(t, t1)
    # final_truth = final_truth[-1,]
    final_truth2 <- as.data.table(final_truth)
    final_truth2$SVLEN <- as.integer(as.character(final_truth2$SVLEN))
    final_truth2$SV_TYPE <- as.character(final_truth2$SV_TYPE)
    truth_file <- final_truth2
    truth_file$"#CHROM" <- as.character(truth_file$"#CHROM")
  }
  
  #FORMATING SV VCF FILES
  VC_tab <- insilico_sv_tab
  # VC_tab <- as.character(VC_tab)
  
  if( grepl( "vcf",VC_tab, fixed = TRUE)      ){
    VC_tab_formated <- read.table(VC_tab)
    colnames(VC_tab_formated) <- c("CHR1",  "POS" ,    "ID"  ,    "REF" ,    "ALT"  ,   "QUAL" ,   "FILTER" , "INFO" ,   "FORMAT",  "NORMAL", "TUMOUR")
    if(grepl( "oicr",VC_tab, fixed = TRUE) ){
      #oicr
      #debere añadir en un momento dado una columna SV_TYPE
      #mejor empezar por este que parece más facil y nos puede ayudar a hacer el otro
      gg <- as.character(VC_tab_formated$INFO)
      
      ggg <- str_split_fixed(gg, "SVCLASS=", Inf)
      
      zz <- ggg[,2]
      
      # table(zz)
      VC_tab_formated$INFO <-NULL
      VC_tab_formated$SV_TYPE <- zz
      
      VC_tab_formated$SV_TYPE <- gsub("inversion", "INV", VC_tab_formated$SV_TYPE)
      VC_tab_formated$SV_TYPE <-gsub("tandem-duplication", "DUP", VC_tab_formated$SV_TYPE)
      VC_tab_formated$SV_TYPE <-gsub("deletion", "DEL", VC_tab_formated$SV_TYPE)
      VC_tab_formated$SV_TYPE <-gsub("translocation", "TRANS", VC_tab_formated$SV_TYPE)
      VC_tab_formated$SV_TYPE <-gsub("insertion", "INS", VC_tab_formated$SV_TYPE)
      
      VC_tab_formated_g<- VC_tab_formated[grepl("_1", VC_tab_formated$ID),]
      VC_tab_formated_g$CHR1 <- as.character(VC_tab_formated_g$CHR1)
      
      SV_tab_file <- VC_tab_formated_g
      # table(SV_tab_file$SV_TYPE)
      # str(SV_tab_file)
      # print("vc file bien hecho")
    } else {
      
      
      VC_tab_formated <- VC_tab_formated[VC_tab_formated$FILTER == "PASS",]
      VC_tab_formated <-VC_tab_formated[grepl("*o",VC_tab_formated$ID),]
      VC_tab_formated$ID <- as.character(VC_tab_formated$ID)
      VC_tab_formated$CHR1 <- as.character(VC_tab_formated$CHR1)
      VC_tab_formated$ALT <- as.character(VC_tab_formated$ALT)
      VC_tab_formated$FILTER <- as.character(VC_tab_formated$FILTER)
      
      inv_a <-VC_tab_formated[grepl("*]$",VC_tab_formated$ALT),]
      inv_a$SV_TYPE <- "INV"
      inv_b <- VC_tab_formated[grepl("^\\Q[\\E",VC_tab_formated$ALT),]
      inv_b$SV_TYPE <- "INV"
      
      dup_b <- VC_tab_formated[grepl("^\\Q]\\E",VC_tab_formated$ALT),]
      dup_b$SV_TYPE <- "DUP"
      
      del_a <- VC_tab_formated[grepl("\\Q[\\E$",VC_tab_formated$ALT),]
      del_a$SV_TYPE <- "DEL"
      
      
      SV_tab_file <- rbind(inv_a, inv_b, del_a, dup_b)
    }
    
    
  }else {
    
    VC_tab_formated <- fread(VC_tab, header = T, sep = "\t", stringsAsFactors = F, na.strings = c(".", ".,."))
    
    if("MAPQ_DELLY2" %in% colnames(VC_tab_formated)) { 
      
      VC_tab_formated$SV_TYPE <- paste(VC_tab_formated$SVTYPE_DELLY2, VC_tab_formated$SVTYPE_BRASS, VC_tab_formated$SVTYPE_SVABA, sep="_")
      
      VC_tab_formated$SV_TYPE <- gsub("inversion", "INV", VC_tab_formated$SV_TYPE)
      VC_tab_formated$SV_TYPE <-gsub("tandem-duplication", "DUP", VC_tab_formated$SV_TYPE)
      VC_tab_formated$SV_TYPE <-gsub("deletion", "DEL", VC_tab_formated$SV_TYPE)
      VC_tab_formated$SV_TYPE <-gsub("translocation", "TRANS", VC_tab_formated$SV_TYPE)
      VC_tab_formated$SV_TYPE <-gsub("insertion", "INS", VC_tab_formated$SV_TYPE)
      
      VC_tab_formated$SV_TYPE <- gsub("NA_DEL", "DEL_DEL", VC_tab_formated$SV_TYPE)
      VC_tab_formated$SV_TYPE <- gsub("NA_INV", "INV_INV", VC_tab_formated$SV_TYPE)
      VC_tab_formated$SV_TYPE <- gsub("NA_DUP", "DUP_DUP", VC_tab_formated$SV_TYPE)
      VC_tab_formated$SV_TYPE <- gsub("NA_TRANS", "TRANS_TRANS", VC_tab_formated$SV_TYPE)
      VC_tab_formated$SV_TYPE <- gsub("NA_INS", "INS_INS", VC_tab_formated$SV_TYPE)
      
      new_VC <- str_split_fixed(VC_tab_formated$SV_TYPE, "_", Inf)
      new_VC <- as.data.frame(new_VC)
      colnames(new_VC) <- c("SV_TYPE2", "C", "D")
      VC_tab_formated2 <- cbind(VC_tab_formated, new_VC$SV_TYPE2)
      VC_tab_formated2$SV_TYPE <- NULL
      names(VC_tab_formated2)[names(VC_tab_formated2) == "V2"] <- "SV_TYPE"
      
      VC_tab_formated2_PASS <- VC_tab_formated2[VC_tab_formated2$CUSTOM_FILTER == "PASS"]
      VC_tab_formated2_PASS$SV_TYPE <- as.character(VC_tab_formated2_PASS$SV_TYPE)
      SV_tab_file <- VC_tab_formated2_PASS
      
      # table(VC_tab_formated2_PASS$SV_TYPE)
      
    }else{ #TSV from charite
      
      colnames(VC_tab_formated)[1] <- "CHR1"
      colnames(VC_tab_formated)[9] <- "SV_TYPE"
      new_VC <- as.data.frame(VC_tab_formated)
      SV_tab_file <- new_VC
      SV_tab_file$SV_TYPE  }
    
    
  }
  
  
  for (type in SV_types_in_sample ){ #c("DEL", "INV", "DUP")
    # type = "INV"
    # print(SV_types_in_sample)
    # print(type)
    
    # SV_tab_file[SV_tab_file$SV_TYPE == "DUP",]
    variants <- SV_tab_file[SV_tab_file$SV_TYPE == type,]
    # print(type)
    # print(paste("variants", type, nrow(variants)))
    # str(truth_file)
    truth_variants <- truth_file[truth_file$SV_TYPE == type,]
    # print(paste("truth variants", type, nrow(truth_variants)))
    merged_detected_all <- vector()
    status_vector_all <- vector() #"21"
    for (chr in c("1", "10" ,"11", "12" ,"13" ,"14" ,"15" ,"16", "17", "18" ,"19" , "2", "20" ,"21" ,"22" ,"3" , "4" , "5" , "6" , "7" , "8" , "9" , "X") ){   #unique(truth_variants$`#CHROM`) c("1", "10" ,"11", "12" ,"13" ,"14" ,"15" ,"16", "17", "18" ,"19" , "2", "20" ,"21" ,"22" ,"3" , "4" , "5" , "6" , "7" , "8" , "9" , "X")
      variants_chr <- variants[variants$CHR1 == chr,]
      truth_variants_chr <- truth_variants[truth_variants$"#CHROM" == chr,]
      status_vector <- numeric(nrow(variants_chr))
      merged_detected <- data.table(test = logical(nrow(truth_variants_chr)))
      # print(paste("New chr number: ",chr))
      # print(paste("Truth_variants here", type))
      # print(nrow(truth_variants_chr))
      # print(paste("variants_to_test", type))
      # print(nrow(variants_chr))
      if(nrow(variants_chr) > 0) {
        for (i in 1:nrow(variants_chr)) {
          # print(paste("nrow(variants): ",nrow(variants_chr) ))
          # print(paste("row numero: ", i, "chromosoma: ",chr))
          if (nrow(truth_variants_chr) > 0 ) {
            # detected <- ((truth_variants_chr[,2]-1000) <= c(variants_chr[i,2]) )  &  ( (truth_variants_chr[,2]+truth_variants_chr[,14]) >= c(variants_chr[i,2]) )
            detected <- ((truth_variants_chr[,2]-100) <= c(variants_chr[i,2]) )  &  ( (truth_variants_chr[,2]+100) >= c(variants_chr[i,2]) )
            #  print(detected)
            # print(paste("CHR number: ", chr))
            if(  TRUE %in% detected ) {
              status_vector[i] <- TRUE 
              # print(status_vector)
            } 
            else {status_vector[i] <- FALSE  }
            merged_detected <- cbind(merged_detected, detected) } #Esto al igual que los appends se podrían cambiar por la estrategia de status[i] usando esa misma i como numero de columna
          # print(merged_detected)
          # print(status_vector)
          else { status_vector[i] <- FALSE } } } 
      else {
        detected <- logical(nrow(truth_variants_chr))
        merged_detected <- cbind(merged_detected, detected)  }  
      
      status_vector_all <- append(status_vector_all, status_vector)
      # print(paste("status_vector_all", status_vector_all))
      # print(status_vector_all)
      # status_vector_all tiene que tener tantas posiciones como nrow(variants)
      merged_detected_all <- append(merged_detected_all, rowSums(merged_detected))
      # print(paste("merged_detected_all", merged_detected_all))
      # print("merged_detected_all")
      # print(merged_detected_all)
      #merged_detected_all tiene que tener tantas posiciones como nrow(truth_variants)
    }
    assign(paste0("TP_", type), variants[as.logical(status_vector_all),] )
    assign(paste0("FP_", type), variants[as.logical(abs(status_vector_all-1)),])
    assign(paste0("FN_", type), truth_variants[merged_detected_all == 0] )
    
    
    
  } 
  
  '%not_in%' <- Negate('%in%')
  if( "INS" %not_in% SV_types_in_sample){
    resultt <- list(FP_DEL, FP_INV, FP_DUP, TP_DEL, TP_DUP, TP_INV, FN_DEL, FN_DUP, FN_INV)
    final_table <- NULL
    final_table <- rbind( c(nrow(resultt[[4]]) / (nrow(resultt[[4]]) + nrow(resultt[[7]])),  nrow(resultt[[4]]) / (nrow(resultt[[4]]) + nrow(resultt[[1]]))   ,nrow(resultt[[1]]), nrow(resultt[[4]]),nrow(resultt[[7]]),  2 * (nrow(resultt[[4]]) / (nrow(resultt[[4]]) + nrow(resultt[[7]]))) * (nrow(resultt[[4]]) / (nrow(resultt[[4]]) + nrow(resultt[[1]]))) /          ((nrow(resultt[[4]]) / (nrow(resultt[[4]]) + nrow(resultt[[7]]))) + (nrow(resultt[[4]]) / (nrow(resultt[[4]]) + nrow(resultt[[1]]))))
    ), 
    c(nrow(resultt[[5]]) / (nrow(resultt[[5]]) + nrow(resultt[[8]])),  nrow(resultt[[5]]) / (nrow(resultt[[5]]) + nrow(resultt[[2]]))   ,nrow(resultt[[2]]), nrow(resultt[[5]]),nrow(resultt[[8]]),  2 * (nrow(resultt[[5]]) / (nrow(resultt[[5]]) + nrow(resultt[[8]]))) * (nrow(resultt[[5]]) / (nrow(resultt[[5]]) + nrow(resultt[[2]]))) /          ((nrow(resultt[[5]]) / (nrow(resultt[[5]]) + nrow(resultt[[8]]))) + (nrow(resultt[[5]]) / (nrow(resultt[[5]]) + nrow(resultt[[2]]))))
    ),
    c(nrow(resultt[[6]]) / (nrow(resultt[[6]]) + nrow(resultt[[9]])),  nrow(resultt[[6]]) / (nrow(resultt[[6]]) + nrow(resultt[[3]]))  ,nrow(resultt[[3]]), nrow(resultt[[6]]),nrow(resultt[[9]]), 2 * (nrow(resultt[[6]]) / (nrow(resultt[[6]]) + nrow(resultt[[9]]))) * (nrow(resultt[[6]]) / (nrow(resultt[[6]]) + nrow(resultt[[3]]))) /          ((nrow(resultt[[6]]) / (nrow(resultt[[6]]) + nrow(resultt[[9]]))) + (nrow(resultt[[6]]) / (nrow(resultt[[6]]) + nrow(resultt[[3]]))))
    ))
    colnames(final_table) <- c("Recall","Precision","FP","TP","FN","F1_Score")
    rownames(final_table) <- c("DEL", "DUP", "INV")
  } else {
    resultt <- list(FP_DEL, FP_INV, FP_DUP, TP_DEL, TP_DUP, TP_INV, FN_DEL, FN_DUP, FN_INV, FP_INS, TP_INS, FN_INS)
    
    final_table <- rbind( c(nrow(resultt[[4]]) / (nrow(resultt[[4]]) + nrow(resultt[[7]])),  nrow(resultt[[4]]) / (nrow(resultt[[4]]) + nrow(resultt[[1]]))   ,nrow(resultt[[1]]), nrow(resultt[[4]]),nrow(resultt[[7]]),  2 * (nrow(resultt[[4]]) / (nrow(resultt[[4]]) + nrow(resultt[[7]]))) * (nrow(resultt[[4]]) / (nrow(resultt[[4]]) + nrow(resultt[[1]]))) /          ((nrow(resultt[[4]]) / (nrow(resultt[[4]]) + nrow(resultt[[7]]))) + (nrow(resultt[[4]]) / (nrow(resultt[[4]]) + nrow(resultt[[1]]))))
    ), 
    c(nrow(resultt[[5]]) / (nrow(resultt[[5]]) + nrow(resultt[[8]])),  nrow(resultt[[5]]) / (nrow(resultt[[5]]) + nrow(resultt[[2]]))   ,nrow(resultt[[2]]), nrow(resultt[[5]]),nrow(resultt[[8]]),  2 * (nrow(resultt[[5]]) / (nrow(resultt[[5]]) + nrow(resultt[[8]]))) * (nrow(resultt[[5]]) / (nrow(resultt[[5]]) + nrow(resultt[[2]]))) /          ((nrow(resultt[[5]]) / (nrow(resultt[[5]]) + nrow(resultt[[8]]))) + (nrow(resultt[[5]]) / (nrow(resultt[[5]]) + nrow(resultt[[2]]))))
    ),
    c(nrow(resultt[[6]]) / (nrow(resultt[[6]]) + nrow(resultt[[9]])),  nrow(resultt[[6]]) / (nrow(resultt[[6]]) + nrow(resultt[[3]]))  ,nrow(resultt[[3]]), nrow(resultt[[6]]),nrow(resultt[[9]]), 2 * (nrow(resultt[[6]]) / (nrow(resultt[[6]]) + nrow(resultt[[9]]))) * (nrow(resultt[[6]]) / (nrow(resultt[[6]]) + nrow(resultt[[3]]))) /          ((nrow(resultt[[6]]) / (nrow(resultt[[6]]) + nrow(resultt[[9]]))) + (nrow(resultt[[6]]) / (nrow(resultt[[6]]) + nrow(resultt[[3]]))))
    ),
    c(nrow(resultt[[11]]) / (nrow(resultt[[11]]) + nrow(resultt[[12]])),  nrow(resultt[[11]]) / (nrow(resultt[[11]]) + nrow(resultt[[10]]))  ,nrow(resultt[[10]]), nrow(resultt[[11]]),nrow(resultt[[12]]), 2 * (nrow(resultt[[11]]) / (nrow(resultt[[11]]) + nrow(resultt[[12]]))) * (nrow(resultt[[11]]) / (nrow(resultt[[11]]) + nrow(resultt[[10]]))) /          ((nrow(resultt[[11]]) / (nrow(resultt[[11]]) + nrow(resultt[[12]]))) + (nrow(resultt[[11]]) / (nrow(resultt[[11]]) + nrow(resultt[[10]])))) 
    ))
    colnames(final_table) <- c("Recall","Precision","FP","TP","FN","F1_Score")
    rownames(final_table) <- c("DEL", "DUP", "INV", "INS")
    
  }
  
  return(final_table)
}  




#Set your own directory
setwd('~/Downloads')


truth_insilico_1_sv_tab <-"./truth.SV.synthetic.challenge.set1.vcf"
insilico_1_sv_tab_hmf <- "./insilico_1_sv_hmf.vcf"
SV_types_in_sample <- c("DEL", "INV", "DUP")
insilico_1_result_sv_hmf <- metrics_calculator_tab(insilico_1_sv_tab_hmf,truth_insilico_1_sv_tab, SV_types_in_sample)
# insilico_1_result_sv_hmf


truth_insilico_2_sv_tab <-"./truth.SV.synthetic.challenge.set2.vcf"
insilico_2_sv_tab_hmf <- "./insilico_2_sv_hmf.vcf"
SV_types_in_sample <- c("DEL", "INV", "DUP", "INS")
insilico_2_result_sv_hmf <- metrics_calculator_tab(insilico_2_sv_tab_hmf,truth_insilico_2_sv_tab, SV_types_in_sample)
# insilico_2_result_sv_hmf


truth_insilico_3_sv_tab <-"./truth.SV.synthetic.challenge.set3.vcf"
insilico_3_sv_tab_hmf <- "./insilico_3_sv_hmf.vcf"
SV_types_in_sample <- c("DEL", "INV", "DUP", "INS")
insilico_3_result_sv_hmf <- metrics_calculator_tab(insilico_3_sv_tab_hmf,truth_insilico_3_sv_tab, SV_types_in_sample)
# insilico_3_result_sv_hmf



truth_insilico_1_sv_tab <-"./truth.SV.synthetic.challenge.set1.vcf"
insilico_1_sv_tab_oicr <- "./insilico_1_sv_oicr.vcf"
SV_types_in_sample <- c("DEL", "INV", "DUP")
insilico_1_result_sv_oicr <- metrics_calculator_tab(insilico_1_sv_tab_oicr,truth_insilico_1_sv_tab, SV_types_in_sample)
# insilico_1_result_sv_oicr


truth_insilico_2_sv_tab <-"./truth.SV.synthetic.challenge.set2.vcf"
insilico_2_sv_tab_oicr <- "./insilico_2_sv_oicr.vcf"
SV_types_in_sample <- c("DEL", "INV", "DUP", "INS")
insilico_2_result_sv_oicr <- metrics_calculator_tab(insilico_2_sv_tab_oicr,truth_insilico_2_sv_tab, SV_types_in_sample)
# insilico_2_result_sv_oicr



truth_insilico_3_sv_tab <-"./truth.SV.synthetic.challenge.set3.vcf"
insilico_3_sv_tab_oicr <- "./insilico_3_sv_oicr.vcf"
SV_types_in_sample <- c("DEL", "INV", "DUP", "INS")
insilico_3_result_sv_oicr <- metrics_calculator_tab(insilico_3_sv_tab_oicr,truth_insilico_3_sv_tab, SV_types_in_sample)
# insilico_3_result_sv_oicr


insilico_1_sv_tab_charite <- "./insilico_1_sv_charite.tsv"
truth_insilico_1_sv_tab <-"./truth.SV.synthetic.challenge.set1.vcf"
SV_types_in_sample <- c("DEL", "INV", "DUP")
insilico_1_result_sv_charite <- metrics_calculator_tab(insilico_1_sv_tab_charite,truth_insilico_1_sv_tab, SV_types_in_sample)
# insilico_1_result_sv_charite

insilico_2_sv_tab_charite <- "./insilico_2_sv_charite.tsv"
truth_insilico_2_sv_tab <-"./truth.SV.synthetic.challenge.set2.vcf"
SV_types_in_sample <- c("DEL", "INV", "DUP", "INS")
insilico_2_result_sv_charite <- metrics_calculator_tab(insilico_2_sv_tab_charite,truth_insilico_2_sv_tab, SV_types_in_sample)
# insilico_2_result_sv_charite

insilico_3_sv_tab_charite <- "./insilico_3_sv_charite.tsv"
truth_insilico_3_sv_tab <-"./truth.SV.synthetic.challenge.set3.vcf"
SV_types_in_sample <- c("DEL", "INV", "DUP", "INS")
insilico_3_result_sv_charite <- metrics_calculator_tab(insilico_3_sv_tab_charite,truth_insilico_3_sv_tab, SV_types_in_sample)
# insilico_3_result_sv_charite



insilico_1_sv_tab_bsc <- "./insilico_1_sv_bsc.tsv"
truth_insilico_1_sv_tab <-"./truth.SV.synthetic.challenge.set1.vcf"
SV_types_in_sample <- c("DEL", "INV", "DUP")
insilico_1_result_sv_bsc <- metrics_calculator_tab(insilico_1_sv_tab_bsc,truth_insilico_1_sv_tab, SV_types_in_sample)
# insilico_1_result_sv_bsc

insilico_2_sv_tab_bsc <- "./insilico_2_sv_bsc.tsv"
truth_insilico_2_sv_tab <-"./truth.SV.synthetic.challenge.set2.vcf"
SV_types_in_sample <- c("DEL", "INV", "DUP", "INS")
insilico_2_result_sv_bsc <- metrics_calculator_tab(insilico_2_sv_tab_bsc,truth_insilico_2_sv_tab, SV_types_in_sample)
# insilico_2_result_sv_bsc

insilico_3_sv_tab_bsc <- "./insilico_3_sv_bsc.tsv"
truth_insilico_3_sv_tab <-"./truth.SV.synthetic.challenge.set3.vcf"
SV_types_in_sample <- c("DEL", "INV", "DUP", "INS")
insilico_3_result_sv_bsc <- metrics_calculator_tab(insilico_3_sv_tab_bsc,truth_insilico_3_sv_tab, SV_types_in_sample)
# insilico_3_result_sv_bsc

insilico_1_result_sv_bsc
insilico_1_result_sv_oicr
insilico_1_result_sv_charite
insilico_1_result_sv_hmf

insilico_2_result_sv_bsc
insilico_2_result_sv_oicr
insilico_2_result_sv_charite
insilico_2_result_sv_hmf

insilico_3_result_sv_bsc
insilico_3_result_sv_oicr
insilico_3_result_sv_charite
insilico_3_result_sv_hmf












# truth_tab <- truth_insilico_sv_tab
# truth_tab_formated <- read.table(truth_tab)
# colnames(truth_tab_formated) <- c("#CHROM",  "POS" ,    "ID"  ,    "REF" ,    "ALT"  ,   "QUAL" ,   "FILTER" , "INFO" ,   "FORMAT",  "SPIKEIN")
# truth_tab_formated <- truth_tab_formated[!grepl("MSK", truth_tab_formated$ALT),]
# 

truth.3.del.vcf <- read.table("./truth.3.del.vcf")
colnames(truth.3.del.vcf) <- c("#CHROM",  "POS" ,    "ID"  ,    "REF" ,    "ALT"  ,   "QUAL" ,   "FILTER" , "INFO" ,   "FORMAT",  "SPIKEIN")

truth.3.dup.vcf <- read.table("./truth.3.dup.vcf")
colnames(truth.3.dup.vcf) <- c("#CHROM",  "POS" ,    "ID"  ,    "REF" ,    "ALT"  ,   "QUAL" ,   "FILTER" , "INFO" ,   "FORMAT",  "SPIKEIN")

truth.3.ins.vcf <- read.table("./truth.3.ins.vcf")
colnames(truth.3.ins.vcf) <- c("#CHROM",  "POS" ,    "ID"  ,    "REF" ,    "ALT"  ,   "QUAL" ,   "FILTER" , "INFO" ,   "FORMAT",  "SPIKEIN")


curie_del.vcf <- read.table("./curie_del.vcf",sep = "\t",)
colnames(curie_del.vcf) <- c("#CHROM",  "POS" ,    "ID"  ,    "REF" ,    "ALT"  ,   "QUAL" ,   "FILTER" , "INFO" ,   "FORMAT",  "SPIKEIN", "ww")

curie_dup.vcf <- read.table("./curie_dup.vcf",sep = "\t",)
colnames(curie_dup.vcf) <- c("#CHROM",  "POS" ,    "ID"  ,    "REF" ,    "ALT"  ,   "QUAL" ,   "FILTER" , "INFO" ,   "FORMAT",  "SPIKEIN", "ww")

curie_ins.vcf <- read.table("./curie_ins.vcf",sep = "\t",)
colnames(curie_ins.vcf) <- c("#CHROM",  "POS" ,    "ID"  ,    "REF" ,    "ALT"  ,   "QUAL" ,   "FILTER" , "INFO" ,   "FORMAT",  "SPIKEIN", "ww")





((curie_del.vcf[1,1] == truth.3.del.vcf[1,1] ) & ((curie_del.vcf[1,2] >= truth.3.del.vcf[1,2]-100) &   (curie_del.vcf[1,2] <= truth.3.del.vcf[1,2]+100))  )
truth.3.del.vcf[1,2]-100
curie_del.vcf[1,2]

#DELETIONS
for (variant in 1:563){
  for (true_variant in 1:705){
    if ((curie_del.vcf[variant,1] == truth.3.del.vcf[true_variant,1]) & (  (curie_del.vcf[variant,2] >= truth.3.del.vcf[true_variant,2]-100) &   (curie_del.vcf[variant,2] <= truth.3.del.vcf[true_variant,2]+100))  ){
      curie_del.vcf[variant,11] <- "TP"
      truth.3.del.vcf[true_variant,10] <- "Founded"
      break
    }
    else {
      curie_del.vcf[variant,11] <- "FP"
    }
  }
}

#INSERTIONS
for (variant in 1:107){
  for (true_variant in 1:688){
    if ((curie_ins.vcf[variant,1] == truth.3.ins.vcf[true_variant,1]) & (  (curie_ins.vcf[variant,2] >= truth.3.ins.vcf[true_variant,2]-100) &   (curie_ins.vcf[variant,2] <= truth.3.ins.vcf[true_variant,2]+100))  ){
      curie_ins.vcf[variant,11] <- "TP"
      truth.3.ins.vcf[true_variant,10] <- "Founded"
      break
    }
    else {
      curie_ins.vcf[variant,11] <- "FP"
    }
  }
}
#DUPLICATIONS
for (variant in 1:635){
  for (true_variant in 1:746){
    if ((curie_dup.vcf[variant,1] == truth.3.dup.vcf[true_variant,1]) & (  (curie_dup.vcf[variant,2] >= truth.3.dup.vcf[true_variant,2]-100) &   (curie_dup.vcf[variant,2] <= truth.3.dup.vcf[true_variant,2]+100))  ){
      curie_dup.vcf[variant,11] <- "TP"
      truth.3.dup.vcf[true_variant,10] <- "Founded"
      break
    }
    else {
      curie_dup.vcf[variant,11] <- "FP"
    }
  }
}

#Calcular métricas con lo de arriba easy y poner en el documento
table(curie_del.vcf[,11])
table(truth.3.del.vcf[,10])
# FP  TP 
# 44 519 
# ./. Founded 
# 186     519 

# Recall = 519 / (106 + 519) = 0.8304
# Precision = 519 / (44 + 519) = 0.9218472
# F-Score = 2 * Precision * Recall / (Precision + Recall) = 2 * 0.8304 * 0.9218472 / ( 0.8304 + 0.9218472 ) = 0.8737374
table(curie_dup.vcf[,11])
table(truth.3.dup.vcf[,10])
# FP  TP 
# 21 614 
# ./. Founded 
# 132     614 

# Recall = 614 / (132 + 614) = 0.8230563
# Precision = 614 / (21 + 614) = 0.9669291
# F-Score = 2 * 0.8230563 * 0.9669291 / ( 0.8230563 + 0.9669291 ) = 0.8892107
table(curie_ins.vcf[,11])
table(truth.3.ins.vcf[,10])
# FP TP 
# 31 76 
# ./. Founded 
# 613      75 

# Recall = 75 / (613 + 75) = 0.1090116
# Precision = 75 / (31 + 75) = 0.7075472
# F-Score = 2 * 0.1090116 * 0.7075472 / ( 0.1090116 + 0.7075472 ) = 0.1889168






truth.3.inv.vcf <- read.table("./truth.3.inv.vcf")
colnames(truth.3.inv.vcf) <- c("#CHROM",  "POS" ,    "ID"  ,    "REF" ,    "ALT"  ,   "QUAL" ,   "FILTER" , "INFO" ,   "FORMAT",  "SPIKEIN")

curie_inv.vcf <- read.table("./curie_bnd.vcf",sep = "\t",)
colnames(curie_inv.vcf) <- c("#CHROM",  "POS" ,    "ID"  ,    "REF" ,    "ALT"  ,   "QUAL" ,   "FILTER" , "INFO" ,   "FORMAT",  "SPIKEIN", "ww")


inv_a <-curie_inv.vcf[grepl("*]$",curie_inv.vcf$ALT),]
inv_a$SV_TYPE <- "INV"
inv_b <- curie_inv.vcf[grepl("^\\Q[\\E",curie_inv.vcf$ALT),]
inv_b$SV_TYPE <- "INV"

curie_inv.vcf <- rbind(inv_a, inv_b)
#Creo que no es necesario juntarlos, asique divido entre 2 los FP y los TP



#INVERSIONS
for (variant in 1:3366){
  for (true_variant in 1:747){
    if ((curie_inv.vcf[variant,1] == truth.3.inv.vcf[true_variant,1]) & (  (curie_inv.vcf[variant,2] >= truth.3.inv.vcf[true_variant,2]-100) &   (curie_inv.vcf[variant,2] <= truth.3.inv.vcf[true_variant,2]+100))  ){
      curie_inv.vcf[variant,11] <- "TP"
      truth.3.inv.vcf[true_variant,10] <- "Founded"
      break
    }
    else {
      curie_inv.vcf[variant,11] <- "FP"
    }
  }
}

#Calcular métricas con lo de arriba easy y poner en el documento
table(curie_inv.vcf[,11])
table(truth.3.inv.vcf[,10])

# FP   TP 
# 2212 1154 
# 1106 577
# 1106 594
# ./. Founded 
# 153     594 

# Recall = 594 / (153 + 594) = 0.7951807
# Precision = 594 / (1106 + 594) = 0.3494118
# F-Score = 2 * Precision * Recall / (Precision + Recall) = 2 * 0.7951807 * 0.3494118 / ( 0.7951807 + 0.3494118 ) = 0.4854925








