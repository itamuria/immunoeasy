from_vcf2excel_germinal <- function(name_vcffile = "exome_pass_v011.vcf",
                                    rdatafile = "20201125_vcf_MB14802_s11",save = TRUE)
{
  # setwd("Exoma/EXOMA_EAA001_V300056632_L02/")
  # dir()
  # name_vcffile = "exome_pass_v011.vcf"

  vcffile <- utils::read.table(name_vcffile)


  vcffile$V9 <- as.character(vcffile$V9)

  names9 <- strsplit(vcffile$V9,":")[[1]]
  len9 <- length(names9)

  col10 <- unlist(strsplit(as.character(vcffile$V10),":"))
  col10b <- data.frame(matrix(col10,ncol = len9,byrow = TRUE))
  names(col10b) <- names9

  col10b$AD <- as.character(col10b$AD)
  col10b$PL <- as.character(col10b$PL)

  sp2 <- strsplit(col10b$AD, ",")
  AD1 <-  unlist(lapply(sp2, `[[`, 1))
  AD2 <-  unlist(lapply(sp2, `[[`, 2))

  sp3 <- strsplit(col10b$PL, ",")
  PL1 <-  unlist(lapply(sp3, `[[`, 1))
  PL2 <-  unlist(lapply(sp3, `[[`, 2))
  PL3 <-  unlist(lapply(sp3, `[[`, 3))

  col10_fin <- cbind(col10b, AD1, AD2, PL1, PL2,PL3)

  # Read header
  d <- readLines(name_vcffile, n = 1000)
  rm_lines <- grep("INFO=",d)
  d <- d[rm_lines]

  d <- gsub("##INFO=<ID=","",d)
  sp2 <- strsplit(d, ",")
  info_names <-  unlist(lapply(sp2, `[[`, 1))

  vcffile[,info_names] <- ""

  # Funtzioa
  vcffile$V8 <- as.character(vcffile$V8)
  sp3 <- strsplit(vcffile$V8, ";")

  locname <- function(x, info)
  {
    all_string <- x[grep(paste0("^",info,"="),x)]
    value_only <- gsub(paste0(info,"="), "", all_string, )
    value_only <- value_only[1]
    # vcffile[h,info_names[g]] <- value_only[1]
    if(length(value_only) == 0) value_only <- NA
    return(value_only)
  }

  # Loop


  for(g in 1:length(info_names))
  {
    print(paste0(g, "/", length(info_names)))
    vcffile[,info_names[g]] <- unlist(lapply(sp3, locname, info = info_names[g]))

  }

  # un1 <- unlist(lapply(sp3, locname, info = info_names[g]))
  # table(un1)

  # Biltzen

  dftot <- cbind(vcffile[,-c(8:10)],col10_fin)
  names(dftot)[1:7] <- c("Chr","Pos","ID","Ref","Alt","QUAL","Filt")

  numvar <- c("Pos","ID","DP","esp6500siv2_ea",
              "SIFT_score","Polyphen2_HDIV_score","Polyphen2_HVAR_score","LRT_score",
              "MutationTaster_score","MutationAssessor_score","FATHMM_score","PROVEAN_score",
              "VEST3_score","CADD_raw","CADD_phred","DANN_score","fathmm-MKL_coding_score",
              "MetaSVM_score","MetaLR_score","integrated_fitCons_score",
              "integrated_confidence_value", "GERPpp_RS",
              "phyloP7way_vertebrate","phyloP20way_mammalian","phastCons7way_vertebrate",
              "phastCons20way_mammalian","SiPhy_29way_logOdds","ExAC_ALL","ExAC_AFR","ExAC_AMR",
              "ExAC_EAS","ExAC_FIN","ExAC_NFE","ExAC_OTH","ExAC_SAS",
              "DP","AD1","AD2","PL1",
              "PL2","PL3")

  for(q in numvar)
  {
    dftot[,q] <- as.numeric(as.character(dftot[,q]))
  }

  # for(q in c(8:23,31,seq(38,64,2),66:71,75,76,78:82))
  # {
  #   dftot[,q] <- as.numeric(as.character(dftot[,q]))
  # }

  if(save)
  {
    save(dftot, file = paste0(rdatafile, ".RData"))
    write.csv(dftot, file = paste0(rdatafile, ".csv"))
  }
}


# Somaticas
from_vcf2excel_somatica <- function(name_vcffile = "exome_pass_v011.vcf",
                                    rdatafile = "20201125_vcf_Ex77801_s10",save = TRUE)
{
  # setwd("Exoma/EXOMA_EAA001_V300056632_L02/")
  # dir()
  # name_vcffile = "exome_pass_v011.vcf"

  vcffile <- utils::read.table(name_vcffile)


  vcffile$V9 <- as.character(vcffile$V9)

  OBQRC <- grep("OBQRC",vcffile$V9)
  PGT <- grep("PGT",vcffile$V9)

  table(vcffile$V9)

  vcffile$OBQRC <- ifelse(rownames(vcffile)%in%OBQRC,"OBQRC",NA)
  vcffile$PGT <- ifelse(rownames(vcffile)%in%PGT,"PGT",NA)

  vcffile$group <- 0
  vcffile$group[vcffile$OBQRC == "OBQRC" &  is.na(vcffile$PGT)] <- 1
  vcffile$group[vcffile$OBQRC == "OBQRC" &  vcffile$PGT == "PGT"] <- 2
  vcffile$group[is.na(vcffile$OBQRC) &  is.na(vcffile$PGT)] <- 3
  vcffile$group[vcffile$PGT == "PGT" &  is.na(vcffile$OBQRC)] <- 4
  vcffile$ID <- 1:(nrow(vcffile))

  table(vcffile$group)
  name_group <- names(table(vcffile$group))

  # Extract names
  u1 <- unique(vcffile$V9[vcffile$OBQRC == "OBQRC" &  vcffile$PGT == "PGT"])
  # if(!is.na(u1))
  # {
  u1 <- u1[!is.na(u1)]
  names9 <- strsplit(u1,":")[[1]]
  len9 <- length(names9)

  # GROUP 1
  g1_10 <- as.character(vcffile$V10[vcffile$group==1])
  # g1_11 <- as.character(vcffile$V11[vcffile$group==1])
  id_1 <- vcffile$ID[vcffile$group==1]

  col10 <- unlist(strsplit(g1_10,":"))
  col10b <- data.frame(matrix(col10,ncol = 12,byrow = TRUE))
  col10c <- cbind(col10b,rep(NA,nrow(col10b)),rep(NA,nrow(col10b)),rep(NA,nrow(col10b)))
  names(col10c) <- names9


  # GROUP 2
  g2_10 <- as.character(vcffile$V10[vcffile$group==2])
  id_2 <- vcffile$ID[vcffile$group==2]

  col210 <- unlist(strsplit(g2_10,":"))
  col210b <- data.frame(matrix(col210,ncol = 15,byrow = TRUE))
  col210c <- col210b
  names(col210c) <- names9



  # GROUP 3
  g3_10 <- as.character(vcffile$V10[vcffile$group==3])
  id_3 <- vcffile$ID[vcffile$group==3]

  col310 <- unlist(strsplit(g3_10,":"))
  col310b <- data.frame(matrix(col310,ncol = 8,byrow = TRUE))
  col310c <- cbind(col310b[,1:8],data.frame(matrix(NA, nrow = nrow(col310b),ncol = 7, byrow = TRUE)))
  names(col310c) <- names9


  # GROUP 4
  g4_10 <- as.character(vcffile$V10[vcffile$group==4])
  id_4 <- vcffile$ID[vcffile$group==4]

  col410 <- unlist(strsplit(g4_10,":"))
  col410b <- data.frame(matrix(col410,ncol = 11,byrow = TRUE))
  col410c <- cbind(col410b[,1:8],data.frame(matrix(NA, nrow = nrow(col410b),ncol = 4, byrow = TRUE)),col410b[,9:11])
  names(col410c) <- names9


  id_t <- c(id_1,id_2,id_3,id_4)
  dfcol <- rbind(col10c,col210c,col310c,col410c)
  # dfcol3 <- rbind(col11c,col211c,col311c,col411c)
  dfcol2 <- cbind(id_t,dfcol)
  dfcol2 <- dfcol2[order(dfcol2$id_t),]
  # names(dfcol2)[13:23]<- paste0(names(dfcol2)[2:12],"_tum")

  dfcol2$AD <- as.character(dfcol2$AD)
  dfcol2$F1R2 <- as.character(dfcol2$F1R2)
  dfcol2$F2R1 <- as.character(dfcol2$F2R1)

  add <- strsplit(dfcol2$AD, ",")
  nadd <- unlist(lapply(add, length))
  ad1a <-  as.numeric(unlist(lapply(add, `[[`, 1)))
  ad1b <-  as.numeric(unlist(lapply(add, `[[`, 2)))
  ad <- data.frame(ad1a, ad1b, nadd)
  names(ad) <- c("AD_ref","AD_alt", "AD_n")

  F1R2d <- strsplit(dfcol2$F1R2, ",")
  nF1R2 <- unlist(lapply(F1R2d, length))
  F1R21a <-  as.numeric(unlist(lapply(F1R2d, `[[`, 1)))
  F1R21b <-  as.numeric(unlist(lapply(F1R2d, `[[`, 2)))
  F1R2 <- data.frame(F1R21a, F1R21b, nF1R2)
  names(F1R2) <- c("F1R2_ref","F1R2_alt","F1R2_n")

  F2R1d <- strsplit(dfcol2$F2R1, ",")
  nF2R1 <- unlist(lapply(F2R1d, length))
  F2R11a <-  as.numeric(unlist(lapply(F2R1d, `[[`, 1)))
  F2R11b <-  as.numeric(unlist(lapply(F2R1d, `[[`, 2)))
  F2R1 <- data.frame(F2R11a, F2R11b, nF2R1)
  names(F2R1) <- c("F2R1_ref","F2R1_alt", "F2R1_n")


  # FINAL FORMAT
  formatfin <- cbind(dfcol2[,c(1,3)],ad,dfcol2[,c(4:6)],F1R2,dfcol2[ ,7],F2R1,dfcol2[,c(8:ncol(dfcol2))])


  # Read header
  d <- readLines(name_vcffile, n = 1000)
  rm_lines <- grep("INFO=",d)
  d <- d[rm_lines]

  d <- gsub("##INFO=<ID=","",d)
  sp2 <- strsplit(d, ",")
  info_names <-  unlist(lapply(sp2, `[[`, 1))

  vcffile[,info_names] <- ""

  # Funtzioa
  vcffile$V8 <- as.character(vcffile$V8)
  sp3 <- strsplit(vcffile$V8, ";")

  locname <- function(x, info)
  {
    all_string <- x[grep(paste0("^",info,"="),x)]
    value_only <- gsub(paste0(info,"="), "", all_string, )
    value_only <- value_only[1]
    # vcffile[h,info_names[g]] <- value_only[1]
    if(length(value_only) == 0) value_only <- NA
    return(value_only)
  }

  # Loop


  for(g in 1:length(info_names))
  {
    print(paste0(g, "/", length(info_names)))
    vcffile[,info_names[g]] <- unlist(lapply(sp3, locname, info = info_names[g]))

  }

  # un1 <- unlist(lapply(sp3, locname, info = info_names[g]))
  # table(un1)

  # Biltzen

  dftot <- cbind(vcffile[,-c(8:10)],formatfin)
  names(dftot)[1:7] <- c("Chr","Pos","ID","Ref","Alt","QUAL","Filt")

  numvar <- c("Pos","group","ID","CONTQ","DP","ECNT","GERMQ","esp6500siv2_ea",
              "SIFT_score","Polyphen2_HDIV_score","Polyphen2_HVAR_score","LRT_score",
              "MutationTaster_score","MutationAssessor_score","FATHMM_score","PROVEAN_score",
              "VEST3_score","CADD_raw","CADD_phred","DANN_score","fathmm-MKL_coding_score",
              "MetaSVM_score","MetaLR_score","integrated_fitCons_score",
              "integrated_confidence_value", "GERPpp_RS",
              "phyloP7way_vertebrate","phyloP20way_mammalian","phastCons7way_vertebrate",
              "phastCons20way_mammalian","SiPhy_29way_logOdds","ExAC_ALL","ExAC_AFR","ExAC_AMR",
              "ExAC_EAS","ExAC_FIN","ExAC_NFE","ExAC_OTH","ExAC_SAS",
              "AD_ref","AD_alt","AD_n","AF","DP","F1R2_ref","F1R2_alt","F1R2_n",
              "F2R1_ref","F2R1_alt","F2R1_n",
              "OBF","OBP","OBQ","OBQRC")

  for(q in numvar)
  {
    dftot[,q] <- as.numeric(as.character(dftot[,q]))
  }

  if(save)
  {
    save(dftot, file = paste0(rdatafile, ".RData"))
    write.csv(dftot, file = paste0(rdatafile, ".csv"))
  }
}
