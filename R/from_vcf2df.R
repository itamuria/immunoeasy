#' From vcf to excel by mutect38__2df
#'
#' @param name_vcffile The name of the vcf file
#' @param excel If we want to export to excel (boolean)
#' @param excel_file Name of the output excel file.

#' @return data frame with ordered information
#' @export
#'
#' @examples
#' \dontrun{
#' mutect38__2df (name_vcffile)
#' }
#'
mutect38_2df <- function(name_vcffile, excel = FALSE, excel_file = "20200306_Mutect.xlsx")
{

  vcffile <- utils::read.table(name_vcffile)
  # separate format
  vcffile$V9 <- as.character(vcffile$V9)
  fox <- grep("FOXOG",vcffile$V9)
  pgt <- grep("PGT",vcffile$V9)

  table(vcffile$V9)

  vcffile$fox <- ifelse(rownames(vcffile)%in%fox,"Fox",NA)
  vcffile$pgt <- ifelse(rownames(vcffile)%in%pgt,"Pgt",NA)

  vcffile$group <- 0
  vcffile$group[vcffile$fox == "Fox" &  is.na(vcffile$pgt)] <- 1
  vcffile$group[vcffile$fox == "Fox" &  vcffile$pgt == "Pgt"] <- 2
  vcffile$group[is.na(vcffile$fox) &  is.na(vcffile$pgt)] <- 3
  vcffile$group[vcffile$pgt == "Pgt" &  is.na(vcffile$fox)] <- 4
  vcffile$ID <- 1:(nrow(vcffile))

  table(vcffile$group)
  name_group <- names(table(vcffile$group))

  # Extract names
  u1 <- unique(vcffile$V9[vcffile$fox == "Fox" &  vcffile$pgt == "Pgt"])
  # if(!is.na(u1))
  # {
  u1 <- u1[!is.na(u1)]
  names9 <- strsplit(u1,":")[[1]]
  len9 <- length(names9)

  # GROUP 1
  g1_10 <- as.character(vcffile$V10[vcffile$group==1])
  g1_11 <- as.character(vcffile$V11[vcffile$group==1])
  id_1 <- vcffile$ID[vcffile$group==1]

  col10 <- unlist(strsplit(g1_10,":"))
  col10b <- data.frame(matrix(col10,ncol = 9,byrow = TRUE))
  col10c <- cbind(col10b[,1:6],rep(NA,nrow(col10b)),rep(NA,nrow(col10b)),col10b[,7:9])
  names(col10c) <- names9

  col11 <- unlist(strsplit(g1_11,":"))
  col11b <- data.frame(matrix(col11,ncol = 9,byrow = TRUE))
  col11c <- cbind(col11b[,1:6],rep(NA,nrow(col11b)),rep(NA,nrow(col11b)),col11b[,7:9])
  names(col11c) <- names9

  # GROUP 2
  g2_10 <- as.character(vcffile$V10[vcffile$group==2])
  g2_11 <- as.character(vcffile$V11[vcffile$group==2])
  id_2 <- vcffile$ID[vcffile$group==2]

  col210 <- unlist(strsplit(g2_10,":"))
  col210b <- data.frame(matrix(col210,ncol = 11,byrow = TRUE))
  col210c <- col210b
  names(col210c) <- names9

  col211 <- unlist(strsplit(g2_11,":"))
  col211b <- data.frame(matrix(col211,ncol = 11,byrow = TRUE))
  col211c <- col211b
  names(col211c) <- names9

  # GROUP 3
  g3_10 <- as.character(vcffile$V10[vcffile$group==3])
  g3_11 <- as.character(vcffile$V11[vcffile$group==3])
  id_3 <- vcffile$ID[vcffile$group==3]

  col310 <- unlist(strsplit(g3_10,":"))
  col310b <- data.frame(matrix(col310,ncol = 8,byrow = TRUE))
  col310c <- cbind(col310b[,1:5],rep(NA,nrow(col310b)),rep(NA,nrow(col310b)),rep(NA,nrow(col310b)),col310b[,6:8])
  names(col310c) <- names9

  col311 <- unlist(strsplit(g3_11,":"))
  col311b <- data.frame(matrix(col311,ncol = 8,byrow = TRUE))
  col311c <- cbind(col311b[,1:5],rep(NA,nrow(col311b)),rep(NA,nrow(col311b)),rep(NA,nrow(col311b)),col311b[,6:8])
  names(col311c) <- names9

  # GROUP 4
  g4_10 <- as.character(vcffile$V10[vcffile$group==4])
  g4_11 <- as.character(vcffile$V11[vcffile$group==4])
  id_4 <- vcffile$ID[vcffile$group==4]

  col410 <- unlist(strsplit(g4_10,":"))
  col410b <- data.frame(matrix(col410,ncol = 10,byrow = TRUE))
  col410c <- cbind(col410b[,1:5],rep(NA,nrow(col410b)),col410b[,6:10])
  names(col410c) <- names9

  col411 <- unlist(strsplit(g4_11,":"))
  col411b <- data.frame(matrix(col411,ncol = 10,byrow = TRUE))
  col411c <- cbind(col411b[,1:5],rep(NA,nrow(col411b)),col411b[,6:10])
  names(col411c) <- names9

  id_t <- c(id_1,id_2,id_3,id_4)
  dfcol <- rbind(col10c,col210c,col310c,col410c)
  dfcol3 <- rbind(col11c,col211c,col311c,col411c)
  dfcol2 <- cbind(id_t,dfcol,dfcol3)
  dfcol2 <- dfcol2[order(dfcol2$id_t),]
  names(dfcol2)[13:23]<- paste0(names(dfcol2)[2:12],"_tum")

  ad1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$AD),",")),ncol = 2,byrow = TRUE))
  names(ad1) <- c("AD_ref","AD_alt")
  gss1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$QSS),",")),ncol = 2,byrow = TRUE))
  names(gss1) <- c("QSS_ref","QSS_alt")

  ad1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$AD_tum),",")),ncol = 2,byrow = TRUE))
  names(ad1b) <- c("AD_ref_tum","AD_alt_tum")
  gss1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$QSS_tum),",")),ncol = 2,byrow = TRUE))
  names(gss1b) <- c("QSS_ref_tum","QSS_alt_tum")

  # FINAL FORMAT
  formatfin <- cbind(dfcol2[,c(1,2)],ad1,dfcol2[,c(4:8)],gss1,dfcol2[,c(11,13)],ad1b,dfcol2[,c(15:20)],gss1b,dfcol2[,c(22,23)])

  ### INFO
  f8 <- vcffile$V8

  s1 <- tidyr::separate(data = vcffile, col = "V8", into = paste0("a",1:58), sep = ";")
  s2 <- s1[,c(8:65)]
  s2[grep("RPA",s2$a6),names(s2)] <- s2[grep("RPA",s2$a6),names(s2)[-c(6:8)]]

  snames <- unlist(strsplit(as.character(s2[1,]),"="))
  snames2 <- snames[seq(1,110,by = 2)]

  for(h in c(1:41,43:58))
  {
    print(h)
    s2[,h] <- gsub(paste0(snames2[h],"="),"",s2[,h])
    names(s2)[h] <- snames2[h]
  }
  names(s2)[42]<- "GERP"

  dftot <- cbind(vcffile[,1:7],s2,formatfin)
  names(dftot)[1:7] <- c("Chr","Pos","V3","Ref","Alt","V6","Filt")

  for(q in c(2,68:73,75:77,79:84,87:90))
  {
    dftot[,q] <- as.numeric(as.character(dftot[,q]))
  }

  dftot$Cob_DNA_nor <- dftot$AD_ref+dftot$AD_alt
  dftot$Cob_DNA_tum <- dftot$AD_ref_tum+dftot$AD_alt_tum

  dftot <- dftot[,!is.na(names(dftot))]

  if(excel) openxlsx::write.xlsx(dftot, file = excel_file)

  return(dftot)

}

#' From vcf to excel by Strelka
#'
#' @param name_vcffile The name of the vcf file
#' @param excel If we want to export to excel (boolean)
#' @param excel_file Name of the output excel file.

#' @return data frame with ordered information
#' @export
#'
#' @examples
#' \dontrun{
#' strelka_2_excel (name_vcffile)
#' }
#'
# cuidado con GT, no quitar o cambiar codigo despues de cambiar
strelka_2_excel <- function(name_vcffile, excel = FALSE, excel_file = "20200306_strelka.xlsx")
{
  vcffile <- utils::read.table(name_vcffile)
  # separate format
  vcffile$V9 <- as.character(vcffile$V9)
  vcffile$ID <- 1:(nrow(vcffile))

  # Extract names
  u1 <- unique(vcffile$V9)
  u1 <- u1[!is.na(u1)]
  names9 <- strsplit(u1,":")[[1]]
  # names9 <- names9[-1]
  # names9 <- names9[-1]
  len9 <- length(names9)

  # GROUP 1
  g1_10 <- as.character(vcffile$V10)
  # nchar(g1_10)

  g1_11 <- as.character(vcffile$V11)
  id_1 <- vcffile$ID

  col10 <- unlist(strsplit(g1_10,":"))
  col10b <- data.frame(matrix(col10,ncol = len9,byrow = TRUE))
  col10c <- col10b
  names(col10c) <- names9

  col11 <- unlist(strsplit(g1_11,":"))
  col11b <- data.frame(matrix(col11,ncol = len9,byrow = TRUE))
  col11c <- col11b
  names(col11c) <- names9


  id_t <- id_1
  dfcol <- rbind(col10c)
  dfcol3 <- rbind(col11c)
  dfcol2 <- cbind(id_t,dfcol,dfcol3)
  dfcol2 <- dfcol2[order(dfcol2$id_t),]
  nnc <- ncol(dfcol2)
  nnc2 <- (nnc-1) / 2
  nnc3 <- nnc2 * 2 + 1
  names(dfcol2)[(nnc2+2):nnc3]<- paste0(names(dfcol2)[2:(nnc2+1)],"_2")

  ad1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$AU),",")),ncol = 2,byrow = TRUE))
  names(ad1) <- c("AU_T1","AU_T2")
  gss1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$CU),",")),ncol = 2,byrow = TRUE))
  names(gss1) <- c("CU_T1","CU_T2")
  gss2 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$GU),",")),ncol = 2,byrow = TRUE))
  names(gss2) <- c("GU_T1","GU_T2")
  gss3 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$TU),",")),ncol = 2,byrow = TRUE))
  names(gss3) <- c("TU_T1","TU_T2")

  ad1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$AU_2),",")),ncol = 2,byrow = TRUE))
  names(ad1b) <- c("AU_T1_2","AU_T2_2")
  gss1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$CU_2),",")),ncol = 2,byrow = TRUE))
  names(gss1b) <- c("CU_T1_2","CU_T2_2")
  gss2b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$GU_2),",")),ncol = 2,byrow = TRUE))
  names(gss2b) <- c("GU_T1_2","GU_T2_2")
  gss3b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$TU),",")),ncol = 2,byrow = TRUE))
  names(gss3b) <- c("TU_T1_2","TU_T2_2")

  # FINAL FORMAT
  formatfin <- cbind(dfcol2[,c(1:5)],ad1,gss1,gss2,gss3,dfcol2[,c(10:13)],ad1b,gss1b,gss2b,gss3b)

  ### INFO
  f8 <- vcffile$V8

  s1 <- tidyr::separate(data = vcffile, col = "V8", into = paste0("a",1:60), sep = ";")
  s2 <- s1[,c(8:67)]
  snames <- unlist(strsplit(as.character(s2[1,]),"="))
  snames2 <- snames[c(1,seq(2,118,by = 2))]

  for(h in c(2:48,50:60))
  {
    print(h)
    s2[,h] <- gsub(paste0(snames2[h],"="),"",s2[,h])
    names(s2)[h] <- snames2[h]
  }

  dftot <- cbind(vcffile[,1:7],s2,formatfin)
  names(dftot)[1:7] <- c("Chr","Pos","V3","Ref","Alt","V6","Filt")

  # plot(as.numeric(dftot$AU_T1),as.numeric(dftot$AU_T2))
  # cor(as.numeric(dftot$AU_T1),as.numeric(dftot$AU_T2))
  # cor.test(as.numeric(dftot$AU_T1),as.numeric(dftot$AU_T2))
  #
  # plot(as.numeric(dftot$GU_T1),as.numeric(dftot$GU_T2))
  # cor(as.numeric(dftot$GU_T1),as.numeric(dftot$GU_T2))
  # cor.test(as.numeric(dftot$GU_T1),as.numeric(dftot$GU_T2))

  for(y in 69:ncol(dftot))
  {
    dftot[,y] <- as.numeric(as.character(dftot[,y]))
  }


  dftot$Ref_T1 <- 0
  dftot$Alt_T1 <- 0
  dftot$Ref_T2 <- 0
  dftot$Alt_T2 <- 0

  for(q in 1:nrow(dftot))
  {
    print(paste0(q,"/",nrow(dftot)))
    dftot$Ref_T1[q] <- ifelse(dftot$Ref[q]=="A",dftot[q,"AU_T1"],
                              ifelse(dftot$Ref[q]=="C",dftot[q,"CU_T1"],
                                     ifelse(dftot$Ref[q]=="G",dftot[q,"GU_T1"],
                                            ifelse(dftot$Ref[q]=="T",dftot[q,"TU_T1"],NA))))

    dftot$Ref_T2[q]  <- ifelse(dftot$Ref[q]=="A",dftot[q,"AU_T2"],
                               ifelse(dftot$Ref[q]=="C",dftot[q,"CU_T2"],
                                      ifelse(dftot$Ref[q]=="G",dftot[q,"GU_T2"],
                                             ifelse(dftot$Ref[q]=="T",dftot[q,"TU_T2"],NA))))

    dftot$Alt_T1[q]  <- ifelse(dftot$Alt[q]=="A",dftot[q,"AU_T1"],
                               ifelse(dftot$Alt[q]=="C",dftot[q,"CU_T1"],
                                      ifelse(dftot$Alt[q]=="G",dftot[q,"GU_T1"],
                                             ifelse(dftot$Alt[q]=="T",dftot[q,"TU_T1"],NA))))

    dftot$Alt_T2[q]  <- ifelse(dftot$Alt[q]=="A",dftot[q,"AU_T2"],
                               ifelse(dftot$Alt[q]=="C",dftot[q,"CU_T2"],
                                      ifelse(dftot$Alt[q]=="G",dftot[q,"GU_T2"],
                                             ifelse(dftot$Alt[q]=="T",dftot[q,"TU_T2"],NA))))
  }




  dftot$Ref_T1_2 <- 0
  dftot$Alt_T1_2 <- 0
  dftot$Ref_T2_2 <- 0
  dftot$Alt_T2_2 <- 0

  for(q in 1:nrow(dftot))
  {
    print(paste0(q,"/",nrow(dftot)))
    dftot$Ref_T1_2[q] <- ifelse(dftot$Ref[q]=="A",dftot[q,"AU_T1_2"],
                                ifelse(dftot$Ref[q]=="C",dftot[q,"CU_T1_2"],
                                       ifelse(dftot$Ref[q]=="G",dftot[q,"GU_T1_2"],
                                              ifelse(dftot$Ref[q]=="T",dftot[q,"TU_T1_2"],NA))))

    dftot$Ref_T2_2[q]  <- ifelse(dftot$Ref[q]=="A",dftot[q,"AU_T2_2"],
                                 ifelse(dftot$Ref[q]=="C",dftot[q,"CU_T2_2"],
                                        ifelse(dftot$Ref[q]=="G",dftot[q,"GU_T2_2"],
                                               ifelse(dftot$Ref[q]=="T",dftot[q,"TU_T2_2"],NA))))

    dftot$Alt_T1_2[q]  <- ifelse(dftot$Alt[q]=="A",dftot[q,"AU_T1_2"],
                                 ifelse(dftot$Alt[q]=="C",dftot[q,"CU_T1_2"],
                                        ifelse(dftot$Alt[q]=="G",dftot[q,"GU_T1_2"],
                                               ifelse(dftot$Alt[q]=="T",dftot[q,"TU_T1_2"],NA))))

    dftot$Alt_T2_2[q]  <- ifelse(dftot$Alt[q]=="A",dftot[q,"AU_T2_2"],
                                 ifelse(dftot$Alt[q]=="C",dftot[q,"CU_T2_2"],
                                        ifelse(dftot$Alt[q]=="G",dftot[q,"GU_T2_2"],
                                               ifelse(dftot$Alt[q]=="T",dftot[q,"TU_T2_2"],NA))))
  }



  dftote <- dftot
  dftot[dftot$NT=="hom",c("Ref_T1","Ref_T2","Ref_T1_2","Ref_T2_2")] <- dftote[dftote$NT=="hom",c("Alt_T1","Alt_T2","Alt_T1_2","Alt_T2_2")]
  dftot[dftot$NT=="hom",c("Alt_T1","Alt_T2","Alt_T1_2","Alt_T2_2")] <- dftote[dftote$NT=="hom",c("Ref_T1","Ref_T2","Ref_T1_2","Ref_T2_2")]


  dftot$VAF_T1 <- dftot$Alt_T1_2/(dftot$Alt_T1_2+dftot$Ref_T1_2)
  dftot$VAF_T2 <- dftot$Alt_T2_2/(dftot$Alt_T2_2+dftot$Ref_T2_2)

  # ff <- dftot$Ref_T1-dftot$Alt_T1
  if(excel) openxlsx::write.xlsx(dftot, file = excel_file)
  # aggregate(ff,by=list(dftot$NT),FUN=function(x) c(min(x),mean(x),max(x)))
  return(dftot)

}

#' From vcf to excel by Strelka fast version
#'
#' @param name_vcffile The name of the vcf file
#' @param excel If we want to export to excel (boolean)
#' @param excel_file Name of the output excel file.
#'
#' @return data frame with ordered information
#' @export
#'
#' @examples
#' \dontrun{
#' strelka_2_excel (name_vcffile)
#' }
#'
strelka_2_excel_snvs_fast <- function(name_vcffile, excel = FALSE, excel_file = "20200306_strelka_snvs.xlsx")
{
  vcffile <- utils::read.table(name_vcffile)
  # separate format
  vcffile$V9 <- as.character(vcffile$V9)
  vcffile$ID <- 1:(nrow(vcffile))

  # Extract names
  u1 <- unique(vcffile$V9)
  u1 <- u1[!is.na(u1)]
  names9 <- strsplit(u1,":")[[1]]
  # names9 <- names9[-1]
  # names9 <- names9[-1]
  len9 <- length(names9)

  # GROUP 1
  g1_10 <- as.character(vcffile$V10)
  nchar(g1_10)

  g1_11 <- as.character(vcffile$V11)
  id_1 <- vcffile$ID

  col10 <- unlist(strsplit(g1_10,":"))
  col10b <- data.frame(matrix(col10,ncol = len9,byrow = TRUE))
  col10c <- col10b
  names(col10c) <- names9

  col11 <- unlist(strsplit(g1_11,":"))
  col11b <- data.frame(matrix(col11,ncol = len9,byrow = TRUE))
  col11c <- col11b
  names(col11c) <- names9


  id_t <- id_1
  dfcol <- rbind(col10c)
  dfcol3 <- rbind(col11c)
  dfcol2 <- cbind(id_t,dfcol,dfcol3)
  dfcol2 <- dfcol2[order(dfcol2$id_t),]
  nnc <- ncol(dfcol2)
  nnc2 <- (nnc-1) / 2
  nnc3 <- nnc2 * 2 + 1
  names(dfcol2)[(nnc2+2):nnc3]<- paste0(names(dfcol2)[2:(nnc2+1)],"_2")

  ad1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$AU),",")),ncol = 2,byrow = TRUE))
  names(ad1) <- c("AU_T1_normal","AU_T2_normal")
  gss1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$CU),",")),ncol = 2,byrow = TRUE))
  names(gss1) <- c("CU_T1_normal","CU_T2_normal")
  gss2 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$GU),",")),ncol = 2,byrow = TRUE))
  names(gss2) <- c("GU_T1_normal","GU_T2_normal")
  gss3 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$TU),",")),ncol = 2,byrow = TRUE))
  names(gss3) <- c("TU_T1_normal","TU_T2_normal")

  ad1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$AU_2),",")),ncol = 2,byrow = TRUE))
  names(ad1b) <- c("AU_T1_tumor","AU_T2_tumor")
  gss1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$CU_2),",")),ncol = 2,byrow = TRUE))
  names(gss1b) <- c("CU_T1_tumor","CU_T2_tumor")
  gss2b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$GU_2),",")),ncol = 2,byrow = TRUE))
  names(gss2b) <- c("GU_T1_tumor","GU_T2_tumor")
  gss3b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$TU_2),",")),ncol = 2,byrow = TRUE))
  names(gss3b) <- c("TU_T1_tumor","TU_T2_tumor")

  # FINAL FORMAT
  formatfin <- cbind(dfcol2[,c(1:5)],ad1,gss1,gss2,gss3,dfcol2[,c(10:13)],ad1b,gss1b,gss2b,gss3b)

  ### INFO
  f8 <- vcffile$V8

  s1 <- tidyr::separate(data = vcffile, col = "V8", into = paste0("a",1:60), sep = ";")
  s2 <- s1[,c(8:67)]
  snames <- unlist(strsplit(as.character(s2[1,]),"="))
  snames2 <- snames[c(1,seq(2,118,by = 2))]

  for(h in c(2:48,50:60))
  {
    print(h)
    s2[,h] <- gsub(paste0(snames2[h],"="),"",s2[,h])
    names(s2)[h] <- snames2[h]
  }

  dftot <- cbind(vcffile[,1:7],s2,formatfin)
  names(dftot)[1:7] <- c("Chr","Pos","V3","Ref","Alt","V6","Filt")

  # plot(as.numeric(dftot$AU_T1),as.numeric(dftot$AU_T2))
  # cor(as.numeric(dftot$AU_T1),as.numeric(dftot$AU_T2))
  # cor.test(as.numeric(dftot$AU_T1),as.numeric(dftot$AU_T2))
  #
  # plot(as.numeric(dftot$GU_T1),as.numeric(dftot$GU_T2))
  # cor(as.numeric(dftot$GU_T1),as.numeric(dftot$GU_T2))
  # cor.test(as.numeric(dftot$GU_T1),as.numeric(dftot$GU_T2))

  for(y in 29:ncol(dftot))
  {
    dftot[,y] <- as.numeric(as.character(dftot[,y]))
  }


  dftot$Ref_T1_normal <- 0
  dftot$Alt_T1_normal <- 0
  dftot$Ref_T2_normal <- 0
  dftot$Alt_T2_normal <- 0

  dftot$Ref_T1_normal <- ifelse(dftot$Ref=="A",dftot[,"AU_T1_normal"],
                                ifelse(dftot$Ref=="C",dftot[,"CU_T1_normal"],
                                       ifelse(dftot$Ref=="G",dftot[,"GU_T1_normal"],
                                              ifelse(dftot$Ref=="T",dftot[,"TU_T1_normal"],NA))))

  dftot$Ref_T2_normal  <- ifelse(dftot$Ref=="A",dftot[,"AU_T2_normal"],
                                 ifelse(dftot$Ref=="C",dftot[,"CU_T2_normal"],
                                        ifelse(dftot$Ref=="G",dftot[,"GU_T2_normal"],
                                               ifelse(dftot$Ref=="T",dftot[,"TU_T2_normal"],NA))))

  dftot$Alt_T1_normal  <- ifelse(dftot$Alt=="A",dftot[,"AU_T1_normal"],
                                 ifelse(dftot$Alt=="C",dftot[,"CU_T1_normal"],
                                        ifelse(dftot$Alt=="G",dftot[,"GU_T1_normal"],
                                               ifelse(dftot$Alt=="T",dftot[,"TU_T1_normal"],NA))))

  dftot$Alt_T2_normal  <- ifelse(dftot$Alt=="A",dftot[,"AU_T2_normal"],
                                 ifelse(dftot$Alt=="C",dftot[,"CU_T2_normal"],
                                        ifelse(dftot$Alt=="G",dftot[,"GU_T2_normal"],
                                               ifelse(dftot$Alt=="T",dftot[,"TU_T2_normal"],NA))))


  dftot$Ref_T1_tumor <- 0
  dftot$Alt_T1_tumor <- 0
  dftot$Ref_T2_tumor <- 0
  dftot$Alt_T2_tumor <- 0

  dftot$Ref_T1_tumor <- ifelse(dftot$Ref=="A",dftot[,"AU_T1_tumor"],
                               ifelse(dftot$Ref=="C",dftot[,"CU_T1_tumor"],
                                      ifelse(dftot$Ref=="G",dftot[,"GU_T1_tumor"],
                                             ifelse(dftot$Ref=="T",dftot[,"TU_T1_tumor"],NA))))

  dftot$Ref_T2_tumor  <- ifelse(dftot$Ref=="A",dftot[,"AU_T2_tumor"],
                                ifelse(dftot$Ref=="C",dftot[,"CU_T2_tumor"],
                                       ifelse(dftot$Ref=="G",dftot[,"GU_T2_tumor"],
                                              ifelse(dftot$Ref=="T",dftot[,"TU_T2_tumor"],NA))))

  dftot$Alt_T1_tumor  <- ifelse(dftot$Alt=="A",dftot[,"AU_T1_tumor"],
                                ifelse(dftot$Alt=="C",dftot[,"CU_T1_tumor"],
                                       ifelse(dftot$Alt=="G",dftot[,"GU_T1_tumor"],
                                              ifelse(dftot$Alt=="T",dftot[,"TU_T1_tumor"],NA))))

  dftot$Alt_T2_tumor  <- ifelse(dftot$Alt=="A",dftot[,"AU_T2_tumor"],
                                ifelse(dftot$Alt=="C",dftot[,"CU_T2_tumor"],
                                       ifelse(dftot$Alt=="G",dftot[,"GU_T2_tumor"],
                                              ifelse(dftot$Alt=="T",dftot[,"TU_T2_tumor"],NA))))



  dftote <- dftot
  dftot[dftot$NT=="hom",c("Ref_T1_normal","Ref_T2_normal","Ref_T1_tumor","Ref_T2_tumor")] <- dftote[dftote$NT=="hom",c("Alt_T1_normal","Alt_T2_normal","Alt_T1_tumor","Alt_T2_tumor")]
  dftot[dftot$NT=="hom",c("Alt_T1_normal","Alt_T2_normal","Alt_T1_tumor","Alt_T2_tumor")] <- dftote[dftote$NT=="hom",c("Ref_T1_normal","Ref_T2_normal","Ref_T1_tumor","Ref_T2_tumor")]


  dftot$VAF_T1 <- dftot$Alt_T1_tumor/(dftot$Alt_T1_tumor+dftot$Ref_T1_tumor)
  dftot$VAF_T2 <- dftot$Alt_T2_tumor/(dftot$Alt_T2_tumor+dftot$Ref_T2_tumor)

  # ff <- dftot$Ref_T1-dftot$Alt_T1

  # aggregate(ff,by=list(dftot$NT),FUN=function(x) c(min(x),mean(x),max(x)))

  if(excel) openxlsx::write.xlsx(dftot, file = excel_file)

  return(dftot)

}

#' From vcf to excel by Strelka for indels
#'
#' @param name_vcffile The name of the vcf file
#' @param excel If we want to export to excel (boolean)
#' @param excel_file Name of the output excel file.
#'
#' @return data frame with ordered information
#' @export
#'
#' @examples
#' \dontrun{
#' strelka_2_excel (name_vcffile)
#' }
#'
strelka_indels_2_excel <- function(name_vcffile, excel = FALSE, excel_file = "20200306_strelka_indels.xlsx")
{
  vcffile <- utils::read.table(name_vcffile)
  # separate format
  vcffile$V9 <- as.character(vcffile$V9)
  vcffile$ID <- 1:(nrow(vcffile))

  # Extract names
  u1 <- unique(vcffile$V9)
  u1 <- u1[!is.na(u1)]
  names9 <- strsplit(u1,":")[[1]]
  # names9 <- names9[-1]
  len9 <- length(names9)

  # GROUP 1
  g1_10 <- as.character(vcffile$V10)
  g1_11 <- as.character(vcffile$V11)
  id_1 <- vcffile$ID

  col10 <- unlist(strsplit(g1_10,":"))
  col10b <- data.frame(matrix(col10,ncol = len9,byrow = TRUE))
  col10c <- col10b
  names(col10c) <- names9

  col11 <- unlist(strsplit(g1_11,":"))
  col11b <- data.frame(matrix(col11,ncol = len9,byrow = TRUE))
  col11c <- col11b
  names(col11c) <- names9


  id_t <- id_1
  dfcol <- rbind(col10c)
  dfcol3 <- rbind(col11c)
  dfcol2 <- cbind(id_t,dfcol,dfcol3)
  dfcol2 <- dfcol2[order(dfcol2$id_t),]
  nnc <- ncol(dfcol2)
  nnc2 <- (nnc-1) / 2
  nnc3 <- nnc2 * 2 + 1
  names(dfcol2)[(nnc2+2):nnc3]<- paste0(names(dfcol2)[2:(nnc2+1)],"_2")

  ad1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$TAR),",")),ncol = 2,byrow = TRUE))
  names(ad1) <- c("TAR_T1","TAR_T2")
  gss1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$TIR),",")),ncol = 2,byrow = TRUE))
  names(gss1) <- c("TIR_T1","TIR_T2")
  gss2 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$TOR),",")),ncol = 2,byrow = TRUE))
  names(gss2) <- c("TOR_T1","TOR_T2")

  ad1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$TAR_2),",")),ncol = 2,byrow = TRUE))
  names(ad1b) <- c("TAR_T1_2","TAR_T2_2")
  gss1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$TIR_2),",")),ncol = 2,byrow = TRUE))
  names(gss1b) <- c("TIR_T1_2","TIR_T2_2")
  gss2b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$TOR_2),",")),ncol = 2,byrow = TRUE))
  names(gss2b) <- c("TOR_T1_2","TOR_T2_2")


  # FINAL FORMAT
  formatfin <- cbind(dfcol2[,c(1:3)],ad1,gss1,gss2,dfcol2[,c(7:12)],ad1b,gss1b,gss2b)

  ### INFO
  f8 <- vcffile$V8

  s1 <- tidyr::separate(data = vcffile, col = "V8", into = paste0("a",1:60), sep = ";")
  s2 <- s1[,c(8:67)]
  snames <- unlist(strsplit(as.character(s2[1,]),"="))
  snames2 <- snames[c(1,seq(2,118,by = 2))]

  for(h in c(2:49,51:60))
  {
    print(h)
    s2[,h] <- gsub(paste0(snames2[h],"="),"",s2[,h])
    names(s2)[h] <- snames2[h]
  }

  dftot <- cbind(vcffile[,1:7],s2,formatfin)
  names(dftot)[1:7] <- c("Chr","Pos","V3","Ref","Alt","V6","Filt")

  dftot$TIR_T1_2 <- as.numeric(as.character(dftot$TIR_T1_2))
  dftot$TAR_T1_2 <- as.numeric(as.character(dftot$TAR_T1_2))

  dftot$vaf <- dftot$TIR_T1_2/(dftot$TIR_T1_2+dftot$TAR_T1_2)

  naw <- which(is.na(names(dftot)))
  if(length(naw)>0)
  {
    dftot <- dftot[,-naw]
  }

  pos1 <- which(names(dftot)=="DP")

  for(g in pos1:ncol(dftot))
  {
    dftot[,g] <- as.numeric(as.character(dftot[,g]))
  }

  if(excel) openxlsx::write.xlsx(dftot, file = excel_file)

  return(dftot)

}

#' From vcf to excel by somatic sniper
#'
#' @param name_vcffile The name of the vcf file
#' @param excel If we want to export to excel (boolean)
#' @param excel_file Name of the output excel file.
#'
#' @return data frame with ordered information
#' @export
#'
#' @examples
#' \dontrun{
#' strelka_2_excel (name_vcffile)
#' }
#'
somaticsniper_2_excel <- function(name_vcffile, excel = FALSE, excel_file = "20200306_somaticsniper.xlsx")
{
  vcffile <- utils::read.table(name_vcffile)
  # separate format
  vcffile$V9 <- as.character(vcffile$V9)
  vcffile$ID <- 1:(nrow(vcffile))

  # Extract names
  u1 <- unique(vcffile$V9)
  u1 <- u1[!is.na(u1)]
  names9 <- strsplit(u1,":")[[1]]
  len9 <- length(names9)

  # GROUP 1
  g1_10 <- as.character(vcffile$V10)
  g1_11 <- as.character(vcffile$V11)
  id_1 <- vcffile$ID

  col10 <- unlist(strsplit(g1_10,":"))
  col10b <- data.frame(matrix(col10,ncol = len9,byrow = TRUE))
  col10c <- col10b
  names(col10c) <- names9

  col11 <- unlist(strsplit(g1_11,":"))
  col11b <- data.frame(matrix(col11,ncol = len9,byrow = TRUE))
  col11c <- col11b
  names(col11c) <- names9


  id_t <- id_1
  dfcol <- rbind(col10c)
  dfcol3 <- rbind(col11c)
  dfcol2 <- cbind(id_t,dfcol,dfcol3)
  dfcol2 <- dfcol2[order(dfcol2$id_t),]
  nnc <- ncol(dfcol2)
  nnc2 <- (nnc-1) / 2
  nnc3 <- nnc2 * 2 + 1
  names(dfcol2)[(nnc2+2):nnc3]<- paste0(names(dfcol2)[2:(nnc2+1)],"_2")

  ad1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$DP4),",")),ncol = 4,byrow = TRUE))
  names(ad1) <- c("Norm_DP_ref_forw","Norm_DP_ref_rev","Norm_DP_alt_forw","Norm_DP_alt_rev")
  gss1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$BCOUNT),",")),ncol = 4,byrow = TRUE))
  names(gss1) <- c("NormCOUNT_A","NormCOUNT_C","NormCOUNT_G","NormCOUNT_T")

  ad1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$DP4_2),",")),ncol = 4,byrow = TRUE))
  names(ad1b) <- c("Tum_DP_ref_forw","Tum_DP_ref_rev","Tum_DP_alt_forw","Tum_DP_alt_rev")
  gss1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$BCOUNT_2),",")),ncol = 4,byrow = TRUE))
  names(gss1b) <- c("TumCOUNT_A","TumCOUNT_C","TumCOUNT_G","TumCOUNT_T")

  for(t in 1:4)
  {
    ad1[,t] <- as.numeric(as.character(ad1[,t]))
    ad1b[,t] <- as.numeric(as.character(ad1b[,t]))
  }

  alt1 <- ad1$Norm_DP_alt_forw+ad1$Norm_DP_alt_rev
  alt2 <- ad1b$Tum_DP_alt_forw+ad1b$Tum_DP_alt_rev

  dfcol2$DP <- as.numeric(as.character(dfcol2$DP))
  dfcol2$DP_2 <- as.numeric(as.character(dfcol2$DP_2))


  vaf1 <- alt1/dfcol2$DP
  vaf2 <- alt2/dfcol2$DP_2

  # FINAL FORMAT
  formatfin <- cbind(dfcol2[,c(1:4)],ad1,alt1,vaf1,gss1,dfcol2[,c(7:17)],ad1b,alt2,vaf2,gss1b,dfcol2[,c(20:27)])

  ### INFO
  f8 <- vcffile$V8

  s1 <- tidyr::separate(data = vcffile, col = "V8", into = paste0("a",1:48), sep = ";")
  s2 <- s1[,c(8:58)]
  snames <- unlist(strsplit(as.character(s2[1,]),"="))
  snames2 <- snames[c(seq(2,94,by = 2))]

  for(h in c(1:35,37:47))
  {
    print(h)
    s2[,h+1] <- gsub(paste0(snames2[h],"="),"",s2[,h+1])
    names(s2)[h+1] <- snames2[h]
  }

  names(s2)[37] <- "GERP"
  # s2[,37] <- gsub(paste0("GERP++_RS=.","",s2[,h+1]))
  dftot <- cbind(vcffile[,1:7],s2,formatfin)
  names(dftot)[1:7] <- c("Chr","Pos","V3","Ref","Alt","V6","Filt")


  nakendu <- sapply(dftot, function(e) {
    len1 <- length(e)
    pr90 <- len1 * 90 / 100
    pr90 < sum(is.na(e))
  })
  dftot <- dftot[,!nakendu]

  if(excel) openxlsx::write.xlsx(dftot, file = excel_file)

  return(dftot)

}

#' From vcf to excel by varscan
#'
#' @param name_vcffile The name of the vcf file
#' @param excel If we want to export to excel (boolean)
#' @param excel_file Name of the output excel file.
#'
#' @return data frame with ordered information
#' @export
#'
#' @examples
#' \dontrun{
#' strelka_2_excel (name_vcffile)
#' }
#'
varscan_2_excel <- function(name_vcffile, excel = FALSE, excel_file = "20200306_varscan.xlsx")
{
  vcffile <- utils::read.table(name_vcffile)
  # separate format
  vcffile$V9 <- as.character(vcffile$V9)
  vcffile$ID <- 1:(nrow(vcffile))

  # Extract names
  u1 <- unique(vcffile$V9)
  u1 <- u1[!is.na(u1)]
  names9 <- strsplit(u1,":")[[1]]
  len9 <- length(names9)

  # GROUP 1
  g1_10 <- as.character(vcffile$V10)
  g1_11 <- as.character(vcffile$V11)
  id_1 <- vcffile$ID

  col10 <- unlist(strsplit(g1_10,":"))
  col10b <- data.frame(matrix(col10,ncol = len9,byrow = TRUE))
  col10c <- col10b
  names(col10c) <- names9

  col11 <- unlist(strsplit(g1_11,":"))
  col11b <- data.frame(matrix(col11,ncol = len9,byrow = TRUE))
  col11c <- col11b
  names(col11c) <- names9


  id_t <- id_1
  dfcol <- rbind(col10c)
  dfcol3 <- rbind(col11c)
  dfcol2 <- cbind(id_t,dfcol,dfcol3)
  dfcol2 <- dfcol2[order(dfcol2$id_t),]
  nnc <- ncol(dfcol2)
  nnc2 <- (nnc-1) / 2
  nnc3 <- nnc2 * 2 + 1
  names(dfcol2)[(nnc2+2):nnc3]<- paste0("Tum_",names(dfcol2)[2:(nnc2+1)])

  ad1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$DP4),",")),ncol = 4,byrow = TRUE))
  names(ad1) <- c("Nor_Ref_ST_p","Nor_Ref_ST_n","Nor_Alt_ST_p","Nor_Alt_ST_n")
  # gss1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$BCOUNT),",")),ncol = 4,byrow = TRUE))
  # names(gss1) <- c("BCOUNT_ref","BCOUNT_alt","BCOUNT_3","BCOUNT_4")

  ad1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$Tum_DP4),",")),ncol = 4,byrow = TRUE))
  names(ad1b) <- c("Tum_Ref_ST_p","Tum_Ref_ST_n","Tum_Alt_ST_p","Tum_Alt_ST_n")
  # gss1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$BCOUNT_2),",")),ncol = 4,byrow = TRUE))
  # names(gss1b) <- c("BCOUNT_ref_2","BCOUNT_alt_2","BCOUNT_2_3","BCOUNT_2_4")


  # FINAL FORMAT
  formatfin <- cbind(dfcol2[,c(1:7)],ad1,dfcol2[,c(9:14)],ad1b)

  ### INFO
  f8 <- vcffile$V8

  s1 <- tidyr::separate(data = vcffile, col = "V8", into = paste0("a",1:53), sep = ";")
  s2 <- s1[,c(8,10:60)]
  # s2 <- s2[,-2]
  snames <- unlist(strsplit(as.character(s2[1,]),"="))
  snames2 <- snames[c(seq(1,106,by = 2))]

  azkeNA <- which(is.na(snames2))[1]-1
  for(h in c(1:azkeNA))
  {
    print(h)
    if(length(grep("GERP\\+\\+_RS", snames2[h]))>0)
    {
      s2[,h] <- gsub("GERP\\+\\+_RS","GERP_RS",s2[,h])
      snames2[h] <- gsub("GERP\\+\\+_RS","GERP_RS",snames2[h])
    }

    s2[,h] <- gsub(paste0(snames2[h],"="),"",s2[,h])
    names(s2)[h] <- snames2[h]
  }

  kendu <- sapply(s2,function(x)
  {
    nr <- length(x)
    sr <- sum(is.na(x))
    nr!=sr
  })

  tt <- table(kendu)
  if(length(tt)>1) s2 <- s2[,kendu]

  names(s2)[which(names(s2)=="DP")] <- "DP_Tum_Nor"

  # names(s2)[42] <- "GERP"
  # s2[,37] <- gsub(paste0("GERP++_RS=.","",s2[,h+1]))
  dftot <- cbind(vcffile[,1:7],s2,formatfin)
  names(dftot)[1:7] <- c("Chr","Pos","V3","Ref","Alt","V6","Filt")

  names(dftot)[which(names(dftot)=="RD")] <- "Nor_Ref_D"
  names(dftot)[which(names(dftot)=="DP")] <- "Nor_DP"
  names(dftot)[which(names(dftot)=="AD")] <- "Nor_Alt_D"

  names(dftot)[which(names(dftot)=="Tum_RD")] <- "Tum_Ref_D"
  names(dftot)[which(names(dftot)=="Tum_AD")] <- "Tum_Alt_D"

  w_names <- which(names(dftot)%in%c("Pos","DP_Tum_Nor","SSC","GPV","SPV","id_t","Nor_DP","Nor_Ref_D","Nor_Alt_D","Nor_Ref_ST_p" ,"Nor_Ref_ST_n", "Nor_Alt_ST_p", "Nor_Alt_ST_n",
                                     "Tum_DP","Tum_Ref_D","Tum_Alt_D","Tum_Ref_ST_p","Tum_Ref_ST_n","Tum_Alt_ST_p","Tum_Alt_ST_n"))

  for(z in w_names)
  {
    dftot[,z] <- as.numeric(as.character(dftot[,z]))
  }

  dftot$Tum_Ratio_Alt_st <- dftot$Tum_Alt_ST_p/dftot$Tum_Alt_ST_n

  dftot$FREQ <- as.character(dftot$FREQ)
  dftot$FREQ <- as.numeric(gsub("%","",dftot$FREQ))
  dftot$Tum_FREQ <- as.character(dftot$Tum_FREQ)
  dftot$Tum_FREQ <- as.numeric(gsub("%","",dftot$Tum_FREQ))

  dftot$Tum_FREQ <- dftot$Tum_FREQ/100
  dftot$FREQ <- dftot$FREQ/100

  names(dftot)[which(names(dftot)=="FREQ")] <- "Nor_VAF"
  names(dftot)[which(names(dftot)=="Tum_FREQ")] <- "Tum_VAF"

  if(excel) openxlsx::write.xlsx(dftot, file = excel_file)


  return(dftot)

}

