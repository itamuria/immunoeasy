

library(ggplot2)
library(ggridges)
library(reshape2)
library(openxlsx)
library(limma)
library(UniProt.ws)
library(org.Hs.eg.db)
# library('DNABarcodes')
# library(alakazam)
library(RecordLinkage)

setwd("G:/Mi unidad/NASIR/data/Puros/")
source("from_vcf2df.R")

# gaur <- day_hour("day")
gaur <- 20200703

# Obtain vcf
pac <- c(1,4,8,15,25,39)



# After MHC ---------------------------------------------------------------


gaur <- 20200703

# Obtain vcf
pac <- c(1,4,8,15,25,39)
pac <- c(pac,5:7,9,13,14,16,17,20,26,28:37,41:43)


zein_karpetatik <- "G:/Mi unidad/NASIR/Belen/FiltradoDeAllRank/ClaseII"
zein_karpetara <- "G:/Mi unidad/NASIR/Belen/FiltradoDeAllRank/Resultados2021"
# pac <- 1
all_hla_sum <- NULL

for(pac_num in pac)
{
  first_filter <- grep(paste0("Nas_Pac", pac_num, "_"), dir(zein_karpetatik), value = T)
  second_filter <- grep("hla_core", first_filter, value = T, invert = T)
  # file_intermediate <- paste0(zein_karpetatik,"/Nas_Pac",pac_num,"_",gaur,"_R.RData")
  # save(selected, file = file_intermediate)
  # load(file_intermediate)
  selected <- openxlsx::read.xlsx(paste0(zein_karpetatik, "/",second_filter))

  # hla_clean <- gsub("HLA-","",selected$HLA)
  # # *
  # hla <- paste0(substr(hla_clean,1,1),substr(hla_clean,3,4))
  # selected$HLA_clean <- hla
  # selected$strong_weak <- ifelse(selected$Pr_Rank < 0.5, "Strong",
  #                                ifelse(selected$Pr_Rank >= 0.5 & selected$Pr_Rank < 2, "Weak",NA))
  #
  size <- nrow(selected)

  # summary_table <- c(summary_table, size, sum(df_hla2$Strong), sum(df_hla2$Weak))
  # selected$kmer <- paste0("K",nchar(selected$Peptide), "_", gsub("M_","", selected$Identity))

  # selected$lenWT <- nchar(selected$seq_mut)
  # selected$lenMUT <- nchar(selected$Peptide)

  # selected2 <- selected[, c("HLA","Peptide","lenMUT","seq_mut","lenWT","Core","GeneName","Change","kmer","Aff(nM)","Pr_Rank",
  #                           "wild_aff","wild_rank", "DAI","DAI_rank", "strong_weak")]
  #
  #
  # names(selected2)[c(2,4,8)] <- c("MutSeq","WildTSe1","Identity")
  # selected2 <- selected2[selected2$DAI != 1,]
  #
  # openxlsx::write.xlsx(selected2, file = paste0(zein_karpetara2,"Nas_Pac",pac_num,"_",gaur,"_excel",size,"_Without_selected.xlsx"))

  selected2 <- selected

  thla <- unique(selected2[!is.na(selected2$strong_weak), c("HLA", "Core","strong_weak")])
  thla2 <-  table(thla$HLA, thla$strong_weak)
  openxlsx::write.xlsx(thla2, file = paste0(zein_karpetara,"/Nas_Pac",pac_num,"_",gaur,"_hla_claseII.xlsx"))


  hla_core <- NULL
  for(sil in c(500, 100, 50))
  {
    seli <- selected2[selected2$`Aff(nM)` < sil, ]

    # seli500 <- selected[selected$`Aff(nM)` < 500, ]
    # seli100 <- selected[selected$`Aff(nM)` < 100, ]

    for(hh in unique(seli$HLA))
    {

      ha <- seli[seli$HLA == hh, ]

      for(core in unique(ha$Core))
      {
        net <- ha[ha$Core == core & ha$HLA == hh,]
        if(nrow(net) > 0)
        {

          hla_core <- c(hla_core, hh, core, sil, nrow(net))
        }
        print(paste("Fil ", sil, " - hla ", hh, "core ", core, " - Count ", nrow(net)))

        #                 ifelse(sum(net$strong_weak == "Strong", na.rm = TRUE) > 0, 1, 0),
        #                 ifelse(sum(net$strong_weak == "Weak", na.rm = TRUE) > 0, 1, 0))
      }
    }
  }

  df_hla <- data.frame(matrix(hla_core, byrow = TRUE, ncol = 4))
  names(df_hla) <- c("HLA","Core","Filtro","Counts")
  df_hla$Counts <- as.numeric(as.character(df_hla$Counts))
  df_hla$Filtro <- as.numeric(as.character(df_hla$Filtro))

  hla_core2 <- NULL
  for(sil in c(500, 100, 50))
  {
    seli2 <- df_hla[df_hla$Filtro == sil, ]

    for(hh in unique(seli2$HLA))
    {

      ha2 <- seli2[seli2$HLA == hh, ]
      hla_core2 <- c(hla_core2, hh, sil, nrow(ha2))
    }
  }

  df_hla3 <- data.frame(matrix(hla_core2, byrow = TRUE, ncol = 3))
  names(df_hla3) <- c("HLA","Filtro","Counts")

  m1 <- merge(df_hla3[df_hla3$Filtro == 500,], df_hla3[df_hla3$Filtro == 100,], by = "HLA", all = TRUE)
  m2 <- merge(m1, df_hla3[df_hla3$Filtro == 50,], by = "HLA", all.x = TRUE)
  m3 <- m2[,-c(2,4,6)]
  names(m3) <- c("HLA", "Fil500", "Fil100", "Fil50")
  m3$Fil500 <- as.numeric(as.character( m3$Fil500))
  m3$Fil100 <- as.numeric(as.character( m3$Fil100))
  m3$Fil50 <- as.numeric(as.character( m3$Fil50))

  all_hla_sum <- c(all_hla_sum, paste0("Pac",pac_num), sum(as.numeric(as.character(df_hla3$Counts))))

  openxlsx::write.xlsx(m3, file = paste0(zein_karpetara,"/Nas_Pac",pac_num,"_",gaur,"_hla_core_allrank50010050_claseII.xlsx"))

}

df_all_hla_sum <- data.frame(matrix(all_hla_sum, ncol = 2, byrow = TRUE))
names(df_all_hla_sum) <- c("Paciente", "Counts")
openxlsx::write.xlsx(df_all_hla_sum, file = paste0(zein_karpetara, "/Puro_Clase2_sum_all_hla.xlsx"))




