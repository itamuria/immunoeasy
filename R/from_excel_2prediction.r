#' Extract the aminoacid change
#'
#' @description Having a character vector with several words, extract the right change
#'
#' @param char A character variable with the information to be extracted
#'
#' @return The right change with first aminoacid, position and last aminoacid
#' @export
#'
#' @examples
#' \dontrun{
#' extact_aa_change ("A203R")
#' }
#'
extact_aa_change <- function(char)
{
  t2 <- strsplit(t,"p.")

  len_ch <- length(t2[[1]])
  ch_all <- NULL
  for(h in 2:(len_ch - 1))
  {
    ch1 <- t2[[1]][h]
    t3 <- strsplit(ch1,",|;")[[1]][1]
    ch_all <- c(ch_all, t3)
  }
  ch1 <- t2[[1]][len_ch]
  t3 <- strsplit(ch1,",|;")[[1]][1]
  ch_all <- c(ch_all, t3)
  ch_all <- unique(ch_all)

  out <- list()

  for(k in 1:length(ch_all))
  {
    primera <- substr(ch_all[k],1,1)
    pos_ult <- nchar(ch_all[k])
    ultimo <- substr(ch_all[k],pos_ult,pos_ult)
    numero <- substr(ch_all[k],2,(pos_ult-1))
    out[[k]] <- c(primera, numero, ultimo)
  }

  return(out)
}


#' Extract the sequences
#'
#' @description Having the right aminoacid change, it extract the needed sequences. For each gene we will generate several sequences
#'
#' @param df_all_fin A data frame with an specific structure with GeneName, aminoacid change
#' @param km The size of the sequence that we want to obtain
#' @return The specific sequences of the mutated and not mutated genes
#' @export
#'
#' @examples
#' \dontrun{
#' extracting_from_uniprot (df_all_fin)
#' }
#'

extracting_from_uniprot <- function(df_all_fin, km = c(8:12))
{
  dat <- df_all_fin[df_all_fin$HowManyVC > 1, ]

  head(dat)
  dat$seq_wt <- ""
  dat$seq_mut <- ""
  dat$pos <- 0
  dat$seq_len <- 0

  # checking symbol
  names(dat)[which(names(dat)=="Gene_name")] <- "GenName"
  names(dat)[which(names(dat)=="cDNA_Protein_change")] <- "AAChange.refGene"
  # dat$GenName

  dat$alias <- ""
  for(d in 1:nrow(dat))
  {
    print(dat$GenName[d])
    if(length(grep("\\\\",dat$GenName[d])) == 1)
    {
      dat$GenName[d] <- strsplit(dat$GenName[d], "\\\\")[[1]][1]
    }
    dat$alias[d] <- alias2Symbol(alias = dat$GenName[d], species = "Hs", expand.symbols = FALSE)
  }

  ## Loop

  for(a in 1:nrow(dat)) #
  {
    print("--------------------------------------------------------------------------------------------------------------------------")
    print(paste0("Paciente: NAS", pac_num, "---", a, "/", nrow(dat)))

    # if(pac_num == 4 & a%in%c(36)) break()

    t <- gene_name <- change <- NULL
    our_seq <- our_pos <- our_len <- NA

    t <- dat$AAChange.refGene[a]
    gene_name <- dat$alias[a]
    print(gene_name)

    change <- extraer_cambio_aa(t)
    print(change)

    my.symbols <- c(gene_name)
    selid <- NA
    selid <- select(hs,
                    keys = my.symbols,
                    columns = c("ENTREZID", "SYMBOL"),
                    keytype = "SYMBOL")

    # for(i in 1:3)
    # {
    #   if(i==2) try(print(a),silent=TRUE)
    #   else print(i)
    # }

    keys <- selid$ENTREZID
    columns <- c("PDB","HGNC","GENES","SEQUENCE")
    kt <- "ENTREZ_GENE"
    res <- NA
    res <- tryCatch(select(up, keys, columns, kt),
                    error=function(e) NULL)
    # res <- select(up, keys, columns, kt)
    if(!is.null(res))
    {
      res <- res[!is.na(as.character(res$SEQUENCE)),]



      for(g in 1:nrow(res))
      {
        print(paste0(g, "/", nrow(res)))
        seq <- res$SEQUENCE[g]
        for(s in 1:length(change))
        {
          cc <- comprobar_estructura(change[[s]])
          if(cc == 1)
          {
            pos <- as.numeric(change[[s]][2])
            sust <- substr(seq, (pos-1), (pos+1))
            if(substr(sust,2,2) == change[[s]][1])
            {

              our_seq <- seq
              our_pos <- pos
              dat$pos[a] <- our_pos
              our_len <- nchar(our_seq)
              dat$seq_len[a] <- our_len
              mut_seq <- NA
              mut_seq <- paste0(substr(our_seq, 1, (pos-1)),change[[s]][3],substr(our_seq, (pos+1), our_len))

              print("done")
              break()
            }
          }
          if(!is.na(our_seq)) break()
        }
        if(!is.na(our_seq)) break()
      }

      # 25 (12 x 2 + 1)
      km <- 8:14
      if(!is.na(our_seq))
      {
        for(e in km)
        {
          if(our_pos < e)
          {
            pos_initial <- 1
          } else {
            pos_initial <- our_pos - e + 1
          }
          if(our_len < (our_pos + e))
          {
            pos_fin <- our_len
          } else {
            pos_fin <- our_pos + e - 1
          }

          final_seq <- NA
          final_seq <- substr(our_seq, pos_initial, pos_fin)

          final_seq_mut <- substr(mut_seq, pos_initial, pos_fin)

          dat[a, paste0("seq_km_",e)] <- ifelse(is.na(final_seq), "",final_seq)
          dat[a, paste0("seq_km_mut",e)] <- ifelse(is.na(final_seq_mut), "",final_seq_mut)
        }
      }
    }


  } # loop mut


  # detectar mutaciones cercanas
  sp <- strsplit(dat$variant_key, ":")
  dat$chr <-  unlist(lapply(sp, `[[`, 1))
  sp2 <-  unlist(lapply(sp, `[[`, 2))
  sp3 <- strsplit(sp2, " ")
  dat$position <-  unlist(lapply(sp3, `[[`, 1))
  dat$verificarcercanos <- 0

  for(q in 1:(nrow(dat)-1))
  {
    print(paste0("Mutacion ", q, " de ", nrow(dat)))
    for(y in (q+1):nrow(dat))
    {
      if(dat$chr[q] == dat$chr[y] & dat$position[q] == dat$position[y] )
      {
        dat$verificarcercanos[q] <- dat$variant_key[y]
      }
    }
  }


  # dat$len_seq <- nchar(dat$seq)
  openxlsx::write.xlsx(dat, file = paste0(zein_karpetara,"/Nas_Pac",pac_num,"_",gaur,"_all_var_call_con_seq.xlsx"))
  return(dat)
}


#' Create the information for FASTA file
#'
#' @description For each gene we will create a FASTA file to run in the mhcpan.
#'
#' @param dat A data frame with an specific structure with needed sequences and gene name
#' @param km The size of the sequence that we want to run
#' @param hla HLA file with the specific information about patients
#' @return Several FASTA files
#' @export
#'
#' @examples
#' \dontrun{
#' create_FASTA_4mhcpan (dat, kmer)
#' }
#'

create_FASTA_4mhcpan <- function(dat, kmer, hla)
{
  print("------------------------------------------------Crear FAST--------------------------------")

  kmer <- 8:11

  for(kk in kmer)
  {
    t <- 1
    print(kk)
    var1 <- paste0("> BORRAR ",1)
    var2 <- "SEQ BORRAR"

    namefile <- paste0(zein_karpetara,"/Nas_Pac",pac_num,"_",gaur,"_all_var_call_con_seq_KM_",kk,".fsa")

    write.table(var1, file = namefile, row.names = F, quote = F, col.names = F)
    write.table(var2, file = namefile, append = T, row.names = F, quote = F, col.names = F)


    for(t in 1:nrow(dat))
    {
      write.table(paste0(">WIL_K",kk,"_",dat$GenName[t]), file = namefile, append = T, row.names = F, quote = F, col.names = F)
      write.table(paste0(dat[t,paste0("seq_km_",kk)]), file = namefile, append = T, row.names = F, quote = F, col.names = F)
      write.table(paste0(">MUT_K",kk,"_",dat$GenName[t]), file = namefile, append = T, row.names = F, quote = F, col.names = F)
      write.table(paste0(dat[t,paste0("seq_km_mut",kk)]), file = namefile, append = T, row.names = F, quote = F, col.names = F)
    }

  }


  # Sacar panMHC

  # obtener los HLA

  hla <- openxlsx::read.xlsx("G:/Mi unidad/NASIR/data/HLA_Nasir_v3_combinacionescorrectas.xlsx")

  hla_A1_pac <- paste0("HLA-A",substr(hla$A[hla$Paciente == paste0("HepaN",pac_num,"_normal_1")],1,5))
  hla_A2_pac <- paste0("HLA-A",substr(hla$A[hla$Paciente == paste0("HepaN",pac_num,"_normal_2")],1,5))

  hla_B1_pac <- paste0("HLA-B",substr(hla$B[hla$Paciente == paste0("HepaN",pac_num,"_normal_1")],1,5))
  hla_B2_pac <- paste0("HLA-B",substr(hla$B[hla$Paciente == paste0("HepaN",pac_num,"_normal_2")],1,5))

  hla_C1_pac <- paste0("HLA-C",substr(hla$C[hla$Paciente == paste0("HepaN",pac_num,"_normal_1")],1,5))
  hla_C2_pac <- paste0("HLA-C",substr(hla$C[hla$Paciente == paste0("HepaN",pac_num,"_normal_2")],1,5))

  hla_all <- unique(c(hla_A1_pac, hla_A2_pac, hla_B1_pac, hla_B2_pac, hla_C1_pac, hla_C2_pac))

  paste0(hla_all, collapse = ",")


  konsolfile <- paste0(zein_karpetara,"/Nas_Pac",pac_num,"_",gaur,"_all_var_call_con_seq_KM_ALL.txt")
  write.table("prog=/home/itamayou/Programas/mhcpan/netMHCpan-4.0/netMHCpan", file = konsolfile, row.names = F, quote = F, col.names = F)
  write.table("data=/home/itamayou/Programas/mhcpan/netMHCpan-4.0/test", append = T, file = konsolfile, row.names = F, quote = F, col.names = F)

  for(kkk in kmer)
  {
    konsola <- paste0("$prog -BA -l ", kkk, " -a ",paste0(hla_all, collapse = ",")," -f ${data}/Nas_Pac",pac_num,"_",gaur,"_all_var_call_con_seq_KM_",kkk,".fsa > ${data}/Nas_Pac",pac_num,"_",gaur,"_all_var_call_con_seq_KM_",kkk,".myout")
    print(konsola)
    write.table(konsola, file = konsolfile, append = T, row.names = F, quote = F, col.names = F)
  }


  print("Copiar los ficheros al cluster")


# CLASE 2 -----------------------------------------------------------------

  hla_DRB1_1_pac <- paste0("DRB1_",substr(hla$DRB1[hla$Paciente == paste0("HepaN",pac_num,"_normal_1")],1,5))
  hla_DRB1_1_pac <- gsub(":","",hla_DRB1_1_pac)
  hla_DRB1_2_pac <- paste0("DRB1_",substr(hla$DRB1[hla$Paciente == paste0("HepaN",pac_num,"_normal_2")],1,5))
  hla_DRB1_2_pac <- gsub(":","",hla_DRB1_2_pac)

  hla_DRB3_1_pac <- paste0("DRB3_",substr(hla$DRB3[hla$Paciente == paste0("HepaN",pac_num,"_normal_1")],1,5))
  hla_DRB3_1_pac <- gsub(":","",hla_DRB3_1_pac)
  hla_DRB3_2_pac <- paste0("DRB3_",substr(hla$DRB3[hla$Paciente == paste0("HepaN",pac_num,"_normal_2")],1,5))
  hla_DRB3_2_pac <- gsub(":","",hla_DRB3_2_pac)

  hla_DRB4_1_pac <- paste0("DRB4_",substr(hla$DRB4[hla$Paciente == paste0("HepaN",pac_num,"_normal_1")],1,5))
  hla_DRB4_1_pac <- gsub(":","",hla_DRB4_1_pac)
  hla_DRB4_2_pac <- paste0("DRB4_",substr(hla$DRB4[hla$Paciente == paste0("HepaN",pac_num,"_normal_2")],1,5))
  hla_DRB4_2_pac <- gsub(":","",hla_DRB4_2_pac)

  hla_DRB5_1_pac <- paste0("DRB5_",substr(hla$DRB5[hla$Paciente == paste0("HepaN",pac_num,"_normal_1")],1,5))
  hla_DRB5_1_pac <- gsub(":","",hla_DRB5_1_pac)
  hla_DRB5_2_pac <- paste0("DRB5_",substr(hla$DRB5[hla$Paciente == paste0("HepaN",pac_num,"_normal_2")],1,5))
  hla_DRB5_2_pac <- gsub(":","",hla_DRB5_2_pac)

  ## DQA

  hla_DQA1_1_pac <- paste0("HLA-DQA1",substr(hla$DQA1[hla$Paciente == paste0("HepaN",pac_num,"_normal_1")],1,5))
  hla_DQA1_1_pac <- gsub(":","",hla_DQA1_1_pac)
  hla_DQA1_2_pac <- paste0("HLA-DQA1",substr(hla$DQA1[hla$Paciente == paste0("HepaN",pac_num,"_normal_2")],1,5))
  hla_DQA1_2_pac <- gsub(":","",hla_DQA1_2_pac)

  # DQB

  hla_DQB1_1_pac <- paste0("DQB1",substr(hla$DQB1[hla$Paciente == paste0("HepaN",pac_num,"_normal_1")],1,5))
  hla_DQB1_1_pac <- gsub(":","",hla_DQB1_1_pac)
  hla_DQB1_2_pac <- paste0("DQB1",substr(hla$DQB1[hla$Paciente == paste0("HepaN",pac_num,"_normal_2")],1,5))
  hla_DQB1_2_pac <- gsub(":","",hla_DQB1_2_pac)

  # combinar DPA + DPB
  kom1 <- paste0(hla_DQA1_1_pac,"-",hla_DQB1_1_pac)
  kom2 <- paste0(hla_DQA1_1_pac,"-",hla_DQB1_2_pac)
  kom3 <- paste0(hla_DQA1_2_pac,"-",hla_DQB1_1_pac)
  kom4 <- paste0(hla_DQA1_2_pac,"-",hla_DQB1_2_pac)

  # DPA

  hla_DPA1_1_pac <- paste0("HLA-DPA1",substr(hla$DPA1[hla$Paciente == paste0("HepaN",pac_num,"_normal_1")],1,5))
  hla_DPA1_1_pac <- gsub(":","",hla_DPA1_1_pac)
  hla_DPA1_2_pac <- paste0("HLA-DPA1",substr(hla$DPA1[hla$Paciente == paste0("HepaN",pac_num,"_normal_2")],1,5))
  hla_DPA1_2_pac <- gsub(":","",hla_DPA1_2_pac)

  # DPB

  hla_DPB1_1_pac <- paste0("DPB1",substr(hla$DPB1[hla$Paciente == paste0("HepaN",pac_num,"_normal_1")],1,5))
  hla_DPB1_1_pac <- gsub(":","",hla_DPB1_1_pac)
  hla_DPB1_2_pac <- paste0("DPB1",substr(hla$DPB1[hla$Paciente == paste0("HepaN",pac_num,"_normal_2")],1,5))
  hla_DPB1_2_pac <- gsub(":","",hla_DPB1_2_pac)

  # combinar DPA + DPB
  com1 <- paste0(hla_DPA1_1_pac,"-",hla_DPB1_1_pac)
  com2 <- paste0(hla_DPA1_1_pac,"-",hla_DPB1_2_pac)
  com3 <- paste0(hla_DPA1_2_pac,"-",hla_DPB1_1_pac)
  com4 <- paste0(hla_DPA1_2_pac,"-",hla_DPB1_2_pac)


  alel <- openxlsx::read.xlsx("G:/Mi unidad/NASIR/data/alelos_hla_clase2.xlsx")

  alel_col <- NULL

  if(hla_DRB1_1_pac %in% alel$alelos) alel_col <- c(alel_col, hla_DRB1_1_pac)
  if(hla_DRB1_2_pac %in% alel$alelos) alel_col <- c(alel_col, hla_DRB1_2_pac)
  if(hla_DRB3_1_pac %in% alel$alelos) alel_col <- c(alel_col, hla_DRB3_1_pac)
  if(hla_DRB3_2_pac %in% alel$alelos) alel_col <- c(alel_col, hla_DRB3_2_pac)
  if(hla_DRB4_1_pac %in% alel$alelos) alel_col <- c(alel_col, hla_DRB4_1_pac)
  if(hla_DRB4_2_pac %in% alel$alelos) alel_col <- c(alel_col, hla_DRB4_2_pac)
  if(hla_DRB5_1_pac %in% alel$alelos) alel_col <- c(alel_col, hla_DRB5_1_pac)
  if(hla_DRB5_2_pac %in% alel$alelos) alel_col <- c(alel_col, hla_DRB5_2_pac)

  if(kom1 %in% alel$alelos) alel_col <- c(alel_col, kom1)
  if(kom2 %in% alel$alelos) alel_col <- c(alel_col, kom2)
  if(kom3 %in% alel$alelos) alel_col <- c(alel_col, kom3)
  if(kom4 %in% alel$alelos) alel_col <- c(alel_col, kom4)

  if(com1 %in% alel$alelos) alel_col <- c(alel_col, com1)
  if(com2 %in% alel$alelos) alel_col <- c(alel_col, com2)
  if(com3 %in% alel$alelos) alel_col <- c(alel_col, com3)
  if(com4 %in% alel$alelos) alel_col <- c(alel_col, com4)

  hla_all <- unique(c(hla_DRB1_1_pac, hla_DRB1_2_pac, hla_DRB3_1_pac, hla_DRB3_2_pac, hla_DRB4_1_pac, hla_DRB4_2_pac,
                      hla_DRB5_1_pac, hla_DRB5_2_pac, alel_col))

  hla_all <- hla_all[grep("_NA",hla_all, invert = TRUE)]

  hla_all <- hla_all[hla_all %in% alel$alelos]

  paste0(hla_all, collapse = ",")


  konsolfile <- paste0(mhcfolder,"/Nas_Pac",pac_num,"_",gaur,"_all_var_call_con_seq_KM_ALL.txt")
  write.table("prog=/home/itamayou/Programas/mhcpan/netMHCIIpan-4.0/netMHCIIpan", file = konsolfile, row.names = F, quote = F, col.names = F)
  write.table("data=/home/itamayou/Programas/mhcpan/netMHCIIpan-4.0/test", append = T, file = konsolfile, row.names = F, quote = F, col.names = F)

  for(kkk in kmer)
  {
    konsola <- paste0("$prog -BA -l ", kkk, " -a ",paste0(hla_all, collapse = ",")," -f ${data}/Nas_Pac",pac_num,"_",gaur,"_all_var_call_con_seq_KM_",kkk,".fsa > ${data}/Nas_Pac",pac_num,"_",gaur,"_all_var_call_con_seq_KM_",kkk,".myout")
    print(konsola)
    write.table(konsola, file = konsolfile, append = T, row.names = F, quote = F, col.names = F)
  }

  print("Copiar los ficheros al cluster")

}

#' Extract the information from mhcpan format
#'
#' @description Once mhcpan is run, we will collect all the information in a user-friendly tables
#'
#' @param folder The foler where we can find the .out files
#' @param kmer The size of the sequence that we want to extract
#' @param hla HLA file with the specific information about patients
#' @return Several FASTA files
#' @export
#'
#' @examples
#' \dontrun{
#' create_FASTA_4mhcpan (dat, kmer)
#' }
#'


from_mhcpan_2table <- function(folder, kmer)
{
  # After MHC ---------------------------------------------------------------

  # kkk <- 8
  kmer <- 8:11

  for(s in kmer)
  {
    did <- NULL
    did <- readLines(paste0(zein_karpetara,"/Nas_Pac",pac_num,"_",gaur,"_all_var_call_con_seq_KM_",s,".myout"))
    namefile2 <- paste0(zein_karpetara,"/Nas_Pac",pac_num,"_",gaur,"temp",s,".txt")
    print(namefile2)
    azken_Excela <- paste0(zein_karpetara,"/Nas_Pac",pac_num,"_",gaur,"excel",s,".xlsx")
    print(azken_Excela)

    rm_lines <- grep("--",did)
    did <- did[-rm_lines]
    rm_lines <- grep("#",did)
    did <- did[-rm_lines]
    rm_lines <- grep("Distance to training data",did)
    did <- did[-rm_lines]
    rm_lines <- grep("Number of high binders",did)
    did <- did[-rm_lines]
    rm_lines <- grep("No peptides",did)
    did <- did[-rm_lines]
    rm_lines <- grep("Score Aff",did)
    did <- did[-rm_lines]
    rm_lines <- which(nchar(did) == 0)
    did <- did[-rm_lines]
    did <- gsub(" <= WB","",did)
    did <- gsub(" <= SB","",did)

    write.table(did, file = namefile2, row.names = F, quote = F, col.names = F)

    df_t <- read.table(namefile2)

    col_names <- c("Pos","HLA","Peptide","Core", "Of", "Gp", "Gl", "Ip", "Il","Icore","Identity","Score","Aff(nM)","Pr_Rank")
    names(df_t) <- col_names

    openxlsx::write.xlsx(df_t, file = azken_Excela)

    df_t2 <- df_t[df_t$Pr_Rank <1, ]
  }


  # Denak elkartu -----------------------------------------------------------


  kmer <- 9:11

  s <- 8
  azken_Excela <- paste0(zein_karpetara,"/Nas_Pac",pac_num,"_",gaur,"excel",s,".xlsx")
  temp1 <- openxlsx::read.xlsx(azken_Excela)

  for(s in kmer)
  {
    azken_Excela <- paste0(zein_karpetara,"/Nas_Pac",pac_num,"_",gaur,"excel",s,".xlsx")
    temp2 <- openxlsx::read.xlsx(azken_Excela)
    temp1 <- rbind(temp1, temp2)
  }

  # Sortu zutabe berriak

  temp3 <- temp1[,c("Pos","HLA","Peptide","Identity","Aff(nM)","Pr_Rank")]

  sp <- strsplit(temp3$Identity, "_")
  temp3$Wild_Mut <-  unlist(lapply(sp, `[[`, 1))
  temp3$kmer <-  unlist(lapply(sp, `[[`, 2))
  temp3$GeneName <-  unlist(lapply(sp, `[[`, 3))

  temp <- temp3

  # selecting

  selected <- temp3[(temp3$Pr_Rank < 2) & temp3$Wild_Mut == "MUT",]
  selected$wild_aff <- 0
  selected$wild_rank <- 0
  selected$seq_mut <- NA

  for(s in 1:nrow(selected))
  {
    t_a <- temp3[temp3$HLA == selected$HLA[s] & temp3$kmer == selected$kmer[s] & temp3$GeneName == selected$GeneName[s]
                 & temp3$Pos == selected$Pos[s] & temp3$Wild_Mut == "WIL",]
    selected$wild_aff[s] <- t_a$`Aff(nM)`[1]
    selected$wild_rank[s] <- t_a$Pr_Rank[1]
    selected$seq_mut[s] <- t_a$Peptide[1]

  }

  # selected$dif_aff <- selected$`Aff(nM)` - selected$wild_aff
  # selected$dif_rank <- selected$Pr_Rank - selected$wild_rank

  selected$DAI <- selected$wild_aff / selected$`Aff(nM)`
  selected$DAI_rank <- selected$wild_rank / selected$Pr_Rank

  selected <- selected[,-1]
  # selected <- unique(selected)
  selected <- selected[order(-selected$DAI),]
  selected <- unique(selected)
  selected <- selected[selected$wild_aff > 500, ]



  openxlsx::write.xlsx(selected, file = paste0(zein_karpetara,"/Nas_Pac",pac_num,"_",gaur,"excel",s,"_selected.xlsx"))


}
