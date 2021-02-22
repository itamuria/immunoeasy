#' Filtering mutations with certain criteria
#'
#' @description Mutation list will be filtered based on some criteria
#'
#' @param dataframe The data frame with all the information
#' @param variant_caller Select the variant caller that generated the data frame.
#' @param PrimaryFilter If we want to apply PASS filter. Boolean with True by default.
#' @param FuncRefGene What function of the Reference Gene we want to apply. Default is exonic.
#' @param ExonicFunc What exonic function of the Reference Gene we want to apply. Default value is nonsynonymous_SNV
#' @param Cob_DNA_tum Minimum deep of the tumoral DNA (bigger values  will be considered).
#' @param VAF_tum Minimum alternative frecuency of the tumor (bigger values  will be considered).
#' @param AD_alt_tum Minimum deep of alternative allele of the tumor (bigger values  will be considered).
#' @param AD_alt Maximum deep of alternative allele of the healthy tissue (smaller values and selected value will be considered).
#'
#' @return data frame filtered based on the selected criteria
#' @export
#'
#' @examples
#' \dontrun{
#' filtering_mutations (filename)
#' }
#'
filtering_mutations <- function(dataframe, variant_caller = "Mutect", PrimaryFilter = T, FuncRefGene = "exonic", ExonicFunc = "nonsynonymous_SNV",
                                Cob_DNA_tum = 6, VAF_tum = 0.05, AD_alt_tum = 3, AD_alt = 1)
{

    if(variant_caller == "varscan")
    {
        names(dataframe)[which(names(dataframe) == "Tum_DP")] <- "Cob_DNA_tum"
        names(dataframe)[which(names(dataframe) == "Tum_VAF")] <- "AF_tum"
        names(dataframe)[which(names(dataframe) == "Tum_Alt_D")] <- "AD_alt_tum"
        names(dataframe)[which(names(dataframe) == "Nor_Alt_D")] <- "AD_alt"
    }
    if(variant_caller == "strelka")
    {
        names(dataframe)[which(names(dataframe) == "DP_2")] <- "Cob_DNA_tum"
        names(dataframe)[which(names(dataframe) == "VAF_T1")] <- "AF_tum"
        names(dataframe)[which(names(dataframe) == "Alt_T1_tumor")] <- "AD_alt_tum"
        names(dataframe)[which(names(dataframe) == "Alt_T1_normal")] <- "AD_alt"
    }
    if(variant_caller == "SomSni")
    {
        names(dataframe)[which(names(dataframe) == "DP_2")] <- "Cob_DNA_tum"
        names(dataframe)[which(names(dataframe) == "vaf2")] <- "AF_tum"
        names(dataframe)[which(names(dataframe) == "alt2")] <- "AD_alt_tum"
        dataframe$AD_alt <- dataframe$Norm_DP_alt_forw + dataframe$Norm_DP_alt_rev
        dataframe$Filt <- "PASS"
    }
    sel0_varcall <- dataframe[dataframe$Filt == "PASS",]
    sel1_varcall <- dataframe[dataframe$Filt == "PASS" & dataframe$Func.refGene == FuncRefGene,]
    sel2_varcall <- dataframe[dataframe$Filt == "PASS" & dataframe$Func.refGene == FuncRefGene & dataframe$ExonicFunc.refGene == ExonicFunc,]
    sel3_varcall <- dataframe[dataframe$Filt == "PASS" & dataframe$Func.refGene == FuncRefGene & dataframe$ExonicFunc.refGene == ExonicFunc &
                               dataframe$Cob_DNA_tum > Cob_DNA_tum ,]
    sel4_varcall <- dataframe[dataframe$Filt == "PASS" & dataframe$Func.refGene == FuncRefGene & dataframe$ExonicFunc.refGene == ExonicFunc &
                               dataframe$Cob_DNA_tum > Cob_DNA_tum & dataframe$AF_tum > VAF_tum,]
    sel5_varcall <- dataframe[dataframe$Filt == "PASS" & dataframe$Func.refGene == FuncRefGene & dataframe$ExonicFunc.refGene == ExonicFunc &
                               dataframe$Cob_DNA_tum > Cob_DNA_tum & dataframe$AF_tum > VAF_tum & dataframe$AD_alt_tum > AD_alt_tum,]

    sel6_varcall <- dataframe[dataframe$Filt == "PASS" & dataframe$Func.refGene == FuncRefGene & dataframe$ExonicFunc.refGene == ExonicFunc &
                               dataframe$Cob_DNA_tum > Cob_DNA_tum & dataframe$AF_tum > VAF_tum & dataframe$AD_alt_tum > AD_alt_tum & dataframe$AD_alt <= AD_alt,]

    izenburuak <- c("VarCall","PASS","Exonic","Nonsynonymous","DNACobTum", "AF_Tum","DeepAltTum","Alt_Deep_Normal")
    kop_varcall <- c(nrow(dataframe), nrow(sel0_varcall), nrow(sel1_varcall),nrow(sel2_varcall), nrow(sel3_varcall), nrow(sel4_varcall),nrow(sel5_varcall),nrow(sel6_varcall))


    varcall_sel <- sel6_varcall[,c("Chr","Pos", "Ref","Alt", "Gene.refGene", "AAChange.refGene")]
    varcall_sel$Chr <- gsub("chr","",varcall_sel$Chr)
    varcall_sel$variant_key <- paste0(varcall_sel$Chr, ":", varcall_sel$Pos, " ", varcall_sel$Ref, " > ", varcall_sel$Alt)
    varcall_sel$Mutated_Protein <- ""
    varcall_sel$Radia <- "No"
    varcall_sel$Mutect38 <- "No"
    varcall_sel$Varscan2 <- "No"
    varcall_sel$Strelka <- "No"
    varcall_sel$SomaticSniper <- "No"
    varcall_sel$HowManyVC <- NA

    varcall_sel$Nucleotide_sequence <- NA
    varcall_sel$Long_Peptide_Screened <- NA

    varcall_sel$RNA_varscan2 <- NA
    varcall_sel$Cuartil <- NA
    varcall_sel$FPKM <- NA



    varcall_sel2 <- varcall_sel[,c("variant_key","Gene.refGene","Mutated_Protein","AAChange.refGene","Radia","Mutect38", "Varscan2","Strelka","SomaticSniper",
                               "HowManyVC","Nucleotide_sequence","Long_Peptide_Screened","RNA_varscan2","Cuartil","FPKM")]
    names(varcall_sel2)[c(2,4)] <- c("Gene_name","cDNA_Protein_change")


    return(list(varcall_sel2 = varcall_sel2, kop_varcall = kop_varcall))

}

#' Comparing similarities between variant callers
#'
#' @description Mutation list will be filtered based on some criteria
#'
#' @param dataframe The data frame with all the information
#' @param variant_caller Select the variant caller that generated the data frame.
#' @param PrimaryFilter If we want to apply PASS filter. Boolean with True by default.
#' @param FuncRefGene What function of the Reference Gene we want to apply. Default is exonic.
#' @param ExonicFunc What exonic function of the Reference Gene we want to apply. Default value is nonsynonymous_SNV
#' @param Cob_DNA_tum Minimum deep of the tumoral DNA (bigger values  will be considered).
#' @param VAF_tum Minimum alternative frecuency of the tumor (bigger values  will be considered).
#' @param AD_alt_tum Minimum deep of alternative allele of the tumor (bigger values  will be considered).
#' @param AD_alt Maximum deep of alternative allele of the healthy tissue (smaller values and selected value will be considered).
#'
#' @return data frame filtered based on the selected criteria
#' @export
#'
#' @examples
#' \dontrun{
#' comparing_var_call (filename)
#' }
#'
comparing_var_call <- function(data_filt_mut38, data_filt_varscan,
                               data_filt_strelka, data_filt_SomSni,
                               group = c("variant_key", "Gene_name"),
                               results_folder,
                               data_summ_mut38, data_summ_SomSni,
                               data_summ_strelka, data_summ_varscan)
{
    # Comparing

    key_all <- c(data_filt_mut38[,group],data_filt_varscan[,group],
                 data_filt_strelka[,group], data_filt_SomSni[,group])
    sort(table(key_all))
    un_key <- unique(key_all)

    # merge data frames

    df_all <- rbind(data_filt_mut38,data_filt_varscan)
    df_all2 <- rbind(data_filt_strelka,data_filt_SomSni)
    df_all_fin <- unique(rbind(df_all,df_all2))
    head(df_all_fin)

    # Fill Positions

    for(r in 1:length(un_key))
    {
        if(un_key[r]%in%data_filt_mut38[ , group ])
        {
            df_all_fin[df_all_fin[ , group ] == un_key[r],"Mutect38"] <- "YES"
        }
        if(un_key[r]%in%data_filt_varscan[ , group ])
        {
            df_all_fin[df_all_fin[ , group ] == un_key[r],"Varscan2"] <- "YES"
        }
        if(un_key[r]%in%data_filt_strelka[ , group ])
        {
            df_all_fin[df_all_fin[ , group ] == un_key[r],"Strelka"] <- "YES"
        }
        if(un_key[r]%in%data_filt_SomSni[ , group ])
        {
            df_all_fin[df_all_fin[ , group ] == un_key[r],"SomaticSniper"] <- "YES"
        }
    }

    df_all_fin$HowManyVC <- (ifelse(df_all_fin$Mutect38 == "YES", 1, 0)) +
        (ifelse(df_all_fin$Varscan2 == "YES", 1, 0)) +
        (ifelse(df_all_fin$Strelka == "YES", 1, 0)) +
        (ifelse(df_all_fin$SomaticSniper == "YES", 1, 0))

    subresults <- paste0(results_folder, "/", group)
    dir.create(subresults)

    write.xlsx(df_all_fin, file=paste0(subresults,"/",gaur,"_df_with_",group,".xlsx"))

    save(df_all_fin, file=paste0(subresults,"/",gaur,"_df_with_",group,".RData"))

    df_how <- as.data.frame(table(df_all_fin$HowManyVC))
    names(df_how)[1] <- "HowMany"
    ggplot(df_how, aes(x=HowMany, y=Freq, fill=HowMany)) +
        geom_bar(stat="identity") +
        geom_text(aes(label=Freq), position=position_dodge(width=0.9), size=4)  +
        coord_flip()

    ggsave(paste0(subresults,"/",gaur,"_",group,"_HowManyVariantCallers.jpg"),width = 4.73,
           height = 4.62)

    df_how2 <- as.data.frame(table(df_all_fin[,group]))
    df_how2 <- df_how2[df_how2$Freq > 1,]
    names(df_how2)[1] <- "HowMany"
    df_how2$HowMany <- as.character(df_how2$HowMany)
    df_how2 <- df_how2[order(-df_how2$Freq),]
    df_how2$HowMany <- factor(df_how2$HowMany, levels = df_how2$HowMany)


    ggplot(df_how2, aes(x=HowMany, y=Freq, fill=HowMany)) +
        geom_bar(stat="identity") +
        geom_text(aes(label=Freq), position=position_dodge(width=0.9), size=4)  +
        coord_flip() + theme(legend.position = "none")

    ggsave(paste0(subresults,"/",gaur,"_",group,"_Repeated_",group,".jpg"),
           width = 4.73,
           height = 4.62)


    grep("data_summ",ls(),value = T)

    # Mutect38 <- kop_mut38
    # SomanitSniper <- kop_somsni
    # Strelka <- kop_strelka
    # VarScan2 <- kop_varscan_snv

    izenburuak <- c("VarCall","PASS","Exonic","Nonsynonymous","DNACobTum", "AF_Tum","DeepAltTum","Alt_Deep_Normal")

    df_pyramid <- data.frame(izenburuak,data_summ_mut38,data_summ_SomSni,data_summ_strelka,data_summ_varscan)

    df_pyramid2 <- melt(df_pyramid,id=c("izenburuak"))

    df_pyramid2$izenburuak <- factor(df_pyramid2$izenburuak, levels = c("Alt_Deep_Normal","DeepAltTum","AF_Tum","DNACobTum","Nonsynonymous","Exonic","PASS","VarCall"))

    ggplot(df_pyramid2,aes(izenburuak, value, fill = variable)) +
        geom_col(show.legend = FALSE) +
        labs(x = NULL, y = "n") + geom_text(aes(label=value), position=position_dodge(width=0.9), size=4)  +
        facet_wrap(~variable, ncol = 2, scales = "free") +
        coord_flip()

    ggsave(paste0(subresults,"/",gaur,"_",group,"_Repeated_",group,"_Pyramids.jpg"),width = 4.73,
           height = 4.62)
    # ggsave(paste0(zein_karpetara,"/",gaur,"_all_pac",pac_num,"_Pyramids.pdf"),width = 4.73,
    #        height = 4.62)


    df_pyramid3 <- df_pyramid2[!df_pyramid2$izenburuak%in%c("VarCall","PASS"),]

    ggplot(df_pyramid3,aes(izenburuak, value, fill = variable)) +
        geom_col(show.legend = FALSE) +
        labs(x = NULL, y = "n") + geom_text(aes(label=value), position=position_dodge(width=0.9), size=4)  +
        facet_wrap(~variable, ncol = 2, scales = "free") +
        coord_flip()

    ggsave(paste0(subresults,"/",gaur,"_",group,"_Repeated_",group,"_Pyramids2.jpg"),width = 4.73,
           height = 4.62)
}
