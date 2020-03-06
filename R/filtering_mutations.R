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

    sel0_mut38 <- data_mut38[data_mut38$Filt == "PASS",]
    sel1_mut38 <- data_mut38[data_mut38$Filt == "PASS" & data_mut38$Func.refGene == "exonic",]
    sel2_mut38 <- data_mut38[data_mut38$Filt == "PASS" & data_mut38$Func.refGene == "exonic" & data_mut38$ExonicFunc.refGene == "nonsynonymous_SNV",]
    sel3_mut38 <- data_mut38[data_mut38$Filt == "PASS" & data_mut38$Func.refGene == "exonic" & data_mut38$ExonicFunc.refGene == "nonsynonymous_SNV" &
                               data_mut38$Cob_DNA_tum > 6 ,]
    sel4_mut38 <- data_mut38[data_mut38$Filt == "PASS" & data_mut38$Func.refGene == "exonic" & data_mut38$ExonicFunc.refGene == "nonsynonymous_SNV" &
                               data_mut38$Cob_DNA_tum > 6 & data_mut38$AF_tum > 0.05,]
    sel5_mut38 <- data_mut38[data_mut38$Filt == "PASS" & data_mut38$Func.refGene == "exonic" & data_mut38$ExonicFunc.refGene == "nonsynonymous_SNV" &
                               data_mut38$Cob_DNA_tum > 6 & data_mut38$AF_tum > 0.05 & data_mut38$AD_alt_tum > 3,]

    sel6_mut38 <- data_mut38[data_mut38$Filt == "PASS" & data_mut38$Func.refGene == "exonic" & data_mut38$ExonicFunc.refGene == "nonsynonymous_SNV" &
                               data_mut38$Cob_DNA_tum > 6 & data_mut38$AF_tum > 0.05 & data_mut38$AD_alt_tum > 3 & data_mut38$AD_alt <= 1,]

    izenburuak <- c("VarCall","PASS","Exonic","Nonsynonymous","DNACobTum", "AF_Tum","DeepAltTum","Alt_Deep_Normal")
    kop_mut38 <- c(nrow(data_mut38), nrow(sel0_mut38), nrow(sel1_mut38),nrow(sel2_mut38), nrow(sel3_mut38), nrow(sel4_mut38),nrow(sel5_mut38),nrow(sel6_mut38))


    mut38_sel <- sel6_mut38[,c("Chr","Pos", "Ref","Alt", "Gene.refGene", "AAChange.refGene")]
    mut38_sel$Chr <- gsub("chr","",mut38_sel$Chr)
    mut38_sel$variant_key <- paste0(mut38_sel$Chr, ":", mut38_sel$Pos, " ", mut38_sel$Ref, " > ", mut38_sel$Alt)
    mut38_sel$Mutated_Protein <- ""
    mut38_sel$Radia <- "No"
    mut38_sel$Mutect38 <- "No"
    mut38_sel$Varscan2 <- "No"
    mut38_sel$Strelka <- "No"
    mut38_sel$SomaticSniper <- "No"
    mut38_sel$HowManyVC <- NA

    mut38_sel$Nucleotide_sequence <- NA
    mut38_sel$Long_Peptide_Screened <- NA

    mut38_sel$RNA_varscan2 <- NA
    mut38_sel$Cuartil <- NA
    mut38_sel$FPKM <- NA



    mut38_sel2 <- mut38_sel[,c("variant_key","Gene.refGene","Mutated_Protein","AAChange.refGene","Radia","Mutect38", "Varscan2","Strelka","SomaticSniper",
                               "HowManyVC","Nucleotide_sequence","Long_Peptide_Screened","RNA_varscan2","Cuartil","FPKM")]
    names(mut38_sel2)[c(2,4)] <- c("Gene_name","cDNA_Protein_change")

    sel6_mut38$id <- 1:nrow(sel6_mut38)
    sel6_mut38$NAF_tum <- sel6_mut38$AD_ref_tum/(sel6_mut38$AD_ref_tum+sel6_mut38$AD_alt_tum)

    mut38_vaf <- melt(sel6_mut38[,c("id","AF_tum","NAF_tum")],id=c("id"))


    return(mut38_sel)

}
