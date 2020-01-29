#' From excel to How Many varian callers
#'
#' @description From an excel with mutation detected in each varian caller, the function convert the data in a table with information about in how many variant callers appeared the mutation.
#'
#' @param filename The name of the excel file with the mutations with choromosome, position, gen name and variant caller
#' @param chr_pos Position of the chromosome name in the data frame
#' @param position Position of the mega position in the data frame
#' @param gen_name Position of the Gene name in the data frame
#' @param varian_caller Position of the Variant Caller in the data frame
#' @param VAF Position of the VAF name in the data frame (if is not present, NA will be presented)
#' @param others Position of the other variables in the data frame that will be merged
#' @param var_cal_4 Names of the four variant caller to calculate the How many field
#'
#'
#' @return data frame with How many variant callers per mutation
#' @export
#'
#' @examples
#' \dontrun{
#' varcall2HowMany (filename)
#' }
#'
varcall2HowMany <- function(filename, chr_pos = 1, position = 2, gen_name = 3, varian_caller = 4, VAF = NA, others = NULL,
                            var_cal_4 = c("mutect38","somaticsniper", "strelka","varscan"))
{

  # Load data  ------------------------------------------------------------
  dat <- openxlsx::read.xlsx(filename)

  dat[,chr_pos] <- as.character(gsub(" ","",dat[,chr_pos]))
  dat[,position] <- as.numeric(as.character(gsub(" ","",dat[,position])))
  dat[,gen_name] <- as.character(gsub(" ","",dat[,gen_name]))

  dat$unique <- paste0(dat[,chr_pos],"_",dat[,position],"_",dat[,gen_name])

  un_unique <- unique(dat$unique)
  un_varcall <- unique(dat[,varian_caller])

  df_out <- dat[1,c(ncol(dat),chr_pos,position,gen_name)]
  df_out <- data.frame(matrix(NA,ncol = (5+length(un_varcall)),nrow = length(un_unique)))
  names(df_out) <- c("unique","Chr","Pos","GenName",un_varcall, "HowMany")

  if(!is.na(VAF))
  {
    df_out$VAFmean <- NA
    df_out$VAFsd <- 0
  }

  if(length(others) > 0)
  {
    for(r in 1:length(others))
    {
      df_out[,names(dat)[others[r]]] <- NA
    }
  }

  h <- 0

  for(q in un_unique)
  {
    temp <- dat[dat$unique == q, ]
    h <- h + 1

    df_out$unique[h] <- q
    df_out$Chr[h] <- temp[1,chr_pos]
    df_out$Pos[h] <- temp[1,position]
    df_out$GenName[h] <- temp[1,gen_name]

    for(v in un_varcall)
    {
      if(v%in%temp[,varian_caller]) df_out[h,v] <- v
    }

    if(length(var_cal_4) == 4) df_out$HowMany[h] <- 4-(is.na(df_out[h,var_cal_4[1]])+is.na(df_out[h,var_cal_4[2]])+is.na(df_out[h,var_cal_4[3]])+is.na(df_out[h,var_cal_4[4]]))
    if(length(var_cal_4) == 3) df_out$HowMany[h] <- 3-(is.na(df_out[h,var_cal_4[1]])+is.na(df_out[h,var_cal_4[2]])+is.na(df_out[h,var_cal_4[3]]))
    if(length(var_cal_4) == 2) df_out$HowMany[h] <- 2-(is.na(df_out[h,var_cal_4[1]])+is.na(df_out[h,var_cal_4[2]]))


    if(!is.na(VAF))
    {
      df_out$VAFmean[h] <- mean(temp[,VAF], na.rm = TRUE)
      df_out$VAFsd[h] <- sd(temp[,VAF], na.rm = TRUE)

    }

    if(length(others) > 0)
    {
      for(r in 1:length(others))
      {
        u1 <- unique(temp[,others[r]])
        u2 <- NULL
        for(p in 1:length(u1))
        {
          u2 <- paste0(u2,";",u1[p])
        }
        u3 <- gsub(";NA","",u2)
        u3 <- gsub("^;","",u3)
        df_out[h,names(dat)[others[r]]] <- u3

      }
    }
  }

  return(df_out)
}
