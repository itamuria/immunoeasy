# Read all files and merge them

directory <- "/home/itamayou/SaraLabiano/005_Counts"
file_name <- "20220704_sara_counts.xlsx"
remove_part_name <- ".counts.txt"
files_names <- dir(directory)
setwd(directory)

dat <- read.table(files_names[1], header = T)
dat = dat[,c(1,7)]
names(dat) <- c("ENSG",gsub(remove_part_name, "", files_names[1]))

for(e in 2:length(files_names))
{
  print(paste0(e, ": ", files_names[e]))
  dat_temp <- read.table(files_names[e], header = T)
  dat_temp = dat_temp[,c(1,7)]
  names(dat_temp) <- c("ENSG",gsub(remove_part_name, "", files_names[e]))

  dat <- merge(dat, dat_temp, by = "ENSG")
}

openxlsx::write.xlsx(dat, file = file_name)
