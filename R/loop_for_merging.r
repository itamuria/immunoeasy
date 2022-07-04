# Read all files and merge them

directory <- "/home/itamayou/NASIR_RNA/ConteosTotales/uniendo/htseq_no"
file_name <- "20210602_htseq_no_nasir.xlsx"
remove_part_name <- "_htseq_counts_no_3"
files_names <- dir(directory)

dat <- read.table(files_names[1])
names(dat) <- c("ENSG",gsub(remove_part_name, "", files_names[1]))

for(e in 2:length(files_names))
{
  print(paste0(e, ": ", files_names[e]))
  dat_temp <- read.table(files_names[e])
  names(dat_temp) <- c("ENSG",gsub(remove_part_name, "", files_names[e]))

  dat <- merge(dat, dat_temp, by = "ENSG")
}

openxlsx::write.xlsx(dat, file = file_name)