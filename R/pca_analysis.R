library(factoextra)
res3 <- scale(res2[,-1],scale = TRUE)
res.pca <- prcomp(res3[,-1], scale = TRUE)

fviz_eig(res.pca)

pat <- openxlsx::read.xlsx("G:/Mi unidad/NASIR/MAIN_FOLDER_NASIR/MAIN_TABLE_v002.xlsx",2)

pat2 <- pat$Progresion[pat$TUMOR_RNA %in% res2$Patients]
pat3 <- res2$Mast.cells.activated
pat3 <- res2$Mast.cells.resting

fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_ind(res.pca,
             col.ind = pat2, # Color by the quality of representation
             repel = F     # Avoid text overlapping
)
fviz_pca_ind(res.pca,
             col.ind = pat3, # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
)

fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)


library(factoextra)
# Eigenvalues
eig.val <- get_eigenvalue(res.pca)
eig.val

# Results for Variables
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation
# Results for individuals
res.ind <- get_pca_ind(res.pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation

library(ggplot2)
# plotCIBER <- function(ciber,nrow=NULL) {
#   ggplot(aes(x = Var1, y = value, fill=Var1), data=ciber) + geom_bar(stat="summary", fun.y=mean) + scale_x_discrete(labels = annos) + scale_fill_manual(values=colors, labels=annos, guide=F) + facet_wrap(~experiment,nrow=nrow) + geom_errorbar(stat="summary",fun.ymin = function(x) mean(x)-sd(x), fun.ymax = function(x) mean(x)+sd(x), width=0.2) + theme(axis.text.x =  element_text(angle=90, size=8)) + ylab("% of population") + xlab("")
# }
#
# plotCIBER(results)
#

library(tidybulk)
# https://github.com/stemangiola/tidybulk
