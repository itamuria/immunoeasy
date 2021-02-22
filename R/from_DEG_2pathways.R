#' Search the most enriched pathways
#'
#' @description From DEG we want to obtain the information about the pathways
#'
#' @param res_frame Table with the DEG
#' @param wpgmtfile wpgmtfile from wikipath
#' @return Several excel file with information about enriched pathways
#' @export
#'
#' @examples
#' \dontrun{
#' from_DEG_2pathways (res_frame)
#' }


from_DEG_2pathways <- function(res_frame,wpgmtfile)
{
  #   # Wikipath ------------------------------------------------------------

  library(magrittr)
  library(clusterProfiler)

  # data(geneList, package="DOSE")
  # gene <- names(geneList)[abs(geneList) > 2]

  resframe8 <- res_frame[res_frame$padj < 0.1,]
  resframe8 <- resframe8[!is.na(resframe8$trans),]
  entrezid <- mapIds(org.Hs.eg.db,
                     keys=gsub("\\..*","",rownames(resframe8)),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
  entrezid2 <- unname(entrezid)
  entrezid2 <- entrezid2[!is.na(entrezid2)]



  wpgmtfile <- "wikipathways-20210110-gmt-Homo_sapiens.gmt"
  wp2gene <- read.gmt(wpgmtfile)
  ont <- "term"
  wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
  wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
  wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

  ewp <- enricher(entrezid2, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
  head(ewp)


  # GSEA

  # res_val3 <- res_val3[order(res_val3$log2FoldChange, decreasing = TRUE),]


  resframeGSEA <- res_frame[order(res_frame$log2FoldChange, decreasing = TRUE),]
  resframeGSEA <- resframeGSEA[!is.na(resframeGSEA$trans),]
  dGSEA <- resframeGSEA$trans

  entrezidGSEA <- mapIds(org.Hs.eg.db,
                         keys=gsub("\\..*","",dGSEA),
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")
  entrezid2GSEA <- unname(entrezidGSEA)

  geneListGSEA <- resframeGSEA$log2FoldChange
  names(geneListGSEA) <- entrezid2GSEA

  entrezid2GSEA_sinNa <- !is.na(names(geneListGSEA))
  geneListGSEA_sinNA <- geneListGSEA[entrezid2GSEA_sinNa]
  geneListGSEA_sinNA2 <- geneListGSEA_sinNA[!is.na(geneListGSEA_sinNA)]

  ewp2 <- GSEA(geneListGSEA_sinNA2, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE)
  head(ewp2)

  heatplot(ewp, foldChange=geneListGSEA_sinNA2)
  emapplot(ewp2)
  upsetplot(ewp)

  library(org.Hs.eg.db)
  ewp <- setReadable(ewp, org.Hs.eg.db, keyType = "ENTREZID")
  ewp2 <- setReadable(ewp2, org.Hs.eg.db, keyType = "ENTREZID")
  head(ewp)
  # ridgeplot(ewp)
  gseaplot(ewp2, geneSetID = 1, title = ewp2$Description[1])
  gseaplot2(ewp2, geneSetID = 1, title = ewp2$Description[1])

  gseaplot2(ewp2, geneSetID = 1:3, pvalue_table = TRUE,
            color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")

  openxlsx::write.xlsx(data.frame(ewp), file = paste0("ResultsCounts2/",contraste,"_Wikipath_hipergeometrico.xlsx"))
  openxlsx::write.xlsx(data.frame(ewp2), file = paste0("ResultsCounts2/",contraste,"_Wikipath_gsea.xlsx"))

  # ploting

  library(enrichplot)
  barplot(ewp, showCategory=4)
  dotplot(ewp, showCategory=30) + ggtitle("Wikipath for Progresion")

  dotplot(ewp2, showCategory=30) + ggtitle("Wikipath for Progresion")


  ### Plot heatmap

  df2 <- data.frame(ewp2)
  df2_genes <- df2$core_enrichment[3]

  sel_gen <- unlist(strsplit(df2_genes, "/"))

  normalized_counts2 <- counts(dds, normalized=TRUE)

  sample[,contraste] <- as.character(sample[,contraste])
  TTPlar <- ifelse(is.na(sample[,contraste]), "NA", sample[,contraste])
  group_df <- data.frame(TTPlar[sele2])
  rownames(group_df) <- sample$TUMOR_RNA[sele2]
  names(group_df) <- contraste

  symbol <- mapIds(org.Hs.eg.db,
                   keys=gsub("\\..*","",rownames(normalized_counts2)),
                   column="SYMBOL",
                   keytype="ENSEMBL",
                   multiVals="first")

  symbol <- unname(ifelse(is.na(symbol), names(symbol), symbol))

  heatdf <- normalized_counts2[symbol %in% sel_gen,]

  symbol_df <- data.frame(symbol)
  rownames(symbol_df) <- rownames(heatdf)
  symbol2 <- symbol[symbol %in% sel_gen]

  pheatmap(heatdf, scale = "row", labels_col = sample$ID[sele2], cluster_cols = FALSE,
           annotation_col = group_df, labels_row = symbol2, fontsize_row = 5)


  # # Cell marker -----------------------------------------------------------



  cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
    tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>%
    dplyr::select(cellMarker, geneID) %>%
    dplyr::mutate(geneID = strsplit(geneID, ', '))
  cell_markers

  y <- enricher(entrezid2, TERM2GENE=cell_markers, minGSSize=1)
  DT::datatable(as.data.frame(y))


  openxlsx::write.xlsx(  data.frame(y), file = paste0("ResultsCounts2/",contraste,"_CellMarker.xlsx"))

  # MSigDb


  library(msigdbr)
  msigdbr_show_species()

  m_df <- msigdbr(species = "Homo sapiens")
  head(m_df, 2) %>% as.data.frame

  m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>%
    dplyr::select(gs_name, entrez_gene)
  head(m_t2g)

  em <- enricher(entrezid2, TERM2GENE=m_t2g)
  em2 <- GSEA(geneListGSEA_sinNA2, TERM2GENE = m_t2g)
  head(em)

  library(org.Hs.eg.db)
  em <- setReadable(em, org.Hs.eg.db, keyType = "ENTREZID")
  em2 <- setReadable(em2, org.Hs.eg.db, keyType = "ENTREZID")
  head(em)
  head(em2)

  openxlsx::write.xlsx(data.frame(em), file = paste0("ResultsCounts2/",contraste,"_MSigDb_hipergeometrico.xlsx"))
  openxlsx::write.xlsx(data.frame(em2), file = paste0("ResultsCounts2/",contraste,"_MSigDb_gsea.xlsx"))

  # C-s
  m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>%
    dplyr::select(gs_name, entrez_gene)
  head(m_t2g)

  em3 <- GSEA(geneListGSEA_sinNA2, TERM2GENE = m_t2g)
  em3 <- setReadable(em3, org.Hs.eg.db, keyType = "ENTREZID")
  d9 <- data.frame(em3)
  gseaplot(em3, geneSetID = 7, title = em3$Description[7])

  for(call in c(paste0("C",1:7),"H"))
  {
    print(call)
    m_t2g <- NULL
    m_t2g <- msigdbr(species = "Homo sapiens", category = call) %>%
      dplyr::select(gs_name, entrez_gene)
    em3 <- GSEA(geneListGSEA_sinNA2, TERM2GENE = m_t2g)
    em3 <- setReadable(em3, org.Hs.eg.db, keyType = "ENTREZID")
    # openxlsx::write.xlsx(data.frame(em3), file = paste0("ResultsCounts2/",contraste,"_MSigDb_gsea_",call,".xlsx"))
  }

  # Disease
  library(DOSE)
  x <- enrichDO(gene          = entrezid2,
                ont           = "DO",
                pvalueCutoff  = 0.05,
                pAdjustMethod = "BH",
                universe      = names(geneListGSEA_sinNA2),
                minGSSize     = 5,
                maxGSSize     = 500,
                qvalueCutoff  = 0.05,
                readable      = FALSE)
  head(x)
  x <- setReadable(x, org.Hs.eg.db, keyType = "ENTREZID")

  gene2 <- names(geneListGSEA_sinNA2)[abs(geneListGSEA_sinNA2) < 3]
  ncg <- enrichNCG(gene2)
  head(ncg)

  dgn <- enrichDGN(entrezid2)
  dgn(dgn)
  dgn <- setReadable(dgn, org.Hs.eg.db, keyType = "ENTREZID")
  openxlsx::write.xlsx(data.frame(dgn), file = paste0("ResultsCounts2/",contraste,"_dgn_.xlsx"))

  library(DOSE)

  y <- gseDO(geneListGSEA_sinNA2,
             nPerm         = 100000,
             minGSSize     = 120,
             pvalueCutoff  = 0.2,
             pAdjustMethod = "BH",
             verbose       = FALSE)
  head(y, 3)

  openxlsx::write.xlsx(data.frame(y), file = paste0("ResultsCounts2/",contraste,"_gseDO_gsea.xlsx"))


  ncg1 <- enrichDGN(entrezid2)


  ncg2 <- gseNCG(geneListGSEA_sinNA2,
                 nPerm         = 10000,
                 minGSSize     = 120,
                 pvalueCutoff  = 0.2,
                 pAdjustMethod = "BH",
                 verbose       = FALSE)
  ncg1 <- setReadable(ncg1, 'org.Hs.eg.db')
  ncg2 <- setReadable(ncg2, 'org.Hs.eg.db')
  head(ncg1, 3)
  openxlsx::write.xlsx(data.frame(ncg1), file = paste0("ResultsCounts2/",contraste,"_ncg_enrich.xlsx"))
  openxlsx::write.xlsx(data.frame(ncg2), file = paste0("ResultsCounts2/",contraste,"_ncg_gsea.xlsx"))

  p1 <- dotplot(ncg1, showCategory=30) + ggtitle("dotplot for Enrich")
  p2 <- dotplot(ncg2, showCategory=30) + ggtitle("dotplot for GSEA")
  plot_grid(p1, p2, ncol=2)





  dgn <- gseDGN(geneListGSEA_sinNA2,
                nPerm         = 10000,
                minGSSize     = 120,
                pvalueCutoff  = 0.2,
                pAdjustMethod = "BH",
                verbose       = FALSE)
  dgn <- setReadable(dgn, 'org.Hs.eg.db')
  head(dgn, 3)
  openxlsx::write.xlsx(data.frame(dgn), file = paste0("ResultsCounts2/",contraste,"_dgn_gsea.xlsx"))


  # Gene Ontology Analysis

  library(clusterProfiler)
  data(geneList, package="DOSE")
  gene <- names(geneList)[abs(geneList) > 2]
  gene.df <- bitr(entrezid2, fromType = "ENTREZID",
                  toType = c("ENSEMBL", "SYMBOL"),
                  OrgDb = org.Hs.eg.db)
  head(gene.df)

  ggo <- groupGO(gene     = entrezid2,
                 OrgDb    = org.Hs.eg.db,
                 ont      = "CC",
                 level    = 2,
                 readable = TRUE)

  head(ggo)

  openxlsx::write.xlsx(data.frame(ggo), file = paste0("ResultsCounts2/",contraste,"_ggo.xlsx"))


  ego <- enrichGO(gene          = entrezid2,
                  universe      = names(geneListGSEA_sinNA2),
                  OrgDb         = org.Hs.eg.db,
                  ont           = "CC",
                  pAdjustMethod = "none",
                  readable      = TRUE)
  head(ego)

  ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "CC",
                   pAdjustMethod = "none")

  head(ego2)

  ego <- setReadable(ego, OrgDb = org.Hs.eg.db)
  ego2 <- setReadable(ego2, OrgDb = org.Hs.eg.db)

  openxlsx::write.xlsx(data.frame(ego), file = paste0("ResultsCounts2/",contraste,"_ego.xlsx"))
  openxlsx::write.xlsx(data.frame(ego2), file = paste0("ResultsCounts2/",contraste,"_ego2.xlsx"))


  ego3 <- gseGO(geneList     = geneListGSEA_sinNA2,
                OrgDb        = org.Hs.eg.db,
                ont          = "CC",
                nPerm        = 1000,
                minGSSize    = 100,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                verbose      = FALSE)

  ego3 <- setReadable(ego3, OrgDb = org.Hs.eg.db)
  openxlsx::write.xlsx(data.frame(ego3), file = paste0("ResultsCounts2/",contraste,"_ego3.xlsx"))


  # KEGG

  # data(geneList, package="DOSE")
  # gene <- names(geneList)[abs(geneList) > 2]

  kk <- enrichKEGG(gene         = entrezid2,
                   organism     = 'hsa',
                   pvalueCutoff = 0.05)

  kk <- setReadable(ego3, OrgDb = org.Hs.eg.db)
  openxlsx::write.xlsx(data.frame(kk), file = paste0("ResultsCounts2/",contraste,"_KEGG_enrich.xlsx"))

  kk2 <- gseKEGG(geneList     = geneListGSEA_sinNA2,
                 organism     = 'hsa',
                 nPerm        = 1000,
                 minGSSize    = 120,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)
  head(kk2)
  # kk2 <- setReadable(kk2, OrgDb = org.Hs.eg.db)
  openxlsx::write.xlsx(data.frame(kk2), file = paste0("ResultsCounts2/",contraste,"_KEGG_gsea.xlsx"))

  mkk <- enrichMKEGG(gene = entrezid2,
                     organism = 'hsa')
  # mkk <- setReadable(mkk, OrgDb = org.Hs.eg.db)
  openxlsx::write.xlsx(data.frame(mkk), file = paste0("ResultsCounts2/",contraste,"_KEGG_mkk.xlsx"))

  mkk2 <- gseMKEGG(geneList = geneListGSEA_sinNA2,
                   organism = 'hsa')
  # mkk2 <- setReadable(mkk2, OrgDb = org.Hs.eg.db)
  openxlsx::write.xlsx(data.frame(mkk2), file = paste0("ResultsCounts2/",contraste,"_KEGG_mkk2_gsea.xlsx"))

  # REACTOME

  library(GOSemSim)
  hsGO <- godata('org.Hs.eg.db', ont="MF")

  library(AnnotationHub)
  hub <- AnnotationHub()
  q <- query(hub, "Cricetulus")
  id <- q$ah_id[length(q)]
  Cgriseus <- hub[[id]]

  goSim("GO:0004022", "GO:0005515", semData=hsGO, measure="Jiang")
  goSim("GO:0004022", "GO:0005515", semData=hsGO, measure="Wang")

  go1 = c("GO:0004022","GO:0004024","GO:0004174")
  go2 = c("GO:0009055","GO:0005515")
  mgoSim(go1, go2, semData=hsGO, measure="Wang", combine=NULL)

  mgoSim(go1, go2, semData=hsGO, measure="Wang", combine="BMA")

  geneSim("241", "251", semData=hsGO, measure="Wang", combine="BMA")

  hsGO2 <- godata('org.Hs.eg.db', keytype = "ENTREZID", ont="MF", computeIC=FALSE)
  # genes <- c("CDC45", "MCM10", "CDC20", "NMU", "MMP1")

  mgeneSim(entrezid2[1:4], semData=hsGO2, measure="Wang", combine="BMA", verbose=FALSE)


  library(ReactomePA)

  x <- enrichPathway(gene=entrezid2, pvalueCutoff = 0.05, readable=TRUE)
  head(x)

  openxlsx::write.xlsx(data.frame(x), file = paste0("ResultsCounts2/",contraste,"_ReactomePA_hiper.xlsx"))

  y <- gsePathway(geneListGSEA_sinNA2,
                  pvalueCutoff = 0.2,
                  pAdjustMethod = "BH",
                  verbose = FALSE)
  head(y)

  # library(graphite)
  # viewPathway(" The citric acid (TCA) cycle and respiratory electron transport",
  #             readable = TRUE,
  #             foldChange = geneListGSEA_sinNA2)


  library("pathview")
  hsa04080 <- pathview(gene.data  = geneListGSEA_sinNA2,
                       pathway.id = "hsa04080",
                       species    = "hsa",
                       limit      = list(gene=max(abs(geneListGSEA_sinNA2)), cpd=1))

  # ------------------------------------------------------------------------

}


