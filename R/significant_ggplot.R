library(ggplot2)
library(ggsignif)



t.test(dd5$ENSG00000133059[dd5$cnv != "High High"] ~ dd5$cnv[dd5$cnv != "High High"])
t.test(dd5$ENSG00000133059[dd5$cnv != "High"] ~ dd5$cnv[dd5$cnv != "High"])
t.test(dd5$ENSG00000133059[dd5$cnv != "Low"] ~ dd5$cnv[dd5$cnv != "Low"])

annotation_df <- data.frame(
  start = c("Low", "High", "Low"),
  end = c("High", "High High","High High"),
  y = c(4800, 5000, 6000),
  label = c("0.0003 ***", "2.628e-06 ***","0.00296 ***")
)

ggplot(dd5) +
  aes(x = "", y = ENSG00000133059, fill = cnv) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Accent") +
  theme_classic() +
  ylab("DYSTK expres norm") +
  geom_signif(data = annotation_df,
              aes(xmin = start, xmax = end, annotations = label, y_position = y),
              textsize = 3, vjust = -0.2,
              manual = TRUE
  )


ggplot(dd5) +
  aes(x = "", y = ENSG00000133059, fill = cnv) +
  geom_boxplot() +
  geom_signif(annotation = "0.00296 ***",
              y_position = 6000,
              xmin = 0.75, xmax = 1.25,
              tip_length = c(0.2, 0.004)
  )
