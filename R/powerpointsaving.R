
library(officer)

doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(doc, value = volc1, location = ph_location_fullsize())

doc <- add_slide(doc)
doc <- ph_with(doc, value = volc2, location = ph_location_fullsize())

doc <- add_slide(doc)
doc <- ph_with(doc, value = volc3, location = ph_location_fullsize())

doc <- add_slide(doc)
doc <- ph_with(doc, value = volc4, location = ph_location_fullsize())

doc <- add_slide(doc)
doc <- ph_with(doc, value = plot_pca, location = ph_location_fullsize())

doc <- add_slide(doc)
doc <- ph_with(doc, value = plot_pca2, location = ph_location_fullsize())

print(doc, target = paste0("ResultsCounts/",contraste, "_Volcanos_results.pptx"))
