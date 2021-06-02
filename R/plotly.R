library(plotly)
mtcars %>%
  highlight_key(~cyl) %>%
  plot_ly(
    x = ~wt, y = ~mpg, text = ~cyl, mode = "markers+text",
    textposition = "top", hoverinfo = "x+y"
  ) %>%
  highlight(on = "plotly_hover", off = "plotly_doubleclick")



# ggg ---------------------------------------------------------------------

dot_plot <- base %>%
  summarise(miss = sum(is.na(median))) %>%
  filter(miss > 0) %>%
  add_markers(
    x = ~miss,
    y = ~forcats::fct_reorder(city, miss),
    hoverinfo = "x+y"
  ) %>%
  layout(
    xaxis = list(title = "Number of months missing"),
    yaxis = list(title = "")
  )

subplot(dot_plot, time_series, widths = c(.2, .8), titleX = TRUE) %>%
  layout(showlegend = FALSE) %>%
  highlight(on = "plotly_selected", dynamic = TRUE, selectize = TRUE)

# fff ---------------------------------------------------------------------



dat <- openxlsx::read.xlsx("G:/Mi unidad/NASIR/MAIN_FOLDER_NASIR/MAIN_TABLE_v003_shiny.xlsx", 2)

dat %>%
  highlight_key(~ID) %>%
  plot_ly(
    x = ~BENEFIT2, y = ~Only50All2, text = ~ID, mode = "markers+text",
    textposition = "top", hoverinfo = "text"
  ) %>%
  highlight(on = "plotly_hover", off = "plotly_doubleclick")

