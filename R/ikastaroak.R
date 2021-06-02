
# Genomic book
http://genomicsclass.github.io/book/
https://bookdown.org/
https://awesomeopensource.com/  #Proiektu eta link pila daude, rankeatuak
https://github.com/const-ae/ggsignif   # significant ggplots
https://pjbartlein.github.io/GeogDataAnalysis/index.html # Geographical very complete course
https://pjbartlein.github.io/REarthSysSci/index.html # aurrekoan bertsio berria
https://towardsdatascience.com/the-most-underrated-r-packages-254e4a6516a1 # package hauek bukatu begiratzen

# Links, tools and resources to check
https://github.com/crazyhottommy/getting-started-with-genomics-tools-and-resources
https://m-clark.github.io/

# Survival
https://cran.r-project.org/web/packages/finalfit/vignettes/survival.html


# Github users
https://github.com/crazyhottommy
https://github.com/jokergoo # complex plots
https://github.com/hadley  # ggplot, tidyverse

# Websites
https://www.r-graph-gallery.com/

# plotly
https://plotly-r.com/client-side-linking.html


# addins
install.packages("ggThemeAssist") # ggplot automatiko, pixkat zailagoa
install.packages("esquisse") # ggplot hobeak
install.packages("questionr") # recoding factors, numeric...
install.packages("remedy") # rmkardown zuzenean
remotes::install_github("ThinkR-open/remedy") # rmakrdown formatoak - OSO ONA
install.packages("snakecase") # convierte texto en variable
install.packages("reprex") # prepara el study case para preguntar en foros
install.packages('addinslist')  # lista todos los addins, ondo begiratzeko
https://github.com/daattali/addinslist
install.packages('Datapasta')  # datuak pegatzeko oso interesagarria
install.packages('ViewPipeSteps')  # %>% pausu bakoitza banatzen du pantailan edo taulatan
devtools::install_github('mwip/beautifyR') # Itxura polita emateko textu eta taulei, sin mas
devtools::install_github("milesmcbain/mufflr") # %>% automatikoki jartzen ditu
remotes::install_github("dreamRs/viewxl") #ikusteko dataset bat excelen
devtools::install_github("seasmith/AlignAssign") # Formateatzen laguntzen du nahiko ondo
install.packages("colourpicker") # koloreekin lan egiteko
devtools::install_github("dracodoc/mischelper") # funtzio ugari ditu denborak neurtzeko etabar
remotes::install_github("gadenbuie/regexplain") # regular expression tutorial





a   <- 1:5
bbb <- 6:10
c   <- letters

# This is a commented line
# So is this
a      <- 1:5
b      <- 6:10
copy_a <- a
# More comments

my_df <- head(iris)

plot()

dat <- iris
plot(dat$Sepal.Length)

## Recoding dat$Species into dat$Species_rec
dat$Species_rec <- fct_recode(dat$Species,
  "seto" = "setosa",
  "versi" = "versicolor",
  "virgi" = "virginica"
)

snakecase::to_random_case("Holasdfasf asdóasdf")
snakecase::to_lower_camel_case("Holasdfasf asdóasdf")
snakecase::to_lower_upper_case("Holasdfasf asdóasdf")
snakecase::to_snake_case("Holasdfasf asdóasdf")
snakecase::to_upper_camel_case("Holasdfasf asdóasdf")
snakecase::to_title_case("Holasdfasf asdóasdf")
snakecase::to_mixed_case("Holasdfasf asdóasdf")

library(ggplot2)
library(tidyverse)
diamonds %>%
  select(carat, cut, color, clarity, price) %>%
  group_by(color) %>%
  summarise(n = n(), price = mean(price)) %>%
  arrange(desc(color)) %>%
  print_pipe_steps() -> result

