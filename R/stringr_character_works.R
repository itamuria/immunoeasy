# https://journal.r-project.org/archive/2010-2/RJournal_2010-2_Wickham.pdf
# https://stringr.tidyverse.org/articles/regular-expressions.html

library(stringr)

str_length("abc")
x <- c("abcdef", "ghifjk")
str_sub(x, 3, 3)
str_sub(x, 2, -2)
str_sub(x, 3, 3) <- "X"
x
str_dup(x, c(2, 3))


x <- c("abc", "defghi")
str_pad(x, 10) # default pads on left
str_pad(x, 10, "both")

str_pad(x, 4)

x <- c("Short", "This is a long string")

x %>%
  str_trunc(10) %>%
  str_pad(10, "right")

x <- c("  a   ", "b   ",  "   c")
str_trim(x)
str_trim(x, "left")

jabberwocky <- str_c(
  "`Twas brillig, and the slithy toves ",
  "did gyre and gimble in the wabe: ",
  "All mimsy were the borogoves, ",
  "and the mome raths outgrabe. "
)
cat(str_wrap(jabberwocky, width = 40))

x <- "I like horses."
str_to_upper(x)
#> [1] "I LIKE HORSES."
str_to_title(x)
#> [1] "I Like Horses."

str_to_lower(x)
#> [1] "i like horses."
# Turkish has two sorts of i: with and without the dot
str_to_lower(x, "tr")
#> [1] "ı like horses."

x <- c("y", "i", "k")
str_order(x)
#> [1] 2 3 1

str_sort(x)
#> [1] "i" "k" "y"
# In Lithuanian, y comes between i and k
str_sort(x, locale = "lt")
#> [1] "i" "y" "k"

strings <- c(
  "apple",
  "219 733 8965",
  "329-293-8753",
  "Work: 579-499-7527; Home: 543.355.3679"
)
phone <- "([2-9][0-9]{2})[- .]([0-9]{3})[- .]([0-9]{4})"

# Which strings contain phone numbers?
str_detect(strings, phone)
#> [1] FALSE  TRUE  TRUE  TRUE
str_subset(strings, phone)
#> [1] "219 733 8965"
#> [2] "329-293-8753"
#> [3] "Work: 579-499-7527; Home: 543.355.3679"

# How many phone numbers in each string?
str_count(strings, phone)
#> [1] 0 1 1 2

# Where in the string is the phone number located?
(loc <- str_locate(strings, phone))
#>      start end
#> [1,]    NA  NA
#> [2,]     1  12
#> [3,]     1  12
#> [4,]     7  18
str_locate_all(strings, phone)
#> [[1]]
#>      start end
#>
#> [[2]]
#>      start end
#> [1,]     1  12
#>
#> [[3]]
#>      start end
#> [1,]     1  12
#>
#> [[4]]
#>      start end
#> [1,]     7  18
#> [2,]    27  38

# What are the phone numbers?
str_extract(strings, phone)
#> [1] NA             "219 733 8965" "329-293-8753" "579-499-7527"
str_extract_all(strings, phone)
#> [[1]]
#> character(0)
#>
#> [[2]]
#> [1] "219 733 8965"
#>
#> [[3]]
#> [1] "329-293-8753"
#>
#> [[4]]
#> [1] "579-499-7527" "543.355.3679"
str_extract_all(strings, phone, simplify = TRUE)
#>      [,1]           [,2]
#> [1,] ""             ""
#> [2,] "219 733 8965" ""
#> [3,] "329-293-8753" ""
#> [4,] "579-499-7527" "543.355.3679"

# Pull out the three components of the match
str_match(strings, phone)
#>      [,1]           [,2]  [,3]  [,4]
#> [1,] NA             NA    NA    NA
#> [2,] "219 733 8965" "219" "733" "8965"
#> [3,] "329-293-8753" "329" "293" "8753"
#> [4,] "579-499-7527" "579" "499" "7527"
str_match_all(strings, phone)
#> [[1]]
#>      [,1] [,2] [,3] [,4]
#>
#> [[2]]
#>      [,1]           [,2]  [,3]  [,4]
#> [1,] "219 733 8965" "219" "733" "8965"
#>
#> [[3]]
#>      [,1]           [,2]  [,3]  [,4]
#> [1,] "329-293-8753" "329" "293" "8753"
#>
#> [[4]]
#>      [,1]           [,2]  [,3]  [,4]
#> [1,] "579-499-7527" "579" "499" "7527"
#> [2,] "543.355.3679" "543" "355" "3679"

str_replace(strings, phone, "XXX-XXX-XXXX")
#> [1] "apple"
#> [2] "XXX-XXX-XXXX"
#> [3] "XXX-XXX-XXXX"
#> [4] "Work: XXX-XXX-XXXX; Home: 543.355.3679"
str_replace_all(strings, phone, "XXX-XXX-XXXX")
#> [1] "apple"
#> [2] "XXX-XXX-XXXX"
#> [3] "XXX-XXX-XXXX"
#> [4] "Work: XXX-XXX-XXXX; Home: XXX-XXX-XXXX"

str_split("a-b-c", "-")
#> [[1]]
#> [1] "a" "b" "c"
str_split_fixed("a-b-c", "-", n = 2)
#>      [,1] [,2]
#> [1,] "a"  "b-c"

a1 <- "\u00e1"
a2 <- "a\u0301"
c(a1, a2)
#> [1] "á" "á"
a1 == a2
#> [1] FALSE

str_detect(a1, fixed(a2))
#> [1] FALSE
str_detect(a1, coll(a2))
#> [1] TRUE

i <- c("I", "İ", "i", "ı")
i
#> [1] "I" "İ" "i" "ı"

str_subset(i, coll("i", ignore_case = TRUE))
#> [1] "I" "i"
str_subset(i, coll("i", ignore_case = TRUE, locale = "tr"))
#> [1] "İ" "i"

x <- "This is a sentence."
str_split(x, boundary("word"))
#> [[1]]
#> [1] "This"     "is"       "a"        "sentence"
str_count(x, boundary("word"))
#> [1] 4
str_extract_all(x, boundary("word"))
#> [[1]]
#> [1] "This"     "is"       "a"        "sentence"

str_split(x, "")
#> [[1]]
#>  [1] "T" "h" "i" "s" " " "i" "s" " " "a" " " "s" "e" "n" "t" "e" "n" "c"
#> [18] "e" "."
str_count(x, "")
#> [1] 19



# https://journal.r-project.org/archive/2010-2/RJournal_2010-2_Wickham.pdf


library(stringr)
strings <- c(" 219 733 8965", "329-293-8753 ", "banana", "595 794 7569",
             "387 287 6718", "apple", "233.398.9187 ", "482 952 3315", "239 923 8115",
             "842 566 4692", "Work: 579-499-7527", "$1000", "Home: 543.355.3679")
phone <- "([2-9][0-9]{2})[- .]([0-9]{3})[- .]([0-9]{4})"
# Which strings contain phone numbers?
str_detect(strings, phone)
strings[str_detect(strings, phone)]
# Where in the string is the phone number located?
loc <- str_locate(strings, phone)
loc
# Extract just the phone numbers
str_sub(strings, loc[, "start"], loc[, "end"])
# Or more conveniently:
str_extract(strings, phone)
# Pull out the three components of the match
str_match(strings, phone)
# Anonymise the data
str_replace(strings, phone, "XXX-XXX-XXXX")

library(stringr)
col2hex <- function(col) {
  rgb <- col2rgb(col)
  rgb(rgb["red", ], rgb["green", ], rgb["blue", ], max = 255)
}
# Goal replace colour names in a string with their hex equivalent
strings <- c("Roses are red, violets are blue", "My favourite colour is green")
colours <- str_c("\\b", colors(), "\\b", collapse="|")
# This gets us the colours, but we have no way of replacing them
str_extract_all(strings, colours)
# Instead, let's work with locations
locs <- str_locate_all(strings, colours)
sapply(seq_along(strings), function(i) {
  string <- strings[i]
  loc <- locs[[i]]
  # Convert colours to hex and replace
  hex <- col2hex(str_sub(string, loc[, "start"], loc[, "end"]))
  str_sub(string, loc[, "start"], loc[, "end"]) <- hex
  string
})
