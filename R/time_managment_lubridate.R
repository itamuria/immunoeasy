# https://www.r-statistics.com/2012/03/do-more-with-dates-and-times-in-r-with-lubridate-1-1-0/

library(lubridate)
ymd("20110604"); mdy("06-04-2011"); dmy("04/06/2011")

arrive <- ymd_hms("2011-06-04 12:00:00", tz = "Pacific/Auckland")
leave <- ymd_hms("2011-08-10 14:00:00", tz = "Pacific/Auckland")

second(arrive)
## 0
second(arrive) <- 25
arrive
## "2011-06-04 12:00:25 NZST"
second(arrive) <- 0
wday(arrive)
## 7
wday(arrive, label = TRUE)
## Sat

meeting <- ymd_hms("2011-07-01 09:00:00", tz = "Pacific/Auckland")
## "2011-07-01 09:00:00 NZST"
with_tz(meeting, "America/Chicago")
## "2011-06-30 16:00:00 CDT"

mistake <- force_tz(meeting, "America/Chicago")
## "2011-07-01 09:00:00 CDT"
with_tz(mistake, "Pacific/Auckland")
## "2011-07-02 02:00:00 NZST"












