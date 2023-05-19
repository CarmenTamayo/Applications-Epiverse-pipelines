library(readr)
library(tidyverse)

#Loading outbreak scenario data from one of the workshop's groups
Linelist_group3 <- read_csv("Data/Linelist_group3.csv")

#Formatting data as daily incidence of cases and deaths
daily_cases_deaths <- Linelist_group3[,c(6,8)]
names(daily_cases_deaths)[1:2]=c("date","deaths")
daily_cases_deaths <- daily_cases_deaths %>% mutate(cases = rep(1, nrow(daily_cases_deaths)), .after = "date")
daily_cases_deaths$deaths <- ifelse(is.na(daily_cases_deaths$deaths)==T, 0, 1)
daily_cases_deaths <- aggregate(cbind(cases, deaths) ~ date, data = daily_cases_deaths, FUN = sum)




