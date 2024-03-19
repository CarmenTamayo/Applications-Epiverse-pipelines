
#This script uses simulated outbreak linelist data, that is first aggregated
#into weekly cases and deaths, to showcase how {EpiEstim} can be used to reconstruct
#daily incidence and deaths from aggregated data. These reconstructed estimates
#are then used as an input for {cfr} to calculate static and rolling CFR. 
#Resulting cfr estimates are then compared to those obtained when using linelist data directly,
#as well as with those obtained from assigning cases and deaths to first day of the week

###### Simulating data manually ######

library(epiparameter)
library(bpmodels)
library(incidence2)
library(tidyverse)

# Define parameters

serial_dist <- epiparameter::epidist_db(disease ="COVID-19", epi_dist = "serial_interval", single_epidist = T)
serial_func <- function(x){rlnorm(x,meanlog=2.6,sdlog=0.8)}
onset_to_hosp <- epiparameter::epidist_db(disease = "COVID-19", epi_dist = "onset_to_hospitalisation", single_epidist = T)
onset_to_death <- epiparameter::epidist_db(disease = "COVID-19", epi_dist = "onset_to_death", single_epidist = T)

set.seed(253456276)
out <- chain_sim(n=5, "pois", "size", lambda=1.2, serial = serial_func, tree=T, infinite = 1000)
out <- out %>% filter(n == 2)
names(out)[3]='infector'


#Adding dates with delay distributions#

set_calendar_date <- as.Date("2023-01-01")
out$time_rounded <- round(out$time, digits = 0)
out$onset_date <-  out$time_rounded + set_calendar_date

# Add deaths

death_params <- get_parameters(onset_to_death)
out$death <- NA
out$death <- out$time + rlnorm(nrow(out),meanlog =death_params[1],sdlog = death_params[2])
out$death_rounded <- round(out$death, digits = 0)

#Making it so that not everyone dies

dead_people <- sample(out$death_rounded, replace = F, size = 0.38*length(out$death_rounded))
out$death_rounded[dead_people] <- NA
out$death_date <- out$death_rounded + set_calendar_date

#Final linelist

linelist_sim <- out[, c(7,10)]
linelist_sim <- linelist_sim[c(1:522),]

#Aggregating data- when using incidence2's argument by=7, it appears like the data isn't actually being grouped, it just cuts the data off
#every 7 days, which results in missing data from the original linelist?
#Process repeated for cases and deaths

weekly_cases<-linelist_sim %>%
  mutate(week_onset = floor_date(onset_date, "week"),  # Extract the first day of the week
         total_cases = 1) %>%  # Set each row as one case
  group_by(week_onset) %>%
  summarise(total_cases = sum(total_cases))
weekly_cases <- incidence2::incidence(weekly_cases, "week_onset", counts = "total_cases") %>% complete_dates(by = 7)

weekly_deaths <- linelist_sim %>%
  mutate(week_death = floor_date(death_date, "week"),  # Extract the first day of the week
         total_deaths = 1) %>%  # Set each row as one case
  group_by(week_death) %>%
  summarise(total_deaths = sum(total_deaths))
weekly_deaths <- weekly_deaths %>% filter(!is.na(weekly_deaths$week_death))
weekly_deaths <- incidence2::incidence(weekly_deaths, "week_death", counts = "total_deaths") %>% complete_dates(by = 7)

#Getting complete dates to match with the daily incidence
#These data is also useful for comparison of method vs aggregation of cases to first day of the week

full_incidence <- weekly_cases %>% incidence("date_index", counts = "count") %>% complete_dates()
complete_dates_inc <- full_incidence$date_index
full_deaths <- weekly_deaths %>% incidence("date_index", counts = "count") %>% complete_dates()
complete_dates_death <- full_deaths$date_index

#Getting si for COVID and formatting it to work with EpiEstim
COVID_si <- epidist_db(disease = "COVID", epi_dist = "serial_interval", subset = is_parameterised, single_epidist = T)
COVID_si <- discretise(COVID_si)
wrap_si <- function(si) {
  domain <- seq(1L, to = si$prob_dist$qf(0.999), by = 1L)
  pmf <- si$prob_dist$d(domain)
  pmf[1] <- 0
  pmf <- pmf / sum(pmf)
  pmf
}
si_distr <- wrap_si(COVID_si)


#Using estimate_R with dt=7 to reconstruct daily incidence
library(EpiEstim)
Rt_est_incidence <- estimate_R(incid = weekly_cases$count,
                     dt = 7L,
                     dt_out = 7L,
                     recon_opt = "naive",
                     iter = 10L,
                     tol = 1e-6,
                     grid = list(precision = 0.001, min = -1, max = 1),
                     config = make_config(si_distr = si_distr),
                     method = "non_parametric_si")


I_inc <- as.data.frame(Rt_est_incidence$I[-c(1:6)])
names(I_inc)="cases"

#Now reconstructing also daily deaths

Rt_est_deaths <- estimate_R(incid = weekly_deaths$count,
                               dt = 7L,
                               dt_out = 7L,
                               recon_opt = "naive",
                               iter = 10L,
                               tol = 1e-6,
                               grid = list(precision = 0.001, min = -1, max = 1),
                               config = make_config(si_distr = si_distr),
                               method = "non_parametric_si")


I_deaths <- as.data.frame(Rt_est_deaths$I[-c(1:6)])
names(I_deaths)= "deaths"

#Now using daily incidence to estimate CFR
#First we create the data frame for {cfr} with daily incidence and the corresponding dates

daily_incidence <- bind_cols(complete_dates_inc, I_inc)
names(daily_incidence)[1]= "date"
daily_deaths <- bind_cols(complete_dates_death, I_deaths)
names(daily_deaths)[1]= "date"

daily_data <- merge(daily_incidence, daily_deaths, by = "date", all.x = T, all.y = T)
daily_data <- replace(daily_data, is.na(daily_data), 0)

plot(daily_data$date,daily_data$cases,type="l", col = "blue")
lines(daily_data$date,daily_data$deaths, type = "l", col="red")

#Now, to be able to use cfr for daily incidence, it's necessary to round the incidence and deaths to an integer

daily_data[,2:3] <- round(daily_data[,2:3])

#Getting onset to death delay for COVID

onset_death_covid <- epidist_db(disease = "COVID", epi_dist = "onset to death", single_epidist = T)

#Estimating static cfr
cfr_overall <- cfr_static(daily_data, delay_density = function(x) density(onset_death_covid, x))

#Rolling cfr

cfr_rolling <- cfr_rolling(daily_data, delay_density = function(x) density(onset_death_covid, x))

ggplot(cfr_rolling) +
  geom_ribbon(
    aes(x = date, ymin = severity_low, ymax = severity_high),
    alpha = 0.5, fill = "deepskyblue3") + 
  geom_line(
    aes(x = date, y = severity_mean), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()


#Now comparing to original linelist

linelist_inc <- incidence(linelist_sim, c("onset_date", "death_date")) %>% complete_dates()
data_for_cfr <- prepare_data(linelist_inc, cases_variable = "onset_date", deaths_variable = "death_date")

cfr_static(data_for_cfr, delay_density = function(x) density(onset_death_covid, x))

cfr_rolling_linelist <- cfr_rolling(data_for_cfr, delay_density = function(x) density(onset_death_covid, x))

ggplot(cfr_rolling_linelist) +
  geom_ribbon(
    aes(x = date, ymin = severity_low, ymax = severity_high),
    alpha = 0.5, fill = "deepskyblue3") + 
  geom_line(
    aes(x = date, y = severity_mean), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

##### Comparison to adding cases and deaths to first day of the week #####

full_incidence <- full_incidence[,-2]
names(full_incidence)[c(1:2)] <- c("date","cases")

full_deaths <- full_deaths[,-2]
names(full_deaths)[c(1:2)] <- c("date","deaths")

from_weekly_inc <- merge(full_incidence, full_deaths, all.x = T, all.y = T, by = "date")
from_weekly_inc <- replace(from_weekly_inc, is.na(from_weekly_inc), 0)

cfr_static(from_weekly_inc, delay_density = function(x) density(onset_death_covid, x))

cfr_rolling_from_weekly <- cfr_rolling(from_weekly_inc, delay_density = function(x) density(onset_death_covid, x))

ggplot(cfr_rolling_from_weekly) +
  geom_ribbon(
    aes(x = date, ymin = severity_low, ymax = severity_high),
    alpha = 0.5, fill = "deepskyblue3") + 
  geom_line(
    aes(x = date, y = severity_mean), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()



###### Using simulist ##### (pending fixing bugs)

#First simulating linelist data

library(simulist)
library(epiparameter)
library(cfr)
library(tidyverse)
library(incidence2)

contact_distribution <- epiparameter::epidist(
  disease = "COVID-19", 
  epi_dist = "contact distribution", 
  prob_distribution = "pois", 
  prob_distribution_params = c(mean = 2)
)

contact_interval <- epiparameter::epidist(
  disease = "COVID-19",
  epi_dist = "contact interval",
  prob_distribution = "gamma",
  prob_distribution_params = c(shape = 1, scale = 1)
)

# get onset to hospital admission from {epiparameter} database
onset_to_hosp <- epiparameter::epidist_db(
  disease = "COVID-19",
  epi_dist = "onset to hospitalisation",
  single_epidist = TRUE
)

onset_to_death <- epiparameter::epidist_db(
  disease = "COVID-19",
  epi_dist = "onset to death",
  single_epidist = TRUE
)

linelist <- sim_linelist(
  contact_distribution = contact_distribution,
  contact_interval = contact_interval,
  prob_infect = 0.5,
  onset_to_hosp = onset_to_hosp,
  onset_to_death = onset_to_death
)

#Now converting to weekly data to demonstrate functionality

Marburg_aggregated <- linelist %>% incidence(c("date_onset","date_death"))

#Onset to death for Marburg
set.seed(1)
extract_param(type = "range", values = c(8, 2, 16), distribution = "gamma", samples = 77) 
onset_to_death <- epidist(
  disease = "Marburg",
  epi_dist = "onset to death",
  prob_distribution = "gamma",
  prob_distribution_params = c(shape = 2.09, scale = 4.51)
)











